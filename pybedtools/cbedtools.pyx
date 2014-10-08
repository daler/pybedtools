# cython: profile=True

"""
    bedtools.pyx: A Cython wrapper for the BEDTools BedFile class

    Authors: Aaron Quinlan[1], Brent Pedersen[2]
    Affl:    [1] Center for Public Health Genomics, University of Virginia
             [2]
    Email:  aaronquinlan at gmail dot com
"""
from cython.operator cimport dereference as deref
import sys
import subprocess
from collections import defaultdict

cdef dict LOOKUPS = {
    "gff":  {"chrom": 0, "start": 3, "end": 4, "stop": 4, "strand": 6},
    "vcf":  {"chrom": 0, "start": 1},
    "bed":  {"chrom": 0, "start": 1, "end": 2, "stop": 2, "score": 4, "strand": 5}
}
for ktype, kdict in LOOKUPS.items():
    for k, v in kdict.items():
        kdict[v] = k

# Keys are tuples of start/start, stop/stop, start/stop, stop/start.
# Values are which operators should return True, otherwise False
# < 0 | <= 1 | == 2 | != 3 |  > 4 | >= 5
PROFILES_TRUE = {
                (0, 0, -1, 1): (2, 1, 5),  # a == b, a >= b, a <= b
                # a  ---------
                # b  ---------

                (-1, -1, -1, -1): (0, 1),  # a < b, a <= b
                # a ----
                # b       -----

                (-1, -1, -1, 0): (1,),  # a <= b
                # a ----
                # b     -----  (book-ended)

                (1, 1, 0, 1): (5,),  # a >= b
                # a     -----
                # b ----      (book-ended)

                (1, 1, 1, 1): (4, 5), # a > b, a >= b
                # a       ------
                # b ----

                (0, 1, -1, 1): (5,),  # a >= b
                # a  ------------
                # b  ---------

                (1, 0, -1, 1): (5,),  # a >= b
                # a   -----------
                # b -------------

                (-1, 0, -1, 1): (1,),  # a <= b
                # a -------------
                # b   -----------

                (0, -1, -1, 1): (1,), # a <= b
                # a  ---------
                # b  ------------

                (-1, -1, -1, 1): (1,), # a <= b
                # a -----------
                # b        -----------

                (1, 1, -1, 1): (5,),  # a >= b
                # a        -----------
                # b -----------

                (1, -1, -1, 1): tuple(), # undef
                # a    ----
                # b -----------

                (-1, 1, -1, 1): tuple(), # undef
                # a -----------
                # b    ----

                (-1, 0, -1, 0): (1,),  # a <= b
                # a -----------
                # b           -

                (1, 0, 0, 1): (5,),  # a >= b
                # a           -
                # b -----------

                (0, 0, 0, 0): (1, 2, 5),  # a == b, a <= b, a >= b
                # a -
                # b -  (starts and stops are identical for all features)
            }


class MalformedBedLineError(Exception):
    pass


class BedToolsFileError(Exception):
    pass


cdef class Attributes(dict):
    """
    Class to map between a dict of attrs and fields[8] of a GFF Interval obj.
    """
    cdef str sep, field_sep, _attr_str
    cdef dict _quoted

    def __init__(self, attr_str=""):
        self._attr_str = attr_str

        # in general, GFF files will have either as many '=' as ';'
        # (or ';'-1 if there's no trailing ';')
        n_semi = attr_str.count(';')
        n_eq = attr_str.count('=')
        n_quotes = attr_str.count('"')

        if n_eq > n_semi - 1:
            self.sep, self.field_sep = (';', '=')
        else:
            self.sep, self.field_sep = (';', ' ')

        self._quoted = {}

        # TODO: pathological case . . . detect this as GFF:
        #
        #   class_code=" "
        #
        # and this as GTF:
        #
        #   class_code "="

        # quick exit
        if attr_str == "":
            return

        kvs = map(str.strip, attr_str.strip().split(self.sep))
        for field, value in [kv.split(self.field_sep, 1) for kv in kvs if kv]:
            if value.count('"') == 2:
                self._quoted[field] = True
            self[field] = value.replace('"', '')

    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)

    def __str__(self):
        # stringify all items first
        items = []
        for field, val in dict.iteritems(self):
            try:
                if self._quoted[field]:
                    val = '"' + str(val) + '"'
            except KeyError:
                pass
            items.append((field, val))
        return self.sep.join([self.field_sep.join(kvs) for kvs in items])

cdef class Interval:
    """
    Class to represent a genomic interval.

    Constructor::

        Interval(chrom, start, end, name=".", score=".", strand=".", otherfields=None)

    Class to represent a genomic interval of any format.  Requires at least 3
    args: chrom (string), start (int), end (int).

    `start` is *always* the 0-based start coordinate.  If this Interval is to
    represent a GFF object (which uses a 1-based coordinate system), then
    subtract 1 from the 4th item in the line to get the start position in
    0-based coords for this Interval.  The 1-based GFF coord will still be
    available, albeit as a string, in fields[3].

    `otherfields` is a list of fields that don't fit into the other kwargs, and
    will be stored in the `fields` attribute of the Interval.

    All the items in `otherfields` must be strings for proper conversion to
    C++.

    By convention, for BED files, `otherfields` is everything past the first 6
    items in the line.  This allows an Interval to represent composite features
    (e.g., a GFF line concatenated to the end of a BED line)

    But for other formats (VCF, GFF, SAM), the entire line should be passed in
    as a list for `otherfields` so that we can always check the
    Interval.file_type and extract the fields we want, knowing that they'll be
    in the right order as passed in with `otherfields`.

    Example usage:

        >>> from pybedtools import Interval
        >>> i = Interval("chr1", 22, 44, strand='-')
        >>> i
        Interval(chr1:22-44)

        >>> i.start, i.end, i.strand, i.length
        (22L, 44L, '-', 22L)

    """
    def __init__(self, chrom, start, end, name=".", score=".", strand=".", otherfields=None):
        if otherfields is None:
            otherfields = []
        self._bed = new BED()
        self._bed.chrom = string(chrom)
        self._bed.start = start
        self._bed.end = end
        self._bed.name = string(name)
        self._bed.score = string(score)
        self._bed.strand = string(strand)
        fields = [chrom, str(start), str(end), name, score, strand]
        fields.extend(otherfields)
        self._bed.fields = list_to_vector(fields)
        self._attrs = None

    def __copy__(self):
        return create_interval_from_list(self.fields)

    def __hash__(self):
        return hash("\t".join(self.fields))

    property chrom:
        """ the chromosome of the feature"""
        def __get__(self):
            return self._bed.chrom.c_str()

        def __set__(self, char* chrom):
            self._bed.chrom = string(chrom)
            idx = LOOKUPS[self.file_type]["chrom"]
            self._bed.fields[idx] = string(chrom)

    # < 0 | <= 1 | == 2 | != 3 |  > 4 | >= 5
    def __richcmp__(self, other, int op):

        if (self.chrom != other.chrom) or (self.strand != other.strand):
            if op == 3: return True
            return False

        # check all 4 so that we can handle nesting and partial overlaps.
        profile = (cmp(self.start, other.start),
                   cmp(self.stop, other.stop),
                   cmp(self.start, other.stop),
                   cmp(self.stop, other.start))

        try:
            if PROFILES_TRUE[profile] == tuple():
                raise NotImplementedError('Features are nested -- comparison undefined')

            if op != 3:
                if op in PROFILES_TRUE[profile]:
                    return True
                return False
            else:
                if 2 in PROFILES_TRUE[profile]:
                    return False
                return True
        except KeyError:
            raise ValueError('Currently unsupported comparison -- please '
                             'submit a bug report')

    property start:
        """The 0-based start of the feature."""
        def __get__(self):
            return self._bed.start

        def __set__(self, int start):
            self._bed.start = start
            idx = LOOKUPS[self.file_type]["start"]

            # Non-BED files should have 1-based coords in fields
            if self.file_type != 'bed':
                start += 1

            s = str(start)
            self._bed.fields[idx] = string(s)

    property end:
        """The end of the feature"""
        def __get__(self):
            return self._bed.end

        def __set__(self, int end):
            self._bed.end = end
            e = str(end)
            idx = LOOKUPS[self.file_type]["stop"]
            self._bed.fields[idx] = string(e)

    property stop:
        """ the end of the feature"""
        def __get__(self):
            return self._bed.end

        def __set__(self, int end):
            idx = LOOKUPS[self.file_type]["stop"]
            e = str(end)
            self._bed.fields[idx] = string(e)
            self._bed.end = end

    property strand:
        """ the strand of the feature"""
        def __get__(self):
            return self._bed.strand.c_str()

        def __set__(self, strand):
            idx = LOOKUPS[self.file_type]["strand"]
            self._bed.fields[idx] = string(strand)
            self._bed.strand = string(strand)

    property length:
        """ the length of the feature"""
        def __get__(self):
            return self._bed.end - self._bed.start

    cpdef deparse_attrs(self):

        if self._attrs is None: return

        if self.file_type != "gff":
            raise ValueError('Interval.attrs was not None, but this was a non-GFF Interval')

        cdef char *cstr
        tmp = self._attrs.__str__()
        cstr = tmp
        self._bed.fields[8] = string(cstr)

    property fields:
        def __get__(self):
            self.deparse_attrs()
            return string_vec2list(self._bed.fields)


    property attrs:
        def __get__(self):
            cdef string ftype = self._bed.file_type
            if self._attrs is None:
                if ftype == <char *>"gff":
                    self._attrs = Attributes(self._bed.fields[8].c_str())
                else:
                    self._attrs = Attributes("")
            return self._attrs

        def __set__(self, attrs):
            self._attrs = attrs

    # TODO: make this more robust.
    @property
    def count(self):
        return int(self.fields[-1])

    property name:
        """
        >>> import pybedtools
        >>> vcf = pybedtools.example_bedtool('v.vcf')
        >>> [v.name for v in vcf]
        ['rs6054257', 'chr1:16', 'rs6040355', 'chr1:222', 'microsat1']

        """
        def __get__(self):
            cdef string ftype = self._bed.file_type
            if ftype == <char *>"gff":
                """
                # TODO. allow setting a name_key in the BedTool constructor?
                if self.name_key and self.name_key in attrs:
                    return attrs[self.name_key]
                """
                for key in ("ID", "Name", "gene_name", "transcript_id", \
                            "gene_id", "Parent"):
                    if key in self.attrs:
                        return self.attrs[key]

            elif ftype == <char *>"vcf":
                s = self.fields[2]
                if s in ("", "."):
                    return "%s:%i" % (self.chrom, self.start)
                return s
            elif ftype == <char *>"bed":
                return self._bed.name.c_str()

        def __set__(self, value):
            cdef string ftype = self._bed.file_type
            if ftype == <char *>"gff":
                for key in ("ID", "Name", "gene_name", "transcript_id", \
                            "gene_id", "Parent"):
                    if not key in self.attrs:
                        continue

                    self.attrs[key] = value
                    break

            elif ftype == <char *>"vcf":
                self._bed.fields[2] = string(value)
            else:
                self._bed.name = string(value)
                self._bed.fields[3] = string(value)

    property score:
        def __get__(self):
            return self._bed.score.c_str()

        def __set__(self, value):
            self._bed.score = string(value)
            idx = LOOKUPS[self.file_type]["score"]
            self._bed.fields[idx] = string(value)

    property file_type:
        "bed/vcf/gff"
        def __get__(self):
            return self._bed.file_type.c_str()

        def __set__(self, value):
            self._bed.file_type = string(value)

    # TODO: maybe bed.overlap_start or bed.overlap.start ??
    @property
    def o_start(self):
        return self._bed.o_start

    @property
    def o_end(self):
        return self._bed.o_end

    @property
    def o_amt(self):
        return self._bed.o_end - self._bed.o_start

    def __str__(self):
        """
        Interval objects always print with a newline to mimic a line in a
        BED/GFF/VCF file
        """
        return "\t".join(self.fields) + "\n"

    def __repr__(self):
        return "Interval(%s:%i-%i)" % (self.chrom, self.start, self.end)

    def __dealloc__(self):
        del self._bed

    def __len__(self):
        return self._bed.end - self._bed.start

    def __getitem__(self, object key):
        cdef int i
        cdef string ftype = self._bed.file_type

        self.deparse_attrs()

        if isinstance(key, (int, long)):
            nfields = self._bed.fields.size()
            if key >= nfields:
                raise IndexError('field index out of range')
            elif key < 0:
                key = nfields + key
            return self._bed.fields.at(key).c_str()
        elif isinstance(key, slice):
            indices = key.indices(self._bed.fields.size())
            return [self._bed.fields.at(i).c_str() for i in range(*indices)]

        elif isinstance(key, basestring):
            if ftype == <char *>"gff":
                try:
                    return self.attrs[key]
                except:
                    pass
            return getattr(self, key)

    def __setitem__(self, object key, object value):
        cdef string ft_string
        cdef char* ft_char
        if isinstance(key, (int, long)):
            nfields = self._bed.fields.size()
            if key >= nfields:
                raise IndexError('field index out of range')
            elif key < 0:
                key = nfields + key
            self._bed.fields[key] = string(value)

            ft_string = self._bed.file_type
            ft = <char *>ft_string.c_str()

            if key in LOOKUPS[ft]:
                setattr(self, LOOKUPS[ft][key], value)

        elif isinstance(key, (basestring)):
            setattr(self, key, value)

    cpdef append(self, object value):
        self._bed.fields.push_back(string(value))


cdef Interval create_interval(BED b):
    cdef Interval pyb = Interval.__new__(Interval)
    pyb._bed = new BED(b.chrom, b.start, b.end, b.name,
                       b.score, b.strand, b.fields,
                       b.o_start, b.o_end, b.bedType, b.file_type, b.status)
    pyb._bed.fields = b.fields
    return pyb

cpdef Interval create_interval_from_list(list fields):
    """
    Create an Interval object from a list of strings.

    Constructor::

        create_interval_from_list(fields)

    Given the list of strings, `fields`, automatically detects the format (BED,
    GFF, VCF, SAM) and creates a new Interval object.

    `fields` is a list with an arbitrary number of items (it can be quite long,
    say after a -wao intersection of a BED12 and a GFF), however, the first
    fields must conform to one of the supported formats.  For example, if you
    want the resulting Interval to be considered a GFF feature, then the first
    9 fields must conform to the GFF format.  Similarly, if you want the
    resulting Interval to be considered a BED feature, then the first three
    fields must be chrom, start, stop.

    Example usage:

        >>> # Creates a BED3 feature
        >>> feature = create_interval_from_list(['chr1', '1', '100'])

    """
    cdef Interval pyb = Interval.__new__(Interval)
    orig_fields = fields[:]
    # BED -- though a VCF will be detected as BED if its 2nd field, id, is a
    # digit
    if (fields[1] + fields[2]).isdigit():
        # if it's too short, just add some empty fields.
        if len(fields) < 7:
            fields.extend(["."] * (6 - len(fields)))
            other_fields = []
        else:
            other_fields = fields[6:]

        pyb._bed = new BED(string(fields[0]), int(fields[1]), int(fields[2]), string(fields[3]),
                string(fields[4]), string(fields[5]), list_to_vector(other_fields))
        pyb.file_type = 'bed'

    # VCF
    elif fields[1].isdigit() and not fields[3].isdigit() and len(fields) >= 8:
        pyb._bed = new BED(string(fields[0]), int(fields[1]), int(fields[1]) + 1,
                           string(fields[2]), string(fields[5]), string('.'),
                           list_to_vector(fields))
        pyb.file_type = 'vcf'

    # SAM
    elif ( len(fields) >= 13) and (fields[1] + fields[3]).isdigit():
        strand = '+'
        if int(fields[1]) & 0x10:
            strand = '-'

        # TODO: what should the stop position be?  Here, it's just the start
        # plus the length of the sequence, but perhaps this should eventually
        # do CIGAR string parsing.
        pyb._bed = new BED(string(fields[2]), int(fields[3])-1, int(fields[3]) + len(fields[9]) - 1,
                           string(strand), string(fields[0]), string(fields[1]), list_to_vector(fields))
        pyb.file_type = 'sam'
    # GFF
    elif len(fields) >= 9 and (fields[3] + fields[4]).isdigit():
        pyb._bed = new BED(string(fields[0]), int(fields[3])-1, int(fields[4]), string(fields[2]),
                           string(fields[5]), string(fields[6]), list_to_vector(fields[7:]))
        pyb.file_type = 'gff'
    else:
        raise MalformedBedLineError('Unable to detect format from %s' % fields)
    pyb._bed.fields = list_to_vector(orig_fields)
    return pyb

cdef vector[string] list_to_vector(list li):
    cdef vector[string] s
    cdef int i
    for i in range(len(li)):
        s.push_back(string(li[i]))
    return s

cdef list string_vec2list(vector[string] sv):
    cdef size_t size = sv.size(), i
    return [sv.at(i).c_str() for i in range(size)]

cdef list bed_vec2list(vector[BED] bv):
    cdef size_t size = bv.size(), i
    cdef list l = []
    cdef BED b
    for i in range(size):
        b = bv.at(i)
        l.append(create_interval(b))
    return l


def overlap(int s1, int s2, int e1, int e2):
    return min(e1, e2) - max(s1, s2)


cdef class IntervalIterator:
    cdef object stream
    cdef int _isstring
    def __init__(self, stream):
        self.stream = stream

        # For speed, check int rather than call isinstance().
        # -1 is unset, 0 assumes list/tuple/iterable, and 1 is a string.
        #
        # Also assumes that all items in the iterable `stream` are the same
        # type...this seems like a reasonable assumption.
        self._isstring = -1

    def __iter__(self):
        return self
    def __next__(self):
        while True:
            try:
                line = self.stream.next()
                if self._isstring < 0:
                    self._isstring = int(isinstance(line, basestring))

            # If you only trap StopIteration, for some reason even after
            # raising a new StopIteration it goes back to the top of the
            # while-loop and tries to get the next line again.  This in turn
            # raises a ValueError, which we catch again . . . and again raise
            # another StopIteration.  Not sure why it works, but it does.
            except (StopIteration, ValueError):
                try:
                    self.stream.close()
                except AttributeError:
                    pass
                raise StopIteration
                break

            if self._isstring:
                if line.startswith(('@', '#', 'track', 'browser')):
                    continue
            break

        if self._isstring:
            fields = line.rstrip('\r\n').split('\t')
        else:
            fields = map(str, line)
        return create_interval_from_list(fields)



cdef class IntervalFile:
    cdef BedFile *intervalFile_ptr
    cdef bint _loaded
    cdef bint _open
    cdef str fn
    """
    An IntervalFile provides low-level access to the BEDTools API.

    >>> fn = pybedtools.example_filename('a.bed')
    >>> intervalfile = pybedtools.IntervalFile(fn)

    """
    def __init__(self, intervalFile):
        self.intervalFile_ptr = new BedFile(string(intervalFile))
        self._loaded = 0
        self._open = 0
        self.fn = intervalFile

    def __dealloc__(self):
        del self.intervalFile_ptr

    def __iter__(self):
        return self

    def __next__(self):
        if not self._open:
            result = self.intervalFile_ptr.Open()
            if result == -1:
                raise BedToolsFileError("Error opening file")
            self._open = 1
        cdef BED b = self.intervalFile_ptr.GetNextBed()
        if b.status == BED_VALID:
            return create_interval(b)
        elif b.status == BED_INVALID:
            self.intervalFile_ptr.Close()
            raise StopIteration
        elif b.status == BED_MALFORMED:
            raise MalformedBedLineError("malformed line: %s" % string_vec2list(b.fields))
        else:
            return self.next()

    @property
    def file_type(self):
        if not self.intervalFile_ptr._typeIsKnown:
            try:
                a = iter(self).next()
                file_type = self.intervalFile_ptr.file_type.c_str()
                self.intervalFile_ptr.Close()
                return file_type
            except MalformedBedLineError:
                # If it's a SAM, raise a meaningful exception.  If not, fail.
                with open(self.fn) as fn:
                    interval = create_interval_from_list(fn.readline().strip().split())
                if interval.file_type == 'sam':
                    raise ValueError('IntervalFile objects do not yet natively support SAM. '
                                     'Please convert to BED/GFF/VCF first if you want to '
                                     'use the low-level API of IntervalFile')
                else:
                    raise


    def loadIntoMap(self):
        """
        Prepares file for checking intersections.  Used by other methods like all_hits()
        """
        if self._loaded:
            return
        self.intervalFile_ptr.loadBedFileIntoMap()
        self._loaded = 1

    def rewind(self):
        """
        Jump to the beginning of the file.
        """
        if not self._open:
            self.intervalFile_ptr.Open()
            self._open = 1
        self.intervalFile_ptr.Rewind()

    def seek(self, offset):
        """
        Jump to a specific byte offset in the file
        """
        if not self._open:
            self.intervalFile_ptr.Open()
            self._open = 1
        self.intervalFile_ptr.Seek(offset)


    def all_hits(self, Interval interval, bool same_strand=False, float overlap=0.0):
        """
        :Signature: `IntervalFile.all_hits(interval, same_strand=False, overlap=0.0)`

        Search for the Interval `interval` this file and return **all**
        overlaps as a list.

        `same_strand`, if True, will only consider hits on the same strand as `interval`.

        `overlap` can be used to specify the fraction of overlap between
        `interval` and each feature in the IntervalFile.

        Example usage:

        >>> fn = pybedtools.example_filename('a.bed')

        >>> # create an Interval to query with
        >>> i = pybedtools.Interval('chr1', 1, 10000, strand='+')

        >>> # Create an IntervalFile out of a.bed
        >>> intervalfile = pybedtools.IntervalFile(fn)

        >>> # get stranded hits
        >>> intervalfile.all_hits(i, same_strand=True)
        [Interval(chr1:1-100), Interval(chr1:100-200), Interval(chr1:900-950)]

        """
        cdef vector[BED] vec_b
        self.loadIntoMap()

        if same_strand == False:
            vec_b = self.intervalFile_ptr.FindOverlapsPerBin(deref(interval._bed), overlap)
            try:
                return bed_vec2list(vec_b)
            finally:
                pass
        else:
            vec_b = self.intervalFile_ptr.FindOverlapsPerBin(deref(interval._bed), same_strand, overlap)
            try:
                return bed_vec2list(vec_b)
            finally:
                pass

    # search() is an alias for all_hits
    search = all_hits

    def any_hits(self, Interval interval, bool same_strand=False, float overlap=0.0):
        """
        :Signature: `IntervalFile.any_hits(interval, same_strand=False, overlap=0.0)`

        Return 1 if the Interval `interval` had >=1 hit in this IntervalFile, 0 otherwise.

        `same_strand`, if True, will only consider hits on the same strand as `interval`.

        `overlap` can be used to specify the fraction of overlap between
        `interval` and each feature in the IntervalFile.

        Example usage:

        >>> fn = pybedtools.example_filename('a.bed')

        >>> # create an Interval to query with
        >>> i = pybedtools.Interval('chr1', 1, 10000, strand='+')

        >>> # Create an IntervalFile out of a.bed
        >>> intervalfile = pybedtools.IntervalFile(fn)

        >>> # any stranded hits?
        >>> intervalfile.any_hits(i, same_strand=True)
        1

        """
        found = 0
        self.loadIntoMap()

        if same_strand == False:
            found = self.intervalFile_ptr.FindAnyOverlapsPerBin(deref(interval._bed), overlap)
        else:
            found = self.intervalFile_ptr.FindAnyOverlapsPerBin(deref(interval._bed), same_strand, overlap)

        return found

    def count_hits(self, Interval interval, bool same_strand=False, float overlap=0.0):
        """
        :Signature: `IntervalFile.count_hits(interval, same_strand=False, overlap=0.0)`

        Return the number of overlaps of the Interval `interval` had with this
        IntervalFile.

        `same_strand`, if True, will only consider hits on the same strand as
        `interval`.

        `overlap` can be used to specify the fraction of overlap between
        `interval` and each feature in the IntervalFile.

        Example usage:

        >>> fn = pybedtools.example_filename('a.bed')

        >>> # create an Interval to query with
        >>> i = pybedtools.Interval('chr1', 1, 10000, strand='+')

        >>> # Create an IntervalFile out of a.bed
        >>> intervalfile = pybedtools.IntervalFile(fn)

        >>> # get number of stranded hits
        >>> intervalfile.count_hits(i, same_strand=True)
        3

        """
        self.loadIntoMap()

        if same_strand == False:
            return self.intervalFile_ptr.CountOverlapsPerBin(deref(interval._bed), overlap)
        else:
            return self.intervalFile_ptr.CountOverlapsPerBin(deref(interval._bed), same_strand, overlap)
