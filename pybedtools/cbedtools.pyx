# distutils: language = c++

# String notes:
#
#   Anything that goes in C++ objects should be converted to a C++ <string>
#   type, using the _cppstr() function.  For example: Interval._bed.file_type,
#   or the entries in Interval._bed.fields.
#
#   Any Python accessor methods (Interval.fields, Interval.__getitem__) should
#   then be converted to Python strings using the _pystr() function.
#
#   Cython uses the `str` type as whatever the native Python version uses as
#   str.


from cpython.version cimport PY_MAJOR_VERSION
from libcpp.string cimport string

# Python byte strings automatically coerce to/from C++ strings.

cdef _cppstr(s):
    # Use this to handle incoming strings from Python.
    #
    # C++ uses bytestrings. PY2 strings need no conversion; bare PY3 strings
    # are unicode and so must be encoded to bytestring.
    if isinstance(s, int):
        s = str(s)
    if isinstance(s, unicode):
        s = s.encode('UTF-8')
    return <string> s

cdef _pystr(string s):
    # Use this to prepare a string for sending to Python.
    #
    # Always returns unicode.
    return s.decode('UTF-8', 'strict')

if PY_MAJOR_VERSION < 3:
    integer_types = (int, long)
else:
    integer_types = (int,)

"""
    bedtools.pyx: A Cython wrapper for the BEDTools BedFile class

    Authors: Aaron Quinlan[1], Brent Pedersen[2]
    Affl:    [1] Center for Public Health Genomics, University of Virginia
             [2]
    Email:  aaronquinlan at gmail dot com
"""
from cython.operator cimport dereference as deref
import sys
import six
import subprocess
from collections import defaultdict

cdef dict LOOKUPS = {
    "gff":  {"chrom": 0, "start": 3, "end": 4, "stop": 4, "strand": 6},
    "vcf":  {"chrom": 0, "start": 1},
    "bed":  {"chrom": 0, "start": 1, "end": 2, "stop": 2, "score": 4, "strand": 5}
}
for ktype, kdict in list(LOOKUPS.items()):
    for k, v in list(kdict.items()):
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


class Attributes(dict):
    """
    Class to map between a dict of attrs and fields[8] of a GFF Interval obj.
    """

    def __init__(self, attr_str=""):
        attr_str = str(attr_str)
        self._attr_str = attr_str
        self.sort_keys = False

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

        pairs = []
        if self.sort_keys:
            items.sort()
        for k, v in items:
            pairs.append(self.field_sep.join([k, v]))

        return self.sep.join(pairs) + self.sep

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


    """
    def __init__(self, chrom, start, end, name=".", score=".", strand=".", otherfields=None):
        if otherfields is None:
            otherfields = []
        otherfields = [_cppstr(i) for i in otherfields]
        self._bed = new BED(
            _cppstr(chrom), start, end, _cppstr(name), _cppstr(score),
            _cppstr(strand), otherfields)

        #self._bed.chrom = _cppstr(chrom)
        #self._bed.start = start
        #self._bed.end = end
        #self._bed.name = _cppstr(name)
        #self._bed.score = _cppstr(score)
        #self._bed.strand = _cppstr(strand)
        fields = [_cppstr(chrom), _cppstr(str(start)), _cppstr(str(end)), _cppstr(name), _cppstr(score), _cppstr(strand)]
        fields.extend(otherfields)
        self._bed.fields = fields
        self._attrs = None

    def __copy__(self):
        return create_interval_from_list(self.fields)

    def __hash__(self):
        return hash("\t".join(self.fields))

    property chrom:
        """ the chromosome of the feature"""
        def __get__(self):
            return _pystr(self._bed.chrom)

        def __set__(self, chrom):
            chrom = _cppstr(chrom)
            self._bed.chrom = chrom
            idx = LOOKUPS[self.file_type]["chrom"]
            self._bed.fields[idx] = _cppstr(chrom)

    # < 0 | <= 1 | == 2 | != 3 |  > 4 | >= 5
    def __richcmp__(self, other, int op):
        if (self.chrom != other.chrom) or (self.strand != other.strand):
            if op == 3: return True
            return False

        def cmp(x, y):
            if x < y:
                return -1
            if x == y:
                return 0
            if x > y:
                return 1


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
            self._bed.fields[idx] = _cppstr(str(start))

    property end:
        """The end of the feature"""
        def __get__(self):
            return self._bed.end

        def __set__(self, int end):
            self._bed.end = end
            idx = LOOKUPS[self.file_type]["stop"]
            self._bed.fields[idx] = _cppstr(str(end))

    property stop:
        """ the end of the feature"""
        def __get__(self):
            return self._bed.end

        def __set__(self, int end):
            idx = LOOKUPS[self.file_type]["stop"]
            self._bed.fields[idx] = _cppstr(str(end))
            self._bed.end = end

    property strand:
        """ the strand of the feature"""
        def __get__(self):
            return _pystr(self._bed.strand)

        def __set__(self, strand):
            idx = LOOKUPS[self.file_type]["strand"]
            self._bed.fields[idx] = _cppstr(strand)
            self._bed.strand = _cppstr(strand)

    property length:
        """ the length of the feature"""
        def __get__(self):
            return self._bed.end - self._bed.start

    cpdef deparse_attrs(self):

        if not self._attrs: return

        if self.file_type != "gff":
            raise ValueError('Interval.attrs was not None, but this was a non-GFF Interval')

        s = self._attrs.__str__()
        self._bed.fields[8] = _cppstr(s)

    property fields:
        def __get__(self):
            self.deparse_attrs()
            items = []
            for i in self._bed.fields:
                if isinstance(i, int):
                    items.append(i)
                else:
                    items.append(_pystr(i))
            return items

    property attrs:
        def __get__(self):
            if self._attrs is None:
                ft = _pystr(self._bed.file_type)
                if ft == 'gff':
                    self._attrs = Attributes(_pystr(self._bed.fields[8]))
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
            value = None
            if ftype == <string>"gff":
                """
                # TODO. allow setting a name_key in the BedTool constructor?
                if self.name_key and self.name_key in attrs:
                    return attrs[self.name_key]
                """
                for key in ("ID", "Name", "gene_name", "transcript_id", \
                            "gene_id", "Parent"):
                    if key in self.attrs:
                        value = self.attrs[key]
                        break

            elif ftype == <string>"vcf":
                s = self.fields[2]
                if s in ("", "."):
                    value = "%s:%i" % (self.chrom, self.start)
                else:
                    value = _pystr(s)
            elif ftype == <string>"bed":
                value = _pystr(self._bed.name)

            return value

        def __set__(self, value):
            cdef string ftype = self._bed.file_type

            if ftype == <string>"gff":
                for key in ("ID", "Name", "gene_name", "transcript_id", \
                            "gene_id", "Parent"):
                    if not key in self.attrs:
                        continue

                    # If it's incoming from Python it's unicode, so store that directly
                    # in the attributes (since an Attribute object works on
                    # unicode)...
                    self.attrs[key] = value
                    break

            # Otherwise use _cppstr() because we're storing it in _bed.fields.
            elif ftype == <string>"vcf":
                self._bed.fields[2] = _cppstr(value)
            else:
                self._bed.name = _cppstr(value)
                self._bed.fields[3] = _cppstr(value)

    property score:
        def __get__(self):
            return _pystr(self._bed.score)

        def __set__(self, value):
            value = _cppstr(value)
            self._bed.score = value
            idx = LOOKUPS[self.file_type]["score"]
            self._bed.fields[idx] = value

    property file_type:
        "bed/vcf/gff"
        def __get__(self):
            return _pystr(self._bed.file_type)

        def __set__(self, value):
            self._bed.file_type = _cppstr(value)

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
        items = []
        for i in self.fields:
            if isinstance(i, int):
                i = str(i)
            items.append(i)

        return '\t'.join(items) + '\n'

    def __repr__(self):
        return "Interval(%s:%i-%i)" % (self.chrom, self.start, self.end)

    def __dealloc__(self):
        del self._bed

    def __len__(self):
        return self._bed.end - self._bed.start

    def __getitem__(self, object key):
        cdef int i
        ftype = _pystr(self._bed.file_type)

        self.deparse_attrs()

        if isinstance(key, (int, long)):
            nfields = self._bed.fields.size()
            if key >= nfields:
                raise IndexError('field index out of range')
            elif key < 0:
                key = nfields + key
            return _pystr(self._bed.fields.at(key))
        elif isinstance(key, slice):
            indices = key.indices(self._bed.fields.size())
            return [_pystr(self._bed.fields.at(i)) for i in range(*indices)]

        elif isinstance(key, str):
            if ftype == "gff":
                try:
                    return self.attrs[key]
                except KeyError:
                    pass
            # We don't have to convert using _pystr() because the __get__
            # methods do that already.
            return getattr(self, key)

    def __setitem__(self, object key, object value):
        if isinstance(key, (int, long)):
            nfields = self._bed.fields.size()
            if key >= nfields:
                raise IndexError('field index out of range')
            elif key < 0:
                key = nfields + key
            self._bed.fields[key] = _cppstr(value)

            ft = _pystr(self._bed.file_type)
            if key in LOOKUPS[ft]:
                setattr(self, LOOKUPS[ft][key], value)

        elif isinstance(key, (basestring)):
            setattr(self, key, value)

    cpdef append(self, object value):
        self._bed.fields.push_back(_cppstr(value))

    def __nonzero__(self):
        return True


cdef Interval create_interval(BED b):
    cdef Interval pyb = Interval.__new__(Interval)
    pyb._bed = new BED(b.chrom, b.start, b.end, b.name,
                       b.score, b.strand, b.fields,
                       b.o_start, b.o_end, b.bedType, b.file_type, b.status)
    pyb._bed.fields = b.fields
    return pyb

# TODO: optimization: Previously we had (fields[1] + fields[2]).isdigit() when
# checking in create_interval_from_list for filetype heuruistics. Is there
# a performance hit by checking instances?
cdef isdigit(s):
    if isinstance(s, integer_types):
        return True
    return s.isdigit()


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

    # TODO: this function is used a lot, and is doing a bit of work. We should
    # have an optimized version that is directly provided the filetype.

    cdef Interval pyb = Interval.__new__(Interval)
    orig_fields = fields[:]
    # BED -- though a VCF will be detected as BED if its 2nd field, id, is a
    # digit

    # SAM
    if (
        (len(fields) >= 11)
        and isdigit(fields[1])
        and isdigit(fields[3])
        and isdigit(fields[4])
        and (fields[5] not in ['.', '+', '-'])
    ):
        # TODO: what should the stop position be?  Here, it's just the start
        # plus the length of the sequence, but perhaps this should eventually
        # do CIGAR string parsing.
        if int(fields[1]) & 0x04:
            # handle unmapped reads
            chrom = _cppstr("*")
            start = 0
            stop = 0
        else:
            chrom = _cppstr(fields[2])
            start = int(fields[3]) - 1
            stop = int(fields[3]) + len(fields[9]) - 1
        name = _cppstr(fields[0])
        score = _cppstr(fields[1])
        if int(fields[1]) & 0x10:
            strand = _cppstr('-')
        else:
            strand = _cppstr('+')

        # Fields is in SAM format
        fields[3] = str(start + 1)

        pyb._bed = new BED(
            chrom,
            start,
            stop,
            strand,
            name,
            score,
            list_to_vector(fields))
        pyb.file_type = _cppstr('sam')


    elif isdigit(fields[1]) and isdigit(fields[2]):
        # if it's too short, just add some empty fields.
        if len(fields) < 7:
            fields.extend([".".encode('UTF-8')] * (6 - len(fields)))
            other_fields = []
        else:
            other_fields = fields[6:]

        pyb._bed = new BED(
            _cppstr(fields[0]),
            int(fields[1]),
            int(fields[2]),
            _cppstr(fields[3]),
            _cppstr(fields[4]),
            _cppstr(fields[5]),
            list_to_vector(other_fields))
        pyb.file_type = _cppstr('bed')

    # VCF
    elif isdigit(fields[1]) and not isdigit(fields[3]) and len(fields) >= 8:
        pyb._bed = new BED(
            _cppstr(fields[0]),
            int(fields[1]) - 1,
            int(fields[1]),
            _cppstr(fields[2]),
            _cppstr(fields[5]),
            _cppstr('.'),
            list_to_vector(fields))
        pyb.file_type = b'vcf'


    # GFF
    elif len(fields) >= 9 and isdigit(fields[3]) and isdigit(fields[4]):
        pyb._bed = new BED(
            _cppstr(fields[0]),
            int(fields[3])-1, int(fields[4]),
            _cppstr(fields[2]),
            _cppstr(fields[5]),
            _cppstr(fields[6]),
            list_to_vector(fields[7:]))
        pyb.file_type = _cppstr('gff')
    else:
        raise MalformedBedLineError('Unable to detect format from %s' % fields)

    if pyb.start > pyb.end:
        raise MalformedBedLineError("Start is greater than stop")
    pyb._bed.fields = list_to_vector(orig_fields)
    return pyb

cdef vector[string] list_to_vector(list li):
    cdef vector[string] s
    cdef int i
    for i in range(len(li)):
        _s = li[i]
        s.push_back(_cppstr(_s))
    return s

cdef list string_vec2list(vector[string] sv):
    cdef size_t size = sv.size(), i
    return [_pystr(sv.at(i)) for i in range(size)]

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
    cdef int _itemtype
    def __init__(self, stream):
        self.stream = stream

        # For speed, check int rather than call isinstance().
        # -1 is unset, 0 assumes list/tuple/iterable, and 1 is a string.
        #
        # Also assumes that all items in the iterable `stream` are the same
        # type...this seems like a reasonable assumption.
        self._itemtype = -1

    def __dealloc__(self):
        try:
            self.stream.close()
        except AttributeError:
            pass

    def __iter__(self):
        return self

    def __next__(self):
        while True:
            if hasattr(self.stream, 'closed'):
                if self.stream.closed:
                    raise StopIteration
            try:
                line = next(self.stream)
            except StopIteration:
                if hasattr(self.stream, 'close'):
                    self.stream.close()
                raise StopIteration

            if self._itemtype < 0:
                if isinstance(line, Interval):
                    self._itemtype = 2
                elif isinstance(line, basestring):
                    self._itemtype = 1
                else:
                    self._itemtype = 0

            if self._itemtype == 1:
                if line.startswith(('@', '#', 'track', 'browser')) or len(line.strip()) == 0:
                    continue
            break

        # Iterable of Interval objects
        if self._itemtype == 2:
            return line

        # Iterable of strings, in which case we need to split
        elif self._itemtype == 1:
            fields = line.rstrip('\r\n').split('\t')

        # Otherwise assume list/tuple/iterable of fields
        else:
            fields = list(line)

        # TODO: optimization: create_interval_from_list should have a version
        # that accepts C++ string instances
        return create_interval_from_list(fields)



cdef class IntervalFile:
    cdef BedFile *intervalFile_ptr
    cdef bint _loaded
    cdef bint _open
    cdef string _fn
    """
    An IntervalFile provides low-level access to the BEDTools API.

    >>> fn = pybedtools.example_filename('a.bed')
    >>> intervalfile = pybedtools.IntervalFile(fn)

    """
    def __init__(self, intervalFile):
        self.intervalFile_ptr = new BedFile(_cppstr(intervalFile))
        self._loaded = 0
        self._open = 0
        self._fn = _cppstr(intervalFile)

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
            self.intervalFile_ptr.Close()
            raise MalformedBedLineError("malformed line: %s" % string_vec2list(b.fields))
        else:
            return next(self)

    @property
    def fn(self):
        return _pystr(self._fn)

    @property
    def file_type(self):
        if not self.intervalFile_ptr._typeIsKnown:
            try:
                a = six.advance_iterator(iter(self))
                file_type = _pystr(self.intervalFile_ptr.file_type)
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
