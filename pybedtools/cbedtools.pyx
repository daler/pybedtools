# cython: profile=True

"""
    bedtools.pyx: A Cython wrapper for the BEDTools BedFile class

    Authors: Aaron Quinlan[1], Brent Pedersen[2]
    Affl:    [1] Center for Public Health Genomics, University of Virginia
             [2]
    Email:  aaronquinlan at gmail dot com
"""
include "cbedtools.pxi"
from cython.operator cimport dereference as deref
import sys
import subprocess


class MalformedBedLineError(Exception):
    pass


cpdef parse_attributes(str attr_str):
    """
    parse the attribute string from gff or gtf into a dictionary
    """
    cdef str sep, field_sep
    cdef dict _attributes = {}
    sep, field_sep = (";", "=") if "=" in attr_str else (";", " ")
    kvs = map(str.strip, attr_str.strip().split(sep))
    for field, value in [kv.split(field_sep) for kv in kvs if kv]:
        _attributes[field] = value.replace('"', '')
    return _attributes


cdef class Attributes:
    """
    Class to map between a dict of attrs and fields[8] of a GFF Interval obj.
    """
    cdef str sep, field_sep, _attr_str
    cdef object _interval_obj
    cdef dict _attr_dict

    def __init__(self, interval_obj, attr_str=""):
        self._attr_str = attr_str
        self._interval_obj = interval_obj
        self._attr_dict = {}

        # quick exit
        if attr_str == "":
            return

        self.sep, self.field_sep = (";", "=") if "=" in attr_str else (";", " ")
        kvs = map(str.strip, attr_str.strip().split(self.sep))
        for field, value in [kv.split(self.field_sep) for kv in kvs if kv]:
            self._attr_dict[field] = value.replace('"', '')

    def __setitem__(self, key, value):
        """
        Sets both the key/item in self.dict *as well as* the interval object's
        attrs field if it's a GFF Interval
        """
        self._attr_dict[key] = value
        if self._interval_obj.file_type == 'gff':
            self._interval_obj[8] = str(self)
        else:
            raise ValueError('Setting attributes not supported for non-GFF-like Intervals')

    def __getitem__(self, key):
        return self._attr_dict[key]

    def __str__(self):
        # stringify all items first
        items = [(i, str(j)) for i,j in self._attr_dict.items()]
        return self.sep.join([self.field_sep.join(kvs) for kvs in items])

    def __repr__(self):
        return repr(self._attr_dict)

cdef class Interval:
    """
    >>> from pybedtools import Interval
    >>> i = Interval("chr1", 22, 44, strand='-')
    >>> i
    Interval(chr1:22-44)

    >>> i.start, i.end, i.strand, i.length
    (22L, 44L, '-', 22L)

    """
    cdef BED *_bed
    cdef object _attrs

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
        if self.chrom != other.chrom:
            if op == 3: return True
            return 0
        n = self.name
        if n is not None and n == other.name:
            v = 0
        else:
            v = cmp(self.start, other.start) or (self.end, other.end)
        if op in (0, 1, 2): return v
        if op in (4, 5): return -v
        return not v

    property start:
        """ the 0-based start of the feature"""
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
        """ the end of the feature"""
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

    property fields:
        def __get__(self):
            return string_vec2list(self._bed.fields)

    property attrs:
        def __get__(self):
            cdef string ftype = self._bed.file_type
            if self._attrs is None:
                if ftype == <char *>"gff":
                    self._attrs = Attributes(self, self._bed.fields[8].c_str())
                else:
                    self._attrs = Attributes(self, "")
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
                attrs = parse_attributes(self._bed.fields[8].c_str())
                for key in ("ID", "Name", "gene_name", "transcript_id", \
                            "gene_id", "Parent"):
                    if key in attrs:
                        return attrs[key]

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
                attrs = parse_attributes(self._bed.fields[8].c_str())
                for key in ("ID", "Name", "gene_name", "transcript_id", \
                            "gene_id", "Parent"):
                    if not key in attrs:
                        continue
                    attrs[key] = value
                    attr_str = self._bed.fields[8].c_str()
                    field_sep, quote = ("=", "") if "=" in attr_str \
                                                 else (" ", '"')
                    attr_str = ";".join(["%s%s%s%s%s" % \
                         (k, field_sep, quote, v, quote) \
                             for k, v in attrs.iteritems()])

                    self._bed.fields[8] = string(attr_str)
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

    @property
    def fields(self):
        return string_vec2list(self._bed.fields)

    def __str__(self):
        return "\t".join(self.fields)

    def __repr__(self):
        return "Interval(%s:%i-%i)" % (self.chrom, self.start, self.end)

    def __dealloc__(self):
        del self._bed

    def __len__(self):
        return self._bed.end - self._bed.start

    def __getitem__(self, object key):
        cdef int i
        cdef string ftype = self._bed.file_type

        if isinstance(key, (int, long)):
            nfields = self._bed.fields.size()
            if key >= nfields:
                raise IndexError('field index out of range')
            elif key < 0:
                key = nfields + key
            return self._bed.fields.at(key).c_str()
        elif isinstance(key, slice):
            return [self._bed.fields.at(i).c_str() for i in \
                    range(key.start or 0,
                          key.stop or self._bed.fields.size(),
                          key.step or 1)]

        elif isinstance(key, basestring):
            if ftype == <char *>"gff":
                attrs = parse_attributes(self._bed.fields[8].c_str())
                try:
                    return attrs[key]
                except:
                    return getattr(self, key)
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
    *fields* is a list with an arbitrary number of items (it can be quite long,
    say after a -wao intersection of a BED12 and a GFF).

    We need to inspect *fields* to make sure that BED class gets the right
    thing.  This means detecting BED or GFF, and looking at how many fields
    there are.  BED constructor gets the first 6 fields; the rest should be a
    list (converted to a vector)
    """
    cdef Interval pyb = Interval.__new__(Interval)
    orig_fields = fields[:]

    # BED
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
    elif ( len(fields) >= 13) and (fields[1] + fields[3]).isdigit():
        strand = '+'
        if int(fields[1]) & 0x10:
            strand = '-'
        pyb._bed = new BED(string(fields[2]), int(fields[3])-1, int(fields[3]) + len(fields[9]) - 1,
                           string(strand), string(fields[0]), string(fields[1]), list_to_vector(fields))
        pyb.file_type = 'sam'
    # GFF
    elif len(fields) >= 9 and (fields[3] + fields[4]).isdigit():
        pyb._bed = new BED(string(fields[0]), int(fields[3])-1, int(fields[4]), string(fields[2]),
                           string(fields[5]), string(fields[6]), list_to_vector(fields[7:]))
        pyb.file_type = 'gff'
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
            self.intervalFile_ptr.Open()
            self._open = 1
        cdef BED b = self.intervalFile_ptr.GetNextBed()
        if b.status == BED_VALID:
            return create_interval(b)
        elif b.status == BED_INVALID:
            raise StopIteration
        elif b.status == BED_MALFORMED:
            raise MalformedBedLineError("malformed line: %s" % string_vec2list(b.fields))
        else:
            return self.next()

    @property
    def file_type(self):
        if not self.intervalFile_ptr._typeIsKnown:
            a = iter(self).next()
        return self.intervalFile_ptr.file_type.c_str()

    def loadIntoMap(self):
        """
        Prepares file for checking intersections.  Used by other methods like all_hits()
        """
        if self._loaded:
            return
        self.intervalFile_ptr.loadBedFileIntoMap()
        self._loaded = 1

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
