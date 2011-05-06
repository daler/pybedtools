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

cpdef parse_attributes(str attr_str):
    """
    parse the attribute string from gff or gtf into a dictionary
    # copied from genomicfeatures
    #>>> parse_attributes('ID=thaliana_1_465_805;match=scaffold_801404.1;rname=thaliana_1_465_805') == {'rname': 'thaliana_1_465_805', 'ID': 'thaliana_1_465_805', 'match': 'scaffold_801404.1'}
    #True
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
        return self.sep.join([self.field_sep.join(kvs) for kvs in self._attr_dict.items()])

    def __repr__(self):
        return repr(self._attr_dict)

cdef class Interval:
    """
    >>> from pybedtools import Interval
    >>> i = Interval("chr1", 22, 44, '-')
    >>> i
    Interval(chr1:22-44)

    >>> i.start, i.end, i.strand, i.length
    (22L, 44L, '-', 22L)

    """
    cdef BED *_bed
    cdef object _attrs

    def __init__(self, chrom, start, end, strand=None):
        if strand is None:
            self._bed = new BED(string(chrom), start, end)
        else:
            self._bed = new BED(string(chrom), start, end, string(strand))
        self._attrs = None

    property chrom:
        """ the chromosome of the feature"""
        def __get__(self):
            return self._bed.chrom.c_str()
        def __set__(self, char* chrom):
            self._bed.chrom = string(chrom)
            idx = LOOKUPS[self.file_type]["chrom"]
            self._bed.fields[idx] = string(chrom)

    property start:
        """ the 0-based start of the feature"""
        def __get__(self):
            return self._bed.start
        def __set__(self, int start):
            self._bed.start = start
            idx = LOOKUPS[self.file_type]["start"]

            # Non-BED files should have 1-based coords in fields
            if self.file_type != 'bed': start += 1

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
                    if key in attrs: return attrs[key]

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
                    if not key in attrs: continue
                    attrs[key] = value
                    attr_str = self._bed.fields[8].c_str()
                    field_sep, quote = ("=", "") if "=" in attr_str \
                                                 else (" ", '"')
                    attr_str = ";".join(["%s%s%s%s%s" % \
                         (k, field_sep, quote, v, quote) \
                             for k,v in attrs.iteritems()])

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
            elif key < 0: key = nfields + key
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
        if isinstance(key, (int,long)):
            nfields = self._bed.fields.size()
            if key >= nfields:
                raise IndexError('field index out of range')
            elif key < 0: key = nfields + key
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
    cdef Interval pyb = Interval.__new__(Interval)
    # BED
    if (fields[1] + fields[2]).isdigit():
        # if it's too short, just add some empty fields.
        if len(fields) < 7:
            fields.extend(["."] * (6 - len(fields)))
            fields.extend([[]])
        pyb._bed = new BED(string(fields[0]), int(fields[1]), int(fields[2]), string(fields[3]),
                       string(fields[4]), string(fields[5]), list_to_vector(fields[6]))
    # GFF
    elif len(fields) == 9 and (fields[3] + fields[4]).isdigit():
        pyb._bed = new BED(string(fields[0]), int(fields[3]), int(fields[4]), string(fields[2]),
                           string(fields[5]), string(fields[6]), list_to_vector(fields[7:]))
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
    return min(e1,e2) - max(s1,s2)


cdef class IntervalFile:
    cdef BedFile *intervalFile_ptr
    cdef bint _loaded
    cdef bint _open

    def __init__(self, intervalFile):
        self.intervalFile_ptr = new BedFile(string(intervalFile))
        self._loaded = 0
        self._open   = 0

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
            raise ValueError("malformed line")
        else:
            return self.next()

    @property
    def file_type(self):
        if not self.intervalFile_ptr._typeIsKnown:
            a = iter(self).next()
        return self.intervalFile_ptr.file_type.c_str()

    def loadIntoMap(self):
        if self._loaded: return
        self.intervalFile_ptr.loadBedFileIntoMap()
        self._loaded = 1

    def all_hits(self, Interval interval, bool same_strand=False, float overlap=0.0):
        """
        Search for the "bed" feature in this file and ***return all overlaps***
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
        Search for the "bed" feature in this file and return
        whether (True/False) >= 1 overlaps are found.
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
        Search for the "bed" feature in this file and return the *** count of hits found ***
        """
        self.loadIntoMap()

        if same_strand == False:
           return self.intervalFile_ptr.CountOverlapsPerBin(deref(interval._bed), overlap)
        else:
           return self.intervalFile_ptr.CountOverlapsPerBin(deref(interval._bed), same_strand, overlap)
