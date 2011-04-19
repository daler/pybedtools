from cpython cimport bool
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        string()
        string(char *)
        char *c_str()

"""
Create Cython definitions for the Interval API defined in Interval.h
"""
cdef extern from "bedFile.h":
    cdef enum BedLineStatus:
        BED_INVALID = -1
        BED_HEADER = 0
        BED_BLANK = 1
        BED_VALID = 2

    ctypedef unsigned int   CHRPOS
    ctypedef bint BOOL

    cdef cppclass BED:
        string chrom
        CHRPOS start
        CHRPOS end
        string name
        string score
        string strand
        vector[string] otherFields
        CHRPOS o_start  # the start of an overlap with another interval
        CHRPOS o_end    # the end of an overlap with another interval
        unsigned short bedType
        BOOL isGff
        BOOL isVcf
        BedLineStatus status
        vector[string] fields

        # constructors
        BED()
        BED(string chrom, CHRPOS start, CHRPOS end, string name,
             string score, string strand, vector[string] otherFields,
             CHRPOS o_start, CHRPOS o_end,
             unsigned short bedType, BOOL isGff, BOOL isVcf, BedLineStatus status)

        BED(string chrom, CHRPOS start, CHRPOS end)
        BED(string chrom, CHRPOS start, CHRPOS end, string strand)
        BED(string chrom, CHRPOS start, CHRPOS end, string name,
             string score, string strand, vector[string] otherFields)

        # methods
        string reportBed()


    cdef cppclass BedFile:
        BedFile(string)
        void Open()
        void Close()
        BED  GetNextBed()
        void loadBedFileIntoMap()

        ### "all" ###
        # this version doesn't care if the strands match.
        vector[BED] FindOverlapsPerBin(BED bed, float overlapFraction)
        # if forceStrand is true, require that the strands match,
        vector[BED] FindOverlapsPerBin(BED bed, bool forceStrand, float overlapFraction)

        ### "any" ###
        int FindAnyOverlapsPerBin(BED bed, float overlapFraction)
        # if forceStrand is true, require that the strands match,
        int FindAnyOverlapsPerBin(BED bed, bool forceStrand, float overlapFraction)


        ### "count" ###
        int CountOverlapsPerBin(BED bed, float overlapFraction)
        # if forceStrand is true, require that the strands match,
        int CountOverlapsPerBin(BED bed, bool forceStrand, float overlapFraction)
        bint _isVcf
        bint _isGff
        bint _typeIsKnown
