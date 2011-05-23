from cbedtools import Interval

cpdef center(object feature, int width=100):
    """
    Return the *width* bp from the center of a feature.  If a feature is
    smaller than *width*, then return the entire feature.
    """
    if len(feature) < width:
        return feature
    cdef int start = feature.start
    cdef int stop = feature.stop
    cdef int center = start + (stop - start) / 2
    halfwidth = width / 2
    feature.start = center - halfwidth
    if feature.start < 1:
        feature.start = 1
    if halfwidth == 0:
        halfwidth = 1
    feature.stop = center + halfwidth
    return feature


cpdef greater_than(object feature, int size=100):
    """
    Return True if feature length > *size*
    """
    return len(feature) > size


cpdef less_than(object feature, int size=100):
    """
    Return True if feature length < *size*
    """
    return len(feature) < size


cpdef normalized_to_length(object feature, int idx=4, float scalar=0.001):
    """
    Normalizes the value at feature[idx] to the feature's length, in kb.

    *idx*, by default, is the score field for a BED file, but specify any
    integer.

    The value at *idx* will be replaced with its scaled value.

    *scalar* will be multiplied by the value at *idx*, by default this is
    0.001, or per kb.

    Useful for calculating RPKM after running intersect with counts
    """
    feature[idx] = str(float(feature[idx]) * scalar / len(feature))
    return feature


cpdef rename(object feature, str name):
    """
    Forces a rename of all features, e.g., for renaming everything in a file
    'exon'
    """
    feature.name = name
    return feature
