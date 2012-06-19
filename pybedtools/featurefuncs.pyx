# cython: profile=True
#
from cbedtools cimport Interval
from cbedtools import create_interval_from_list


cpdef extend_fields(Interval feature, int n):
    """
    Pads the fields of the feature with "." to a total length of `n` fields,
    """
    fields = feature.fields[:]
    while len(fields) < n:
        fields.append('.')
    return create_interval_from_list(fields)


cpdef center(Interval feature, int width=100):
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


cpdef midpoint(Interval feature):
    """
    Specialized version of `center()` that just returns the single-bp midpoint
    """
    start = feature.start + (feature.stop - feature.start) / 2
    stop = start
    feature.start = start
    feature.stop = stop
    return feature


cpdef greater_than(Interval feature, int size=100):
    """
    Return True if feature length > *size*
    """
    return len(feature) > size


cpdef less_than(Interval feature, int size=100):
    """
    Return True if feature length < *size*
    """
    return len(feature) < size


cpdef normalized_to_length(Interval feature, int idx=4, float scalar=0.001):
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


cpdef rename(Interval feature, str name):
    """
    Forces a rename of all features, e.g., for renaming everything in a file
    'exon'
    """
    feature.name = name
    return feature


cpdef bedgraph_scale(Interval feature, float scalar):
    feature[3] = str(float(feature[3]) * scalar)
    return feature


cdef safe_start_stop(start, stop):
    """
    Ensures that feature start/stop coords are non-negative and that start <
    stop.

    If start is negative, reset to zero.

    If start > stop make start and stop equal to the original start.
    """
    if start < 0:
        start = 0
    if start > stop:
        stop = start
    return start, stop


cpdef TSS(Interval feature, int upstream=500, int downstream=500, add_to_name=None):
    """
    Returns the 5'-most coordinate, plus `upstream` and `downstream` bp; adds
    the string `add_to_name` to the feature's name if provided (e.g., "_TSS")
    """
    return five_prime(feature, upstream, downstream, add_to_name)


cpdef five_prime(Interval feature, int upstream=500, int downstream=500, add_to_name=None):
    """
    Returns the 5'-most coordinate, plus `upstream` and `downstream` bp; adds
    the string `add_to_name` to the feature's name if provided (e.g., "_TSS")
    """
    if feature.strand == '-':
        start = feature.stop - downstream
        stop = feature.stop + upstream
    else:
        start = feature.start - upstream
        stop = feature.start + downstream
    if add_to_name:
        try:
            feature.name += add_to_name
        except AttributeError:
            pass
    feature.start, feature.stop = safe_start_stop(start, stop)
    return feature


cpdef three_prime(Interval feature, int upstream=500, int downstream=500, add_to_name=None):
    """
    Returns the 3'-most coordinate, plus `upstream` and `downstream` bp; adds
    the string `add_to_name` to the feature's name if provided (e.g.,
    "_polyA-site")
    """
    if feature.strand == '-':
        start = feature.start - downstream
        stop = feature.start + upstream
    else:
        start = feature.stop - upstream
        stop = feature.stop + downstream
    if add_to_name:
        try:
            feature.name += add_to_name
        except AttributeError:
            pass
    feature.start, feature.stop = safe_start_stop(start, stop)
    return feature
