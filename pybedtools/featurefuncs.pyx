from cbedtools import Interval

cpdef center(object feature, int width=100):
    """
    Return the *width* bp from the center of a feature
    """
    if len(feature) < width:
        return feature
    cdef int start = feature.start
    cdef int stop = feature.stop
    cdef int center = start + (stop-start)/2
    feature.start = center - width/2
    if feature.start < 1: feature.start = 1
    feature.stop = center + width/2
    return feature



