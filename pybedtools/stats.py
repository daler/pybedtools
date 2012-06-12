import helpers
import pybedtools

def jaccard(x, y, intersect_kwargs=None):
    """
    Returns the naive Jaccard statistic (intersection over union; both in bp)
    for BedTools `x` and `y`.

    Optionally specify `intersect_kwargs` to adjust how the intersection is
    calculated.
    """
    if intersect_kwargs is None:
        intersect_kwargs = {}
    intersection = x.intersect(y, **intersect_kwargs)
    union = x.cat(y, postmerge=True)
    numerator = sum(len(i) for i in intersection)
    denominator = sum(len(i) for i in union)

    # Delete tempfiles immediately to make this function multiprocessing-safe
    helpers.close_or_delete(intersection)
    helpers.close_or_delete(union)

    return float(numerator) / denominator

def random_jaccard(x, y, genome_fn, shuffle_kwargs=None, intersect_kwargs=None):
    """
    Shuffles `x`, using the chromsizes in `genome_fn` and any additional
    `shuffle_kwargs`, then call `jaccard()` with any additional
    `intersect_kwargs`.
    """
    if shuffle_kwargs is None:
        shuffle_kwargs = {}
    z = x.shuffle(g=genome_fn, **shuffle_kwargs)
    result = jaccard(z, y, intersect_kwargs)
    helpers.close_or_delete(z)
    return result

def random_intersection(x, y, genome_fn, shuffle_kwargs, intersect_kwargs):
    z = x.shuffle(g=genome_fn, **shuffle_kwargs)
    zz = z.intersect(y, stream=True, **intersect_kwargs)
    result = len(zz)
    helpers.close_or_delete(z)
    helpers.close_or_delete(zz)
    return result
