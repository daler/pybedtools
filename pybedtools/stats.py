import helpers
import pybedtools


def jaccard(x, y, intersect_kwargs=None):
    """
    Returns the naive Jaccard statistic (intersection over union; both in bp)
    for BedTools `x` and `y`.

    Optionally specify `intersect_kwargs` to adjust how the intersection is
    calculated.

    See comments in source for details on multiprocessing support.
    """
    # Note: This function may run many thousands of times, so we need to be
    # really careful about cleaning up after each operation.
    #
    # It seems that IntervalFile leaks file handles when iterating -- say, when
    # we try to do sum(len(i) for i in bedtool). A workaround is to use an
    # IntervalIterator which cleans up when done.  IntervalIterators are used
    # automatically when specifying using stream=True, so we create a streaming
    # BedTool just before computing the total bp.
    #
    if intersect_kwargs is None:
        intersect_kwargs = {}

    # Streaming so we can iterate over features and then clean up
    intersection = x.intersect(y, stream=True, **intersect_kwargs)

    # To get the union with current BEDTools programs, we first cat the files
    # then do a merge.
    #
    # BedTool.cat() is a little more sophisticated, but returns a file-based
    # BedTool. So a streaming cat is implemented here:
    tmp = pybedtools.BedTool._tmp()
    z = open(tmp, 'w')
    z.write(open(x.fn).read())
    z.write(open(y.fn).read())
    z.close()
    catted = pybedtools.BedTool(tmp)

    # Merge now requires a sort.
    sortd = catted.sort()

    # Streaming so we can iterate over features and then clean up
    union = sortd.merge(stream=True)

    # This is the part that leads to too many open file handles when working
    # with large numbers of files at once (10,000+).  But iterating over
    # a streaming BedTool works fine.
    numerator = sum(len(i) for i in intersection)
    denominator = sum(len(i) for i in union)

    # Be obsessive about deleting tempfiles created and filehandles opened.
    helpers.close_or_delete(intersection)
    helpers.close_or_delete(sortd)
    helpers.close_or_delete(union)
    helpers.close_or_delete(catted)

    return float(numerator) / denominator


def random_jaccard(x, y, genome_fn, shuffle_kwargs=None,
                   intersect_kwargs=None):
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
