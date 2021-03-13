import sys
import multiprocessing
from . import helpers
import pybedtools


def _parallel_wrap(
    orig_bedtool,
    shuffle_kwargs,
    genome_fn,
    method,
    method_args,
    method_kwargs,
    sort=False,
    shuffle=True,
    reduce_func=None,
):
    """
    Given a BedTool object `orig_bedtool`, call its `method` with `args` and
    `kwargs` and then call `reduce_func` on the results.

    See parallel_apply docstring for details

    """

    # note: be careful about cleaning up tempfiles
    if not shuffle:
        to_use = orig_bedtool
    else:
        shuffled = orig_bedtool.shuffle(g=genome_fn, **shuffle_kwargs)
        if sort:
            to_use = shuffled.sort()
            helpers.close_or_delete(shuffled)
        else:
            to_use = shuffled

    result = getattr(to_use, method)(*method_args, **method_kwargs)

    if shuffle:
        helpers.close_or_delete(to_use)

    if reduce_func:
        reduced = reduce_func(result)
        if isinstance(result, pybedtools.BedTool):
            helpers.close_or_delete(result)
        return reduced
    else:
        return result


def parallel_apply(
    orig_bedtool,
    method,
    genome=None,
    genome_fn=None,
    method_args=None,
    method_kwargs=None,
    shuffle_kwargs=None,
    shuffle=True,
    reduce_func=None,
    processes=1,
    sort=False,
    _orig_pool=None,
    iterations=1000,
    debug=False,
    report_iterations=False,
):
    """
    Call an arbitrary BedTool method many times in parallel.

    An example use-case is to generate a null distribution of intersections,
    and then compare this to the actual intersections.

    **Important:** due to a known file handle leak in BedTool.__len__, it's
    best to simply check the number of lines in the file, as in the below
    function. This works because BEDTools programs strip any non-interval lines
    in the results.

    >>> # set up example BedTools
    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')

    >>> # Method of `a` to call:
    >>> method = 'intersect'

    >>> # Kwargs provided to `a.intersect` each iteration
    >>> method_kwargs = dict(b=b, u=True)

    >>> # Function that will be called on the results of
    >>> # `a.intersect(**method_kwargs)`.
    >>> def reduce_func(x):
    ...     return sum(1 for _ in open(x.fn))

    >>> # Create a small artificial genome for this test (generally you'd
    >>> # use an assembly name, like "hg19"):
    >>> genome = dict(chr1=(0, 1000))

    >>> # Do 10 iterations using 1 process for this test (generally you'd
    >>> # use 1000+ iterations, and as many processes as you have CPUs)
    >>> results = pybedtools.parallel.parallel_apply(a, method, genome=genome,
    ... method_kwargs=method_kwargs, iterations=10, processes=1,
    ... reduce_func=reduce_func, debug=True, report_iterations=True)

    >>> # get results
    >>> print(list(results))
    [1, 0, 1, 2, 4, 2, 2, 1, 2, 4]

    >>> # We can compare this to the actual intersection:
    >>> reduce_func(a.intersect(**method_kwargs))
    3

    Alternatively, we could use the `a.jaccard` method, which already does the
    reduction to a dictionary.  However, the Jaccard method requires the input
    to be sorted.  Here, we specify `sort=True` to sort each shuffled BedTool
    before calling its `jaccard` method.

    >>> from pybedtools.parallel import parallel_apply
    >>> a = pybedtools.example_bedtool('a.bed')
    >>> results = parallel_apply(a, method='jaccard', method_args=(b,),
    ... genome=genome, iterations=3, processes=1, sort=True, debug=True)
    >>> for i in results:
    ...     print(sorted(i.items()))
    [('intersection', 12), ('jaccard', 0.0171184), ('n_intersections', 1), ('union', 701)]
    [('intersection', 0), ('jaccard', 0.0), ('n_intersections', 0), ('union', 527)]
    [('intersection', 73), ('jaccard', 0.137996), ('n_intersections', 1), ('union', 529)]

    Parameters
    ----------
    orig_bedtool : BedTool

    method : str
        The method of `orig_bedtool` to run

    method_args : tuple
        Passed directly to getattr(orig_bedtool, method)()

    method_kwargs : dict
        Passed directly to getattr(orig_bedtool, method)()

    shuffle : bool
        If True, then `orig_bedtool` will be shuffled at each iteration and
        that shuffled version's `method` will be called with `method_args` and
        `method_kwargs`.

    shuffle_kwargs : dict
        If `shuffle` is True, these are passed to `orig_bedtool.shuffle()`.
        You do not need to pass the genome here; that's handled separately by
        the `genome` and `genome_fn` kwargs.

    iterations : int
        Number of iterations to perform

    genome : string or dict
        If string, then assume it is the assembly name (e.g., hg19) and get
        a dictionary of chromsizes for that assembly, then converts to
        a filename.

    genome_fn : str
        Mutually exclusive with `genome`; `genome_fn` must be an existing
        filename with the chromsizes.  Use the `genome` kwarg instead if you'd
        rather supply an assembly or dict.

    reduce_func : callable
        Function or other callable object that accepts, as its only argument,
        the results from `orig_bedtool.method()`.  For example, if you care
        about the number of results, then you can use `reduce_func=len`.

    processes : int
        Number of processes to run.  If `processes=1`, then multiprocessing is
        not used (making it much easier to debug).  This argument is ignored if
        `_orig_pool` is provided.

    sort : bool
        If both `shuffle` and `sort` are True, then the shuffled BedTool will
        then be sorted.  Use this if `method` requires sorted input.

    _orig_pool : multiprocessing.Pool instance
        If provided, uses `_orig_pool` instead of creating one.  In this case,
        `processes` will be ignored.

    debug : bool
        If True, then use the current iteration index as the seed to shuffle.

    report_iterations : bool
        If True, then report the number of iterations to stderr.
    """

    shuffle_kwargs = shuffle_kwargs or {}
    method_args = method_args or ()
    if not isinstance(method_args, list) and not isinstance(method_args, tuple):
        raise ValueError(
            "method_args must be a list or tuple, got %s" % type(method_args)
        )
    method_kwargs = method_kwargs or {}

    if genome_fn and genome:
        raise ValueError("only of of genome_fn or genome should be provided")

    if shuffle:
        if not genome_fn:
            if not genome:
                raise ValueError(
                    "shuffle=True, so either genome_fn" " or genome must be provided"
                )
            genome_fn = pybedtools.chromsizes_to_file(genome)

    _parallel_wrap_kwargs = dict(
        orig_bedtool=orig_bedtool,
        shuffle_kwargs=shuffle_kwargs,
        genome_fn=genome_fn,
        method=method,
        method_args=method_args,
        method_kwargs=method_kwargs,
        shuffle=shuffle,
        reduce_func=reduce_func,
        sort=sort,
    )

    def add_seed(i, kwargs):
        if debug and shuffle:
            kwargs["shuffle_kwargs"]["seed"] = i
        return kwargs

    if processes == 1:
        for it in range(iterations):
            yield _parallel_wrap(**add_seed(it, _parallel_wrap_kwargs))
        return

    if _orig_pool:
        p = _orig_pool
    else:
        p = multiprocessing.Pool(processes)

    results = [
        p.apply_async(_parallel_wrap, (), add_seed(it, _parallel_wrap_kwargs))
        for it in range(iterations)
    ]
    for i, r in enumerate(results):
        yield r.get()
        if report_iterations:
            sys.stderr.write("%s\r" % i)
            sys.stderr.flush()
