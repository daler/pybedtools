.. include:: includeme.rst

Randomization
=============
:mod:`pybedtools` provides some basic functionality for assigning some
significance value to the overlap between two BEDfiles.

The strategy is to randomly shuffle a file many times, each time doing an
intersection with another file of interest and counting the number of
intersections (or computing some other statistic on the overlap).  Upon doing
this many times, an empirical distribution is constructed, and the number of
intersections between the original, un-shuffled file is compared to this
empirical distribution to obtain a p-value, or compared to the median of the
distribution to get a score.

There are two methods, :meth:`pybedtools.BedTool.randomintersection` which does the
brute force randomizations, and :meth:`BedTool.randomstats` which compiles
and reports the results from the former method.

Example workflow
----------------
As a somewhat trivial example, we'll intersect the example `a.bed` with
`b.bed`, taking care to set some options that will let it run in a
determinisitic way so that these tests will run.

We will be shuffling `a.bed`, so we'll need to specify the limits of its
chromosomes with :meth:`BedTool.set_chromsizes`.  Here, we set it to an
artifically small chromosome size so that we can get some meaningful
results in reasonable time.  In practice, you would either supply your own
dictionary or use a string assembly name (e.g., `'hg19'`, `'mm9'`, `'dm3'`,
etc).  The genome-handling code will find the chromsizes we've set, so
there's no need to tell `shuffleBed` which genome file to use each time.

.. doctest::

    >>> chromsizes = {'chr1': (0, 1000)}
    >>> a = pybedtools.example_bedtool('a.bed').set_chromsizes(chromsizes)
    >>> b = pybedtools.example_bedtool('b.bed')

We have the option of specifying what kwargs to provide
:meth:`BedTool.shuffle` and :meth:`BedTool.intersect`, which will be called
each iteration.  In this example, we'll tell `shuffleBed` to only shuffle
within the chromsome just to illustrate the kwargs passing. We also need to
specify how many iterations to perform.  In practice, 1000 or 10000 are
good numbers, but for the sake of this example we'll only do 100.

Last, setting `debug=True` means that the random seed will be set in a
predictable manner so that we'll always get the same results for testing.
In practice, make sure you use `debug=False` (the default) to ensure random
results.

Furthermore, using the `processes` kwarg will substantially speed up the
comparison (e.g., `processes=8` to split the randomizations across 8 cores).

.. doctest::

    >>> results = a.randomintersection(b, iterations=100, shuffle_kwargs={'chrom': True}, debug=True)

`results` is a generator of intersection counts where each number is the
number of times the shuffled `a` intersected with `b`.  We need to convert
it to a list in order to look at it:


.. doctest::

    >>> results = list(results)
    >>> len(results)
    100

    >>> print(results[:10])
    [1, 0, 1, 2, 4, 2, 2, 1, 2, 4]

Running thousands of iterations on files with many features will of course
result in more complex results.  We could then take these results and plot
them in matplotlib, or get some statistics on them.

The method :meth:`BedTool.randomstats` does this for us, but requires NumPy
and SciPy to be installed.  This method also calls
:meth:`BedTool.randomintersection` for us, returning the summarized results
in a dictionary.

:meth:`BedTool.randomstats` takes the same arguments as
:meth:`BedTool.randomintersection`:


.. doctest::

    >>> results_dict = a.randomstats(b, iterations=100, shuffle_kwargs={'chrom': True}, debug=True)

The keys to this results dictionary are as follows (some are redundant,
I've found these keys useful for writing out to file):

:iterations: 
    the number of iterations we specified

:actual: 
    the number of intersections between then un-shuffled `a` and `b`

:file_a:
    the filename of `a`

:file_b: 
    the filename of `b`

:<file_a>: 
    the key is actully the filename of `a`, and the value is the number of
    features in `a`

:<file_b>: 
    the key is actually the filename of `b` and the value is the number of
    features in `b`

:self: 
    number of features in `a` (or "self"; same value as for <file_a>)

:other:
    number of features in `b` (or "other"; same value as for <file_b>)

:frac randomized above actual: fraction of iterations that had counts above the actual count

:frac randomized below actual:
    fraction of iterations that had counts below the actual count

:median randomized:
    the median of the distribution of randomized intersections

:normalized:
    the actual count divided by the median; can be considered as a score

:percentile:
    the percentile of actual within the distribution of randomized
    intersections; can be considered an empirical p-value

:upper 97.5th:
    the 97.5th percentile of the randomized distribution

:lower 2.5th:
    the 2.5th percentile of the randomized distribution

For example:

.. doctest::

    >>> keys = ['self', 'other', 'actual', 'median randomized', 'normalized', 'percentile']
    >>> for key in keys:
    ...     print('%s: %s' % (key, results_dict[key]))
    self: 4
    other: 2
    actual: 3
    median randomized: 2.0
    normalized: 1.5
    percentile: 90.0

Contributions toward improving this code or implementing other methods of
statistical testing are very welcome!


Other statistics
----------------
In practice, a comparison between two sets of features (say, two transcription
factors) with 1000 randomizations will have an empirical p-value of < 0.001.
That is, out of all the randomizations performed,  every single one had fewer
intersections than the original.  Of course the resolution of the p-value is
dependent on the number of randomizations:  the lowest nonzero p-value for
10000 iterations will be 0.0001.  Getting a non-zero p-value often requires
doing more randomizations than is practical (several million to tens of
millions).

That's where the enrichment score comes in.  The randomized intersections
typically have a normal distribution, but just in case, we take the median of
the randomized intersections and call this the background or control.  Then we
divide the actual intersections by this median to get an enrichment score.

The advantage to using the enrichment score is that it gives nonzero scores for
more fine-grained comparison among sets of features without performing
impractical amounts of randomization.  The first example of its usage that I'm
aware of is Negre et al. (2010) PLoS Genet 6(1): e1000814,  The downside of
this metric is that the numbers are relative, and have their greatest utility
for making biological conclusions when used in large matrices of pairwise
comparisons.

:meth:`BedTool.randomintersection` and :meth:`BedTool.randomstats` both use the
intersection count method.  That is, for each randomization the calculated
metric is "number of intersection events".  An alternative is to compute the
Jaccard statistic on each iteration, as implemented in
:meth:`BedTool.naive_jaccard`. The Jaccard statistic (or Jaccard similarity) is
the ratio of the intersection over the union, and is introduced in a genomic
intersection context in Favorov et al. (2012) PLoS Comput Biol 8(5): e1002529.
However, this still has the same p-value resolution limitation, so the
actual-divided-by-median approach could be tried here as well.
