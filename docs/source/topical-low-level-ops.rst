.. include:: includeme.rst

Low-level operations
--------------------

We can use the :meth:`BedTool.as_intervalfile` method to return an
:class:`IntervalFile` instance.  This class provides low-level support to
the BEDTools C++ API.

The method :meth:`IntervalFile.all_hits` takes a single :class:`Interval`
as the query and returns a list of all features in the
:class:`IntervalFile` that intersect:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> ivf = a.as_intervalfile()
    >>> query = a[2]
    >>> ivf.all_hits(query)
    [Interval(chr1:100-200), Interval(chr1:150-500)]

Similarly, we can just return if there were *any* hits, a much faster
operation:

.. doctest:: 
    :options: +NORMALIZE_WHITESPACE

    >>> ivf.any_hits(query)
    1

Or count how many hits:

.. doctest:: 
    :options: +NORMALIZE_WHITESPACE

    >>> ivf.count_hits(query)
    2

See the docstrings for :meth:`IntervalFile.all_hits`,
:meth:`IntervalFile.any_hits`, and :meth:`IntervalFile.count_hits` for
more, including stranded hits and restricting hits to a specific overlap.

.. note::

    These methods are now available as :class:`BedTool` methods,
    :meth:`BedTool.all_hits`, :meth:`BedTool.any_hits`, and
    :meth:`BedTool.count_hits`
