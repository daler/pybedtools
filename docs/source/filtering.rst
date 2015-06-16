
.. include:: includeme.rst

.. _filtering:

Filtering
~~~~~~~~~
The :meth:`BedTool.filter` method lets you pass in a function that accepts an
:class:`Interval` as its first argument and returns True for False.  This
allows you to perform "grep"-like operations on :class:`BedTool` objects.  For
example, here's how to get a new :class:`BedTool` containing features from `a`
that are more than 100 bp long:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = a.filter(lambda x: len(x) > 100)
    >>> print(b)
    chr1	150	500	feature3	0	-
    <BLANKLINE>

The :meth:`filter` method will pass its `*args` and `**kwargs` to the function
provided.  So here is a more generic case, where the function is defined once
and different arguments are passed in for filtering on different lengths:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> def len_filter(feature, L):
    ...     "Returns True if feature is longer than L"
    ...     return len(feature) > L

Now we can pass different lengths without defining a new function for each
length of interest, like this:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> a = pybedtools.example_bedtool('a.bed')

    >>> print(a.filter(len_filter, L=10))
    chr1	1	100	feature1	0	+
    chr1	100	200	feature2	0	+
    chr1	150	500	feature3	0	-
    chr1	900	950	feature4	0	+
    <BLANKLINE>

    >>> print(a.filter(len_filter, L=99))
    chr1	100	200	feature2	0	+
    chr1	150	500	feature3	0	-
    <BLANKLINE>

    >>> print(a.filter(len_filter, L=200))
    chr1	150	500	feature3	0	-
    <BLANKLINE>


See :ref:`BedTools as iterators` for more advanced and space-efficient usage
of :meth:`filter` using iterators.

Note that we could have used the built-in Python function, `filter()`, but that
would have returned an iterator that we would have to construct a new
:class:`pybedtools.BedTool` out of.  The :meth:`BedTool.filter` method returns
a ready-to-use :class:`BedTool` object, which allows embedding of
:meth:`BedTool.filter` calls in a chain of commands, e.g.::

    >>> a.intersect(b).filter(lambda x: len(x) < 100).merge()

Fast filtering functions in Cython
----------------------------------

The :mod:`featurefuncs` module contains some ready-made functions written
in Cython that will be faster than pure Python equivalents.  For example,
there are :func:`greater_than` and :func:`less_than` functions, which are
about 70% faster.  In IPython::

    >>> from pybedtools.featurefuncs import greater_than

    >>> len(a)
    310456

    >>> def L(x,width=100):
    ...     return len(x) > 100

    >>> # The %timeit command is from IPython, and won't work
    >>> # in a regular Python script:
    >>> %timeit a.filter(greater_than, 100)
    1 loops, best of 3: 1.74 s per loop

    >>> %timeit a.filter(L, 100)
    1 loops, best of 3: 2.96 s per loop

