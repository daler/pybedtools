.. include:: includeme.rst

Filtering
~~~~~~~~~
The :meth:`filter` method lets you pass in a function that accepts an :class:`Interval` as its first
argument and returns True for False.  For example, to only get features of a certain size:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = a.filter(lambda x: len(x) > 100)
    >>> print b
    chr1	150	500	feature3	0	-
    <BLANKLINE>

The :meth:`filter` method will pass its `*args` and `**kwargs` to the function
provided.  So a more generic case would be the following, where the function is defined once
and different arguments are passed in for filtering on different lengths:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> def len_filter(feature, L):
    ...     return len(feature) > L

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> print a.filter(len_filter, L=10)
    chr1	1	100	feature1	0	+
    chr1	100	200	feature2	0	+
    chr1	150	500	feature3	0	-
    chr1	900	950	feature4	0	+
    <BLANKLINE>

    >>> print a.filter(len_filter, L=99)
    chr1	100	200	feature2	0	+
    chr1	150	500	feature3	0	-
    <BLANKLINE>

    >>> print a.filter(len_filter, L=200)
    chr1	150	500	feature3	0	-
    <BLANKLINE>

The :meth:`filter` method uses a file-based format, where the new
:class:`BedTool` object refers to a new temp file.  You can use a generator
function to create a new :class:`BedTool` if you want to save disk space:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.BedTool( (i for i in a if len_filter(i, L=200)))
    >>> print b
    chr1	150	500	feature3	0	-
    <BLANKLINE>

However, keep in mind that printing `b`, which was created using a
generator expression, has now been consumed -- so printing `b` again will
do nothing:

.. doctest::

   >>> print b
   <BLANKLINE>
   <BLANKLINE>

If you create a :class:`BedTool` with a generator expression, you can
always save it as a file for later use. This is what :meth:`filter` is
doing:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.BedTool( (i for i in a if len_filter(i, L=200))).saveas('len-filtered-b.bed')
    >>> print b
    chr1	150	500	feature3	0	-
    <BLANKLINE>

    >>> print b
    chr1	150	500	feature3	0	-
    <BLANKLINE>

The :mod:`featurefuncs` module contains some ready-made functions written
in Cython that will be faster than pure Python equivalents.  For example,
there are :func:`greater_than` and :func:`less_than` functions, which are
about 70% faster.  In IPython::

    >>> len(a)
    310456

    >>> def L(x,width=100):
    ...     return len(x) > 100

    >>> %timeit a.filter(greater_than, 100)
    1 loops, best of 3: 1.74 s per loop

    >>> %timeit a.filter(L, 100)
    1 loops, best of 3: 2.96 s per loop
