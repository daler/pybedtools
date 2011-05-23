.. _`BedTools as iterators`:

Using BedTool objects as iterators/generators
=============================================

Typically, :mod:`BedTool` objects are used somewhat like handles to individual
files on disk that contain BED lines.  To save disk space, :mod:`BedTool`
objects also have the ability to "stream", much like piping in Unix.  That
is, the data are created only one line at a time in memory, instead of
either creating a list of all data in memory or writing all data to disk.

.. warning::

    You'll need to be careful when using :mod:`BedTool` objects as
    generators, since any operation that reads all the features of a
    :mod:`BedTool` will consume the iterable.

To get a streaming BedTool, use the `stream=True` kwarg.  This
:class:`BedTool` will act a little differently from a standard, file-based
:class:`BedTool`.

.. doctest::

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')
    >>> c = a.intersect(b, stream=True)

    >>> # checking the length consumes the iterator
    >>> len(c)
    3

    >>> # nothing left, so checking length again returns 0
    >>> len(c)
    0

In some cases, a stream may be "rendered" to a temp file.  This is because
BEDTools programs can only accept one input file as `stdin`.  This is typically
the first input (`-i` or `-a`), while the other input (`-b`) must be a file.
Consider this example, where the second intersection needs to convert the
streaming BedTool to a file before sending to `intersectBed`:

.. doctest::

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')

    >>> # first we set up a streaming BedTool:
    >>> c = a.intersect(b, stream=True)

    >>> # But supplying a streaming BedTool as the first unnamed argument
    >>> # means it is being passed as -b to intersectBed, and so must be a file.
    >>> # In this case, `c` is rendered to a tempfile before being passed.
    >>> d = a.intersect(c, stream=True)


Creating a :class:`BedTool` from an iterable
--------------------------------------------
You can create a :class:`BedTool` on the fly from a generator or iterator -- in
fact, this is what the :meth:`BedTool.filter` method does for you:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> print a
    chr1	1	100	feature1	0	+
    chr1	100	200	feature2	0	+
    chr1	150	500	feature3	0	-
    chr1	900	950	feature4	0	+
    <BLANKLINE>

    >>> b = pybedtools.BedTool(f for f in a if f.start > 200)

    >>> # this is the same as using filter:
    >>> c = a.filter(lambda x: x.start > 200)

We need to "render" these BedTools to string before we can check equality
-- consuming them both -- since they are both iterables for which `==` is
not defined:

.. doctest::
    :options: +ELLIPSIS

    >>> b == c
    Traceback (most recent call last):
        ...
    NotImplementedError: Testing equality only supported for BedTools that point to a file

    >>> str(b) == str(c)
    True
