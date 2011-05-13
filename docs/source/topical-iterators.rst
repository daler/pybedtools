.. _`BedTools as iterators`:

Using BedTool objects as iterators
==================================

Typically, :mod:`BedTool` objects are used somewhat like handles to individual
files on disk that contain BED lines.  To save disk space, :mod:`BedTool`
objects have the ability to "stream", much like piping in Unix.

You'll need to be careful when using :mod:`BedTool` objects as generators,
since any operation that reads all the features of a :mod:`BedTool` will
consume the iterable.

To get a streaming BedTool, use the `stream=True` kwarg:

.. doctest::

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')
    >>> c = a.intersect(b, stream=True)

    # checking the length consumes the iterator
    >>> len(c)
    3

    # nothing left, so checking length again returns 0
    >>> len(c)
    0
