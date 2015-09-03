.. include:: includeme.rst

.. _intersections:

Intersections
=============
One common use of BEDTools_ and :mod:`pybedtools` is to perform
intersections.

First, let's create some example :class:`BedTool` instances:

.. doctest::

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')

Then do the intersection with the :meth:`BedTool.intersect` method:

.. doctest::

    >>> a_and_b = a.intersect(b)

`a_and_b` is a new :class:`BedTool` instance.  It now points to a temp file
on disk, which is stored in the attribute `a_and_b.fn`; this temp file contains
the intersection of `a` and `b`.

We can either print the new :class:`BedTool` (which will show ALL features
-- use with caution if you have huge files!) or use the
:meth:`BedTool.head` method to show up to the first N lines (10 by
default).  Here's what `a`, `b`, and `a_and_b` look like:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> a.head()
    chr1    1   100 feature1    0   +
    chr1    100 200 feature2    0   +
    chr1    150 500 feature3    0   -
    chr1    900 950 feature4    0   +

    >>> b.head()
    chr1    155 200 feature5    0   -
    chr1    800 901 feature6    0   +

    >>> a_and_b.head()
    chr1    155 200 feature2    0   +
    chr1    155 200 feature3    0   -
    chr1    900 901 feature4    0   +

The :meth:`BedTool.intersect` method simply wraps the BEDTools_ program
`intersectBed`.  This means that we can pass :meth:`BedTool.intersect` any
arguments that `intersectBed` accepts.  For example, if we want to use the
`intersectBed` switch `-u` (which, according to the BEDTools documentation,
acts as a True/False switch to indicate that we want to see the features in `a`
that overlapped something in `b`), then we can use the keyword argument
`u=True`, like this:


.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> # Intersection using the -u switch
    >>> a_with_b = a.intersect(b, u=True)
    >>> a_with_b.head()
    chr1    100 200 feature2    0   +
    chr1    150 500 feature3    0   -
    chr1    900 950 feature4    0   +

This time, `a_with_b` is another :class:`BedTool` object that points to a
different temp file whose name is stored in `a_with_b.fn`.  You can read
more about the use of temp files in :ref:`temp principle`.  More on
arguments that you can pass to :class:`BedTool` objects in a moment, but
first, some info about saving files.

