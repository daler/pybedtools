.. _bam:

Working with BAM files
======================
Some BEDTools programs support BAM files as input; for example
`intersectBed`, `windowBed`, and others accept a `-abam` argument instead
of `-a` for the first input file.

This section describes the workflow for working with BAM files within
:mod:`pybedtools`.

As an example, let's intersect a BAM file of reads with annotations using
files that ship with :mod:`pybedtools`.  First, we create the
:class:`BedTool` objects:

.. doctest::

    >>> a = pybedtools.example_bedtool('x.bam')
    >>> b = pybedtools.example_bedtool('dm3-chr2L-5M.gff.gz')

If `a` referred to a BED file like `a.bed`, we could just do
`a.intersect(b)` because `a.bed` would be implictly passed as `-a` and the
gzipped GFF file would be passed as `-b`.  In order to use a BAM file,
however, we need to explicitly specify an `abam` kwarg.  In addition, since
Python doesn't allow non-keyword arguments after keyword arguments, we need
to explicitly specify a `b` kwarg.  This should be much clearer with a
simple example:

.. doctest::

    >>> c = a.intersect(abam=a.fn, b=b)

Now `c` points to a new BAM file on disk.  Keep in mind that there is not
yet iterable BAM support in :mod:`pybedtools`, so things like `c.count()`
or iterating over `c` with a `for feature in c: ...` won't work.  For now,
consider using a package like HTSeq_ for access to individual reads in BAM
format.

Alternatively, we can specify the `bed=True` kwarg to convert the
intersected BAM results to BED format, and use those like a normal BED
file:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> d = a.intersect(abam=a.fn, b=b, bed=True)

The resulting BedTool `d` refers to a BED file and can be used like any other:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> d.count()
    341324


    >>> print iter(d).next()
    chr2L	9329	9365	HWUSI-NAME:2:69:512:1017#0	3	-


.. _HTSeq: http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
