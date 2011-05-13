.. _bam:

Working with BAM files
======================
Some BEDTools programs support BAM files as input, for example
`intersectBed` accepts a `-abam` argument.

This section describes the workflow for working with BAM files within
:mod:`pybedtools`.

As an example, let's intersect a BAM file of reads with annotations using
files that ship with :mod:`pybedtools`.  First, we create the
:class:`BedTool` objects:

.. doctest::

    >>> a = pybedtools.example_bedtool('x.bam')
    >>> b = pybedtools.example_bedtool('dm3-chr2L-5M.gff.gz')


The trick is to use only keyword arguments, and not let :mod:`pybedtools`
make assumptions about what should be passed as `-a` and `-b`:

.. doctest::

    >>> c = a.intersect(abam=a.fn, b=b)

There is not yet iterable BAM support in :mod:`pybedtools`, so things like
`c.count()` or iterating over `c` with a `for feature in c: ...` won't
work.  Consider using a package like HTSeq_ for access to individual reads
in BAM format.

However, we can specify the `bed=True` kwarg to convert the intersected BAM
results to BED format, and use those like a normal BED file:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> d = a.intersect(abam=a.fn, b=b, bed=True)
    >>> d.count()
    341324

    >>> print iter(d).next()
    chr2L	9329	9365	HWUSI-NAME:2:69:512:1017#0	3	-


.. _HTSeq: http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
