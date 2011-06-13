.. _bam:

Working with BAM files
======================
Some BEDTools programs support BAM files as input; for example
`intersectBed`, `windowBed`, and others accept a `-abam` argument instead
of `-a` for the first input file.

If you create a :class:`BamTool` out of a BAM file, like:

.. doctest::

    x = pybedtools.example_bedtool('x.bam')

then any methods that support BAM as input via the `-abam` or `-ibam`
arguments will automatically use those, just as `-a` or `-i` arguments will
automatically be used for BED or GFF files.

The output follows the semantics of BEDTools.  That is, for programs like
`intersectBed`, if `abam` is used then the output will be BAM format as
well.  But if the `-bed` argument is passed, then the output will be BED
format.

As an example, let's intersect a BAM file of reads with annotations using
files that ship with :mod:`pybedtools`.  First, we create the
:class:`BedTool` objects:

.. doctest::

    >>> a = pybedtools.example_bedtool('x.bam')
    >>> b = pybedtools.example_bedtool('dm3-chr2L-5M.gff.gz')

The following will return BAM results:

.. doctest::

    >>> c = a.intersect(b)

We can iterate over BAM files to get :class:`Interval` objects just like
iterating over BED or GFF files.  Indexing works, too:

.. doctest::
    :options: +ELLIPSIS +NORMALIZE_WHITESPACE

    >>> for i in c[:2]:
    ...     print i
    HWUSI-NAME:2:69:512:1017#0	16	chr2L	9330	3	36M	*	0	0	TACAAATCTTACGTAAACACTCCAAGCATGAATTCG	Y`V_a_TM[\_V`abb`^^Q]QZaaaaa_aaaaaaa	NM:i:0	NH:i:2	CC:Z:chrX	CP:i:19096815
    HWUSI-NAME:2:91:1201:1113#0	16	chr2L	10213	255	36M	*	0	0	TGTAGAATGCAAAAATTACATTTGTGAGTATCATCA	UV[aY`]\VZ`baaaZa`_aab_`_`a`ab``b`aa	NM:i:0	NH:i:1

    >>> c[0]
    Interval(chr2L:9329-9365)

    >>> c[:10]
    <itertools.islice object at ...>

    >>> cigar_string = i[5]

Note that :mod:`pybedtools` uses the convention that BAM features in plain
text format are considered SAM features, so these SAM features are
**one-based and include the stop coordinate** as illustrated below:

.. doctest::

    >>> c[0].start
    9329L

    >>> c[0][3]
    '9330'


Currently, the stop coordinate is defined as the start coord plus the
length of the sequence; eventually a more sophisticated, CIGAR-aware
approach may be used.  Similarly, the length is defined to be `stop -
start`, again, not CIGAR-aware at the moment.  For more sophisticated
low-level manipulation of BAM features, you might want to consider HTSeq_.

Alternatively, we can specify the `bed=True` kwarg to convert the
intersected BAM results to BED format, and use those like a normal BED
file:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> d = a.intersect(b, bed=True)
    >>> d.head(3)
    chr2L	9329	9365	HWUSI-NAME:2:69:512:1017#0	3	-
    chr2L	9329	9365	HWUSI-NAME:2:69:512:1017#0	3	-
    chr2L	9329	9365	HWUSI-NAME:2:69:512:1017#0	3	-


.. _HTSeq: http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
