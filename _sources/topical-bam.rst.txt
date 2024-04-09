.. include:: includeme.rst

.. _bam:

Working with BAM files
======================
Some BEDTools programs, like `intersecteBed`, support BAM files as input.
From the command line, you would need to specify the `-abam`
argument to do so.  However, :mod:`pybedtools` auto-detects BAM files and
passes the `abam` argument automatically for you.  That means if you create
a :class:`BedTool` out of a BAM file, like this:

.. doctest::

    x = pybedtools.example_bedtool('gdc.bam')

you can intersect it with a BED file without doing anything special:

.. doctest::

    b = pybedtools.example_bedtool('gdc.gff')
    y = x.intersect(b)

The output of this operation follows the semantics of BEDTools.  That is,
for programs like `intersectBed`, if `abam` is used then the output will be
BAM format as well.  But if the `-bed` argument is passed, then the output
will be BED format. Similarly, in :mod:`pybedtools`, if a BAM file is used
to create the :class:`BedTool` then the results will also be in BAM
format.  If the `bed=True` kwarg is passed, then the results be in BED
format.

As an example, let's intersect a BAM file of reads with annotations using
files that ship with :mod:`pybedtools`.  First, we create the
:class:`BedTool` objects:

.. doctest::

    >>> a = pybedtools.example_bedtool('x.bam')
    >>> b = pybedtools.example_bedtool('dm3-chr2L-5M.gff.gz')

The first call below will return BAM results, and the second will return
BED results.

.. doctest::

    >>> bam_results = a.intersect(b)
    >>> str(bam_results.file_type) == 'bam'
    True

    >>> bed_results = a.intersect(b, bed=True)
    >>> str(bed_results.file_type) == 'bed'
    True


We can iterate over BAM files to get :class:`Interval` objects just like
iterating over BED or GFF files.  Indexing works, too:

.. doctest::
    :options: +ELLIPSIS +NORMALIZE_WHITESPACE

    >>> for i in bam_results[:2]:
    ...     print(i)
    HWUSI-NAME:2:69:512:1017#0	16	chr2L	9330	3	36M	*	0	0	TACAAATCTTACGTAAACACTCCAAGCATGAATTCG	Y`V_a_TM[\_V`abb`^^Q]QZaaaaa_aaaaaaa	NM:i:0	NH:i:2	CC:Z:chrX	CP:i:19096815
    <BLANKLINE>
    HWUSI-NAME:2:91:1201:1113#0	16	chr2L	10213	255	36M	*	0	0	TGTAGAATGCAAAAATTACATTTGTGAGTATCATCA	UV[aY`]\VZ`baaaZa`_aab_`_`a`ab``b`aa	NM:i:0	NH:i:1
    <BLANKLINE>

    >>> bam_results[0]
    Interval(chr2L:9329-9365)

    >>> bam_results[:10]
    <itertools.islice object at ...>

    >>> cigar_string = i[5]

There are several things to watch out for here.

First, note that :mod:`pybedtools` uses the convention that BAM features in
plain text format are considered SAM features, so these SAM features are
**one-based and include the stop coordinate** as illustrated below. (Note that
there is some additional complexity here due to supporting Python 2 and
3 simultaneously in this tested documentation)

.. doctest::

    >>> bam_results[0].start
    9329

    >>> isinstance(bam_results[0][3], str)
    True

    >>> print(bam_results[0][3])
    9330


Second, the stop coordinate is defined as the *start coord plus the
length of the sequence*; eventually a more sophisticated, CIGAR-aware
approach may be used.  Similarly, the length is defined to be `stop -
start` -- again, not CIGAR-aware at the moment.  For more sophisticated
low-level manipulation of BAM features, you might want to consider using
HTSeq_.

Third, while we can iterate over a BAM file and manipulate the features as
shown above, *calling BEDTools programs on a BAM-based generator is not
well-supported*.

Specifically::

    >>> a = pybedtools.example_bed('gdc.bam')
    >>> b = pybedtools.example_bed('b.bed')

    >>> # works, gets BAM results
    >>> results = a.intersect(b)

    >>> # make a generator of features in `a`
    >>> a2 = pybedtools.BedTool(i for i in a)

    >>> # this does NOT work
    >>> a2.intersect(b)

When we specified the `bed=True` kwarg above, the intersected BAM results
are converted to BED format.  We can use those like a normal BED file.
Note that since we are viewing BED output, *the start and stops are 0-based*:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> d = a.intersect(b, bed=True)
    >>> d.head(3)
    chr2L	9329	9365	HWUSI-NAME:2:69:512:1017#0	3	-	9329	9365	0,0,0	1	36,	0,
    chr2L	9329	9365	HWUSI-NAME:2:69:512:1017#0	3	-	9329	9365	0,0,0	1	36,	0,
    chr2L	9329	9365	HWUSI-NAME:2:69:512:1017#0	3	-	9329	9365	0,0,0	1	36,	0,

Consistent with BEDTools programs, BAM files are **not** supported as the
second input argument.  In other words, `intersectBed` does not have both
`-abam` and `-bbam` arguments, so :mod:`pybedtools` will not not allow this
either.

However, :mod:`pybedtools` does allow streaming BAM files to be the input of
methods that allow BAM input as the first input. In this [trivial] example, we
can stream the first intersection to save disk space, and then send that
streaming BAM to the next :meth:`BedTool.intersect` call. Since it's not
streamed, the second intersection will be saved as a temp BAM file on disk::


    >>> a.intersect(b, stream=True).intersect(b)

.. _HTSeq: http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
