.. include:: includeme.rst

Saving :class:`BedTool` results
===============================
In general, there are three different ways of saving results from
:class:`BedTool` operations:

Use the :meth:`BedTool.saveas` method
-------------------------------------
The :meth:`BedTool.saveas` method makes a **copy** of the results, so beware
that for large files, this can be time and/or memory-consuming.  However, when
working with a streaming or iterating :class:`BedTool`, this is a great way to
render the results to disk in the middle of a pipeline.

A good example of this is saving the results from a :meth:`BedTool.each` call:


.. doctest::

    >>> from pybedtools.featurefuncs import TSS
    >>> a = pybedtools.example_bedtool('a.bed')
    >>> result = a.each(TSS, upstream=1000, downstream=0)\
    ...     .saveas('upstream_regions.bed')

Use the :meth:`BedTool.moveto` method
-------------------------------------
The :meth:`BedTool.moveto` method does a **move** operation of the results.
This is best used when the results have been written to disk already (perhaps
to a tempfile) but you'd like to give the file a more reasonable/memorable
name.

.. doctest::

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')
    >>> c = a.intersect(b).moveto('intersection_of_a_and_b.bed')


Use the ``output`` keyword argument
-----------------------------------
If you know ahead of time that you want to save the output to a particular
file, use the ``output`` keyword argument to any wrapped :class:`BedTool`
method that returns another :class:`BedTool` object.  This will override the
default behavior of creating a tempfile.

.. doctest::


    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')
    >>> c = a.intersect(b, output='intersection_of_a_and_b.bed')


Working with non-interval output files
--------------------------------------
`BEDTools` commands offer lots of flexibility. This means it is possible to
return results that are not supported interval files like
BED/GFF/GTF/BAM/SAM/VCF.

Consider the following example, which uses :meth:`BedTool.groupby` to get
a 2-column file containing the number of intervals in each featuretype:

.. doctest::

    >>> a = pybedtools.example_bedtool('gdc.gff')
    >>> b = pybedtools.example_bedtool('gdc.bed')
    >>> c = a.intersect(b, c=True)
    >>> d = c.groupby(g=[3], c=10, o=['sum'])

We can read the file created by `d` looks like this:

(note: the latest version of BEDTools, v2.26.0, causes this to fail. This will
be fixed in the next BEDTools release (see
https://github.com/arq5x/bedtools2/issues/453,
https://github.com/arq5x/bedtools2/issues/450,
https://github.com/arq5x/bedtools2/issues/435,
https://github.com/arq5x/bedtools2/issues/436 for details).

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> # bedtools v2.26.0
    >>> print(open(d.fn).read())
    UTR	0
    CDS	2
    intron	4
    CDS	0
    UTR	1
    exon	3
    mRNA	7
    CDS	2
    exon	2
    tRNA	2
    gene	7
    <BLANKLINE>


Trying to iterate over `d` (`[i for i in d]`) or save it (`d.saveas()`) raises
exceptions. This is because:

* `saveas()` is expected to return a `BedTool` object that can be
  used with other `BEDTools` tools. We can't create a `BedTool` object out of
  an unsupported file format like this

* iterating over a `BedTool` object is expected to yield `Interval` objects,
  but these lines can't be converted into the supported formats


To save the output to a filename of your choosing, provide the `output`
argument instead of `saveas()`, like this:

.. doctest::

    >>> # only works with bedtools != v2.26.0
    >>> # d = c.groupby(g=[3], c=10, o=['sum'], output='counts.txt')

To iterate over the lines of the file, you can use standard Python
tools, e.g.:

.. doctest::

    >>> # only works with bedtools != v2.26.0
    >>> # for line in open(d.fn):
    >>> #     featuretype, count = line.strip().split()
