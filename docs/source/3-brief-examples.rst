
.. _BEDTools: http://github.com/arq5x/bedtools


.. _3examples:

Three brief examples
--------------------
Here are three examples to show typical usage of :mod:`pybedtools`.  More
info can be found in the docstrings of :mod:`pybedtools` methods and in the
:ref:`tutorial`.

Example 1: Save a BED file of intersections, with track line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example saves a new BED file of intersections between `a.bed` and
`b.bed`, adding a track line to the output::

    >>> import pybedtools
    >>> a = pybedtools.BedTool('a.bed')
    >>> a.intersect('b.bed').saveas('a-and-b.bed', trackline="track name='a and b' color=128,0,0")

Example 2: Intersections for a 3-way Venn diagram
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example gets values for a 3-way Venn diagram of overlaps.  This
demonstrates operator overloading of :class:`BedTool` objects::


    >>> import pybedtools

    >>> # set up 3 different bedtools
    >>> a = pybedtools.BedTool('a.bed')
    >>> b = pybedtools.BedTool('b.bed')
    >>> c = pybedtools.BedTool('c.bed')

    >>> (a-b-c).count()  # unique to a
    >>> (a+b-c).count()  # in a and b, not c
    >>> (a+b+c).count()  # common to all 
    >>> # ... and so on, for all the combinations.



For more, see the :mod:`pybedtools.scripts.venn_mpl` and
:mod:`pybedtools.scripts.venn_gchart` scripts, which wrap this functionality in
command-line scripts to create Venn diagrams using either matplotlib or Google
Charts API respectively.

Example 3: Count reads in introns and exons, in parallel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example shows how to count the number of reads in introns and exons in
parallel. It is somewhat more involved, but illustrates several additional
features of :mod:`pybedtools` such as:

* BAM file support (for more, see :ref:`bam`)
* indexing into Interval objects (for more, see :ref:`intervals`)
* filtering (for more, see :ref:`filtering`)
* streaming (for more, see :ref:`BedTools as iterators`)
* ability to use parallel processing

The first listing has many explanatory comments, and the second listing shows
the same code with no comments to give more of a feel for :mod:`pybedtools`.

.. literalinclude:: example_3

Here's the same code but with no comments:

.. literalinclude:: example_3_no_comments

For more on using :mod:`pybedtools`, continue on to the :ref:`tutorial` . .
.
