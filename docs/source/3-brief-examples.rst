
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

Example 3: Classify reads into intron/exon
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example is a little more involved, but highlights some useful features
of :mod:`pybedtools`.  Here, we classify reads in a BAM file into exon,
intron, exon and intron (i.e., exonic in one isoform but intronic in
anotheroverlapping isoform) or intergenic classes in a stranded manner.  A
more generalized version of this example can be found in the
:mod:`pybedtools.scripts.classify_reads` script.


.. literalinclude:: example-script

Here is the same script without comments, to give a better feel for
:mod:`pybedtools`:

.. literalinclude:: example-script-nocomments

For more, continue on to the :ref:`tutorial`, and then check out the :ref:`topical`.

