
.. _BEDTools: http://github.com/arq5x/bedtools


.. _3examples:

Three brief examples
--------------------
Here are three examples to show typical usage of :mod:`pybedtools`.  More
info can be found in the docstrings of :mod:`pybedtools` methods and in the
:ref:`tutorial`. 

You can also check out :ref:`shell_comparison` for a simple
example of how :mod:`pybedtools` can improve readability of your code with no
loss of speed compared to bash scripting.

.. note::

    Please take the time to read and understand the conventions
    :mod:`pybedtools` uses to handle files with different coordinate systems
    (e.g., 0-based BED files vs 1-based GFF files) which are described
    :ref:`here <zero_based_coords>`.

    In summary,

    * **Integer** values representing start/stop are *always in 0-based
      coordinates*, regardless of file format.  This means that all
      :class:`Interval` objects can be treated identically, and greatly
      simplifies underlying code.

    * **String** values representing start/stop will use coordinates appropriate
      for the format (1-based for GFF; 0-based for BED).

Example 1: Save a BED file of intersections, with track line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example saves a new BED file of intersections between your files `mydata/snps.bed` and
`mydata/exons.bed`, adding a track line to the output::

    >>> import pybedtools
    >>> a = pybedtools.BedTool('mydata/snps.bed')
    >>> a.intersect('mydata/exons.bed').saveas('snps-in-exons.bed', trackline="track name='SNPs in exons' color=128,0,0")

Example 2: Intersections for a 3-way Venn diagram
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This example gets values for a 3-way Venn diagram of overlaps.  This
demonstrates operator overloading of :class:`BedTool` objects.  It assumes that
you have the files `a.bed`, `b.bed`, and `c.bed` in your current working
directory.  If you'd like to use example files that come with
:mod:`pybedtools`, then replace strings like `'a.bed'` with
`pybedtools.example_filename('a.bed')`, which will retrieve the absolute path
to the example data file.::


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
Charts API respectively.  Also see the :mod:`pybedtools.contrib.venn_maker`
module for a flexible interface to the VennDiagram `R` package.

.. _third example:

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

.. literalinclude:: example_3

For more on using :mod:`pybedtools`, continue on to the :ref:`tutorial` . .  .
