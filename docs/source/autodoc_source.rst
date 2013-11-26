
.. _autodoc:

.. _pybedtools reference:

:mod:`pybedtools` Reference
===========================
The following tables summarize the methods and functions; click on a method or
function name to see the complete documentation.

.. contents::

:class:`BedTool` creation
-------------------------
The main :class:`BedTool` documentation, with a list of all methods in
alphabetical order at the bottom.  For more details, please see :ref:`creating
a BedTool`.

.. autosummary::
    :toctree: autodocs

    pybedtools.BedTool

`BEDTools` wrappers
-------------------
These methods wrap `BEDTools` programs for easy use with Python; you can then
use the other :mod:`pybedtools` functionality for further manipulation and
analysis.

The documentation of each of these methods starts with
:mod:`pybedtools`-specific documentation, possibly followed by an example.
Finally, the `BEDTools` help is copied verbatim from whatever version was
installed when generating these docs.

In general the `BEDTool` wrapper methods adhere to the :ref:`Design principles`:

    * :ref:`temp principle`
    * :ref:`similarity principle`
    * :ref:`version principle`
    * :ref:`default args principle`

A new interface was introduced in BEDTools v2.15 which retains compatibility
with :mod:`pybedtools`.  For clarity, in the table below, both the "old" (e.g.,
`intersectBed`) or "new" (e.g., `bedtools intersect`) versions of calling the
program are indicated. 

.. autosummary::
    :toctree: autodocs

    pybedtools.BedTool.intersect
    pybedtools.BedTool.window
    pybedtools.BedTool.closest
    pybedtools.BedTool.coverage
    pybedtools.BedTool.map
    pybedtools.BedTool.genome_coverage
    pybedtools.BedTool.merge
    pybedtools.BedTool.cluster
    pybedtools.BedTool.complement
    pybedtools.BedTool.subtract
    pybedtools.BedTool.slop
    pybedtools.BedTool.flank
    pybedtools.BedTool.sort
    pybedtools.BedTool.random
    pybedtools.BedTool.shuffle
    pybedtools.BedTool.annotate
    pybedtools.BedTool.multi_intersect
    pybedtools.BedTool.union_bedgraphs
    pybedtools.BedTool.pair_to_bed
    pybedtools.BedTool.pair_to_pair
    pybedtools.BedTool.bam_to_bed
    pybedtools.BedTool.to_bam
    pybedtools.BedTool.bedpe_to_bam
    pybedtools.BedTool.bed6
    pybedtools.BedTool.bam_to_fastq
    pybedtools.BedTool.sequence
    pybedtools.BedTool.mask_fasta
    pybedtools.BedTool.nucleotide_content
    pybedtools.BedTool.multi_bam_coverage
    pybedtools.BedTool.tag_bam
    pybedtools.BedTool.jaccard
    pybedtools.BedTool.reldist
    pybedtools.BedTool.overlap
    pybedtools.BedTool.links
    pybedtools.BedTool.igv
    pybedtools.BedTool.window_maker
    pybedtools.BedTool.groupby
    pybedtools.BedTool.expand

Other :class:`BedTool` methods
------------------------------
These methods are some of the ways in which :mod:`pybedtools` extend the
BEDTools suite.


Feature-by-feature operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Methods that operate on a feature-by-feature basis to modify or filter features
on the fly.

.. autosummary::
    :toctree: autodocs

    pybedtools.BedTool.each
    pybedtools.BedTool.filter
    pybedtools.BedTool.split
    pybedtools.BedTool.truncate_to_chrom
    pybedtools.BedTool.remove_invalid

The :mod:`pybedtools.featurefuncs` module contains some commonly-used functions
that can be passed to :meth:`BedTool.each`:

.. autosummary::
    :toctree:

    pybedtools.featurefuncs.three_prime
    pybedtools.featurefuncs.five_prime
    pybedtools.featurefuncs.TSS
    pybedtools.featurefuncs.extend_fields
    pybedtools.featurefuncs.center
    pybedtools.featurefuncs.midpoint
    pybedtools.featurefuncs.normalized_to_length
    pybedtools.featurefuncs.rename
    pybedtools.featurefuncs.greater_than
    pybedtools.featurefuncs.less_than
    pybedtools.featurefuncs.normalized_to_length
    pybedtools.featurefuncs.rename
    pybedtools.featurefuncs.bedgraph_scale
    pybedtools.featurefuncs.add_color
    pybedtools.featurefuncs.gff2bed
    pybedtools.featurefuncs.bed2gff



Searching for features
~~~~~~~~~~~~~~~~~~~~~~
These methods take a single interval as input and return the intervals of the
BedTool that overlap.

This can be useful when searching across many BED files for a particular
coordinate range -- for example, they can be used identify all binding sites,
stored in many different BED files, that fall within a gene's coordinates.

.. autosummary::
    :toctree: autodocs

    pybedtools.BedTool.all_hits
    pybedtools.BedTool.any_hits
    pybedtools.BedTool.count_hits
    pybedtools.BedTool.tabix_intervals
    pybedtools.BedTool.tabix
    pybedtools.BedTool.bgzip


:class:`BedTool` introspection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
These methods provide information on the :class:`BedTool` object.

If using :meth:`BedTool.head`, don't forget that you can index into
:class:`BedTool` objects, too.

.. autosummary::
    :toctree: autodocs

    pybedtools.BedTool.head
    pybedtools.BedTool.count
    pybedtools.BedTool.field_count
    pybedtools.BedTool.file_type


Randomization helpers
~~~~~~~~~~~~~~~~~~~~~
Helper methods useful for assessing empirical instersection
distributions between interval files.

.. autosummary::
    :toctree: autodocs

    pybedtools.BedTool.parallel_apply
    pybedtools.BedTool.randomstats
    pybedtools.BedTool.randomintersection
    pybedtools.BedTool.randomintersection_bp
    pybedtools.BedTool.random_subset
    pybedtools.BedTool.random_jaccard
    pybedtools.BedTool.random_op

Managing :class:`BedTool` objects on disk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
These methods are used to specify where to save results from :class:`BedTool`
operations.

.. autosummary::
    :toctree: autodocs

    pybedtools.BedTool.saveas
    pybedtools.BedTool.moveto


Misc operations
~~~~~~~~~~~~~~~
Methods that can't quite be categorized into the above sections.

.. autosummary::
    :toctree: autodocs

    pybedtools.BedTool.cat
    pybedtools.BedTool.at
    pybedtools.BedTool.absolute_distance
    pybedtools.BedTool.cut
    pybedtools.BedTool.total_coverage
    pybedtools.BedTool.with_attrs
    pybedtools.BedTool.as_intervalfile
    pybedtools.BedTool.introns
    pybedtools.BedTool.set_chromsizes
    pybedtools.BedTool.print_sequence
    pybedtools.BedTool.save_seqs
    pybedtools.BedTool.seq
    pybedtools.BedTool.liftover
    pybedtools.BedTool.colormap_normalize
    pybedtools.BedTool.relative_distance

Module-level functions
----------------------

Working with example files
~~~~~~~~~~~~~~~~~~~~~~~~~~
:mod:`pybedtools` comes with many example files.  Here are some useful
functions for accessing them.

.. autosummary::
    :toctree: autodocs

    pybedtools.example_bedtool
    pybedtools.list_example_files
    pybedtools.example_filename

Creating :class:`Interval` objects from scratch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:class:`Interval` objects are the core object in :mod:`pybedtools` to represent
a genomic interval, written in Cython for speed.

.. autosummary::
    :toctree: autodocs

    pybedtools.Interval
    pybedtools.create_interval_from_list

:mod:`pybedtools` setup and config
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Use these functions right after importing in order to use custom paths or to
clean up the temp directory.

.. autosummary::
    :toctree: autodocs

    pybedtools.set_bedtools_path
    pybedtools.set_samtools_path
    pybedtools.get_tempdir
    pybedtools.set_tempdir
    pybedtools.cleanup
    pybedtools.debug_mode


Working with "chromsizes" or assembly coordinate files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Many `BEDTools` programs need "genome files" or "chromsizes" files so as to
remain within the coordinates of the assembly you're working on.  These
functions help manage these files.

.. autosummary::
    :toctree: autodocs

    pybedtools.get_chromsizes_from_ucsc
    pybedtools.chromsizes
    pybedtools.chromsizes_to_file


Performing operations in parallel (multiprocessing)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
    :toctree: autodocs

    pybedtools.parallel.parallel_apply

:mod:`pybedtools.contrib`
-------------------------
The :mod:`pybedtools.contrib` module contains higher-level code that leverages
:class:`BedTool` objects for common analyses.


Plotting
~~~~~~~~
Plotting results from BEDTools/pybedtools operations is very useful for
exploring and understanding the tools as well as for teaching purposes.

.. autosummary::
    :toctree: autodocs

    pybedtools.contrib.plotting.Track
    pybedtools.contrib.plotting.TrackCollection
    pybedtools.contrib.plotting.binary_heatmap
    pybedtools.contrib.plotting.binary_summary
    pybedtools.contrib.plotting.BedToolsDemo
    pybedtools.contrib.plotting.ConfiguredBedToolsDemo




Working with bigWig files
~~~~~~~~~~~~~~~~~~~~~~~~~
At this time, :mod:`pybedtools` does not support reading bigWig files, only
creating them via UCSC utilities.

.. autosummary::
    :toctree: autodocs

    pybedtools.contrib.bigwig.bam_to_bigwig
    pybedtools.contrib.bigwig.bedgraph_to_bigwig
    pybedtools.contrib.bigwig.wig_to_bigwig

Working with bigBed files
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
    :toctree: autodocs

    pybedtools.contrib.bigbed.bigbed
    pybedtools.contrib.bigbed.bigbed_to_bed

:class:`MultiClassifier`
~~~~~~~~~~~~~~~~~~~~~~~~
An example use-case of the :class:`MultiClassifier` class would be to determine the
distribution of ChIP-seq peaks in introns/exons/intergenic space.

.. autosummary::
    :toctree: autodocs

    pybedtools.contrib.MultiClassifier
    pybedtools.contrib.MultiClassifier.classify
    pybedtools.contrib.MultiClassifier.print_table

:class:`IntersectionMatrix`
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The :class:`IntersectionMatrix` class makes it easy to intersect a large number
of interval files with each other.

.. autosummary::
    :toctree: autodocs

    pybedtools.contrib.IntersectionMatrix

:mod:`contrib.venn_maker`
~~~~~~~~~~~~~~~~~~~~~~~~~
The :mod:`venn_maker` module helps you make Venn diagrams using the R package
`VennDiagram <http://www.biomedcentral.com/1471-2105/12/35>`_.

Note that Venn diagrams are not good for when you have nested intersections.
See the docs for :func:`pybedtools.contrib.venn_maker.cleaned_intersect` and
its source for more details.

.. autosummary::
    :toctree: autodocs

    pybedtools.contrib.venn_maker
    pybedtools.contrib.venn_maker.venn_maker
    pybedtools.contrib.venn_maker.cleaned_intersect

Scripts
-------
These scripts demonstrate ways of using :mod:`pybedtools` for genomic analyses.

Typically a script will be added here and if the functionality is useful, it is
abstracted out into a more powerful and flexible module.  For example, the
:mod:`pybedtools.contrib.venn_maker` module is a more powerful and flexible way
of making Venn diagrams than the simpler `venn_mpl` and `venn_gchart` scripts
below.

Another example is the :mod:`pybedtools.contrib.IntersectionMatrix` class,
which extends the `intersection_matrix.py` script.  The class stores results
and timestamps in a local sqlite3 database to avoid re-computing up-to-date
results.

.. autosummary::
    :toctree: autodocs

    pybedtools.scripts.pybedtools_demo
    pybedtools.scripts.venn_mpl
    pybedtools.scripts.venn_gchart
    pybedtools.scripts.intersection_matrix
    pybedtools.scripts.peak_pie
    pybedtools.scripts.annotate
    pybedtools.scripts.intron_exon_reads
    pybedtools.scripts.py_ms_example
