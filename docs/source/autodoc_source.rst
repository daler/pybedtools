
.. _autodoc:

.. _pybedtools reference:

.. currentmodule:: pybedtools

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

    pybedtools.bedtool.BedTool

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


.. autosummary::
    :toctree: autodocs

    pybedtools.bedtool.BedTool.intersect
    pybedtools.bedtool.BedTool.window
    pybedtools.bedtool.BedTool.closest
    pybedtools.bedtool.BedTool.coverage
    pybedtools.bedtool.BedTool.map
    pybedtools.bedtool.BedTool.genome_coverage
    pybedtools.bedtool.BedTool.merge
    pybedtools.bedtool.BedTool.cluster
    pybedtools.bedtool.BedTool.complement
    pybedtools.bedtool.BedTool.subtract
    pybedtools.bedtool.BedTool.slop
    pybedtools.bedtool.BedTool.flank
    pybedtools.bedtool.BedTool.sort
    pybedtools.bedtool.BedTool.random
    pybedtools.bedtool.BedTool.shuffle
    pybedtools.bedtool.BedTool.annotate
    pybedtools.bedtool.BedTool.multi_intersect
    pybedtools.bedtool.BedTool.union_bedgraphs
    pybedtools.bedtool.BedTool.pair_to_bed
    pybedtools.bedtool.BedTool.pair_to_pair
    pybedtools.bedtool.BedTool.bam_to_bed
    pybedtools.bedtool.BedTool.to_bam
    pybedtools.bedtool.BedTool.bedpe_to_bam
    pybedtools.bedtool.BedTool.bed6
    pybedtools.bedtool.BedTool.bam_to_fastq
    pybedtools.bedtool.BedTool.sequence
    pybedtools.bedtool.BedTool.mask_fasta
    pybedtools.bedtool.BedTool.nucleotide_content
    pybedtools.bedtool.BedTool.multi_bam_coverage
    pybedtools.bedtool.BedTool.tag_bam
    pybedtools.bedtool.BedTool.jaccard
    pybedtools.bedtool.BedTool.reldist
    pybedtools.bedtool.BedTool.overlap
    pybedtools.bedtool.BedTool.links
    pybedtools.bedtool.BedTool.igv
    pybedtools.bedtool.BedTool.window_maker
    pybedtools.bedtool.BedTool.groupby
    pybedtools.bedtool.BedTool.expand

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

    pybedtools.bedtool.BedTool.each
    pybedtools.bedtool.BedTool.filter
    pybedtools.bedtool.BedTool.split
    pybedtools.bedtool.BedTool.truncate_to_chrom
    pybedtools.bedtool.BedTool.remove_invalid

The :mod:`pybedtools.featurefuncs` module contains some commonly-used functions
that can be passed to :meth:`BedTool.each`:

.. currentmodule:: pybedtools

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

    pybedtools.bedtool.BedTool.all_hits
    pybedtools.bedtool.BedTool.any_hits
    pybedtools.bedtool.BedTool.count_hits
    pybedtools.bedtool.BedTool.tabix_intervals
    pybedtools.bedtool.BedTool.tabix
    pybedtools.bedtool.BedTool.bgzip


:class:`BedTool` introspection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
These methods provide information on the :class:`BedTool` object.

If using :meth:`BedTool.head`, don't forget that you can index into
:class:`BedTool` objects, too.

.. autosummary::
    :toctree: autodocs

    pybedtools.bedtool.BedTool.head
    pybedtools.bedtool.BedTool.count
    pybedtools.bedtool.BedTool.field_count
    pybedtools.bedtool.BedTool.file_type


Randomization helpers
~~~~~~~~~~~~~~~~~~~~~
Helper methods useful for assessing empirical instersection
distributions between interval files.

.. autosummary::
    :toctree: autodocs

    pybedtools.bedtool.BedTool.parallel_apply
    pybedtools.bedtool.BedTool.randomstats
    pybedtools.bedtool.BedTool.randomintersection
    pybedtools.bedtool.BedTool.randomintersection_bp
    pybedtools.bedtool.BedTool.random_subset
    pybedtools.bedtool.BedTool.random_jaccard
    pybedtools.bedtool.BedTool.random_op

Managing :class:`BedTool` objects on disk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
These methods are used to specify where to save results from :class:`BedTool`
operations.

.. autosummary::
    :toctree: autodocs

    pybedtools.bedtool.BedTool.saveas
    pybedtools.bedtool.BedTool.moveto


Misc operations
~~~~~~~~~~~~~~~
Methods that can't quite be categorized into the above sections.

.. autosummary::
    :toctree: autodocs

    pybedtools.bedtool.BedTool.cat
    pybedtools.bedtool.BedTool.at
    pybedtools.bedtool.BedTool.absolute_distance
    pybedtools.bedtool.BedTool.cut
    pybedtools.bedtool.BedTool.total_coverage
    pybedtools.bedtool.BedTool.with_attrs
    pybedtools.bedtool.BedTool.as_intervalfile
    pybedtools.bedtool.BedTool.introns
    pybedtools.bedtool.BedTool.set_chromsizes
    pybedtools.bedtool.BedTool.print_sequence
    pybedtools.bedtool.BedTool.save_seqs
    pybedtools.bedtool.BedTool.seq
    pybedtools.bedtool.BedTool.liftover
    pybedtools.bedtool.BedTool.colormap_normalize
    pybedtools.bedtool.BedTool.relative_distance

Module-level functions
----------------------

Working with example files
~~~~~~~~~~~~~~~~~~~~~~~~~~
:mod:`pybedtools` comes with many example files.  Here are some useful
functions for accessing them.

.. autosummary::
    :toctree: autodocs

    pybedtools.bedtool.example_bedtool
    pybedtools.filenames.list_example_files
    pybedtools.filenames.example_filename

Creating :class:`Interval` objects from scratch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:class:`Interval` objects are the core object in :mod:`pybedtools` to represent
a genomic interval, written in Cython for speed.

.. autosummary::
    :toctree: autodocs

    pybedtools.cbedtools.Interval
    pybedtools.cbedtools.create_interval_from_list

:mod:`pybedtools` setup and config
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Use these functions right after importing in order to use custom paths or to
clean up the temp directory.

.. autosummary::
    :toctree: autodocs

    pybedtools.helpers.set_bedtools_path
    pybedtools.helpers.get_tempdir
    pybedtools.helpers.set_tempdir
    pybedtools.helpers.cleanup
    pybedtools.debug_mode


Working with "chromsizes" or assembly coordinate files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Many `BEDTools` programs need "genome files" or "chromsizes" files so as to
remain within the coordinates of the assembly you're working on.  These
functions help manage these files.

.. autosummary::
    :toctree: autodocs

    pybedtools.helpers.get_chromsizes_from_ucsc
    pybedtools.helpers.chromsizes
    pybedtools.helpers.chromsizes_to_file


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

:mod:`contrib.long_range_interaction`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
    :toctree: autodocs

    pybedtools.contrib.long_range_interaction.tag_bedpe
    pybedtools.contrib.long_range_interaction.cis_trans_interactions
