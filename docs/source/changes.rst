.. include:: includeme.rst

Changelog
=========
Changes in v0.6.6
-----------------
This is a compatibility release, updated for BEDTools v2.20.0.

There is one API change that affects the behavior of overloaded operators (that
is, using `+` and `-` with BedTool objects) when one of the BedTool objects
represents an empty file.

Assume `a` is a BedTool object representing a regular BED file but `b` is
empty.  Previously:

    * a + b = a
    * b + a = b
    * a - b = a
    * b - a = a
    * b - b = b
    * a + a = a

The following changes have been made (indicated in **bold**), which hopefully
make more logical sense:

    * **a + b = b**
    * b + a = b
    * a - b = a
    * **b - a = b**
    * b - b = b
    * a + a = a

Changes in v0.6.5
-----------------
This is a minor bug-fix release:

* Fix for BedTool.all_hits() and any_hits() which will now show hits for
  zero-length features with the same coordinates, like the rest of BEDTools.

* Improved error-handling to avoid Python interpreter crashing in cases when
  a BED file on the filesystem becomes unavailable after a BedTool object has
  been created for it.


Changes in v0.6.4
-----------------

* Full integration with BEDTools v2.18.  This includes some compatibility fixes
  for the new buffered output capabilities of BEDTool `intersect` and wrapping
  the new `bedtools sample` tool.

* Overloaded operators (`+` and `-`) allow empty files as input, even using
  BEDTools v2.18+.

* Travis-CI builds now use BEDTools v2.18+ for tests.

* Fix for :func:`pybedtools.featurefuncs.midpoint` (thanks ny-shao)

* Fix to :meth:`BedTool.randomstats` (thanks Michael Reschen)


Changes in v0.6.3
-----------------

* New :mod:`pybedtools.parallel` module for working with many operations in
  parallel.  See the docs for :func:`pybedtools.parallel.parallel_apply` for
  details.

* :func:`pybedtools.contrib.bigbed.bigbed` for converting to bigBed format,
  along with auto-SQL creation as needed.

* New function :func:`pybedtools.contrib.bigbed.bigbed_to_bed`, so now bigBed
  -> BED and BED -> bigBed interconversions are trivial.

* Support for remote BAMs by passing `remote=True` when creating
  a :class:`BedTool` object

* New method :meth:`BedTool.at` for subsetting a BedTool by a set of (sorted)
  indexes.

* New functions :func:`featurefuncs.gff2bed` and :func:`featurefuncs.bed2gff`
  for use with :meth:`BedTool.each`, for easy converting GFF/GTF to BED

* New function :func:`add_color` for applying matplotlib colormaps to BED
  files; see also new method :meth:`pybedtools.BedTool.colormap_normalize`.

* :class:`pybedtools.plotting.BinaryHeatmap` class for working with results
  from :meth:`BedTool.multi_intersect`.

* :meth:`BedTool.each` now also has some filter capabilities (if provided
  function's return value evaluates to False, feature will be skipped)

* Better detection for samtools (thanks Luca Beltrame)

* Expand BEDToolsError (thanks Ryan Layer)

* Creating a BedTool from a list of intervals now saves to temp file instead of treating
  like a consume-once iterator (#73)

* Various fixes to keyword arg handling to match semantics of BEDTools.

* Command line help and improved docs for the `peak_pie.py` script.

* Fix to GFF attributes (thanks Libor Mořkovský)

* Fix to labels in :mod:`pybedtools.contrib.venn_maker.py` (thanks Luca
  Pinello)

* Make the naive scaling (to million mapped reads) in
  :func:`pybedtools.contrib.bigwig.bam_to_bigwiq` optional.

* Fix for :meth:`BedTool.cat` to handle cases where at least one input is an
  empty file

* Removed SciPy dependency

* Every commit is built with Travis-CI for continuous integration testing of
  changes to source code.

Changes in v0.6.2
-----------------

* Wrapped new tools available in BEDTools 2.17: :meth:`BedTool.jaccard` and
  :meth:`BedTool.reldist` wrap the new `bedtools jaccard` and `bedtools
  reldist` respectively.

* Initial implementations of building blocks for computing statistics,
  :meth:`BedTool.absolute_distance` and :meth:`BedTool.relative_distance`

* :func:`pybedtools.featurefuncs.three_prime`,
  :func:`pybedtools.featurefuncs.five_prime`, and
  :func:`pybedtools.featurefuncs.TSS` modifier functions that can be passed to
  :meth:`BedTool.each`

* :func:`pybedtools.contrib.plotting.binary_heatmap` for visualizing results
  from :meth:`BedTool.multi_intersect`

* Fixed a long-standing issue where streaming :class:`BedTool` objects did not
  close their open file handles (stdout).  When working with many (i.e. tens
  of thousands) files, this caused the operating system to hit its open file
  limit.  This is now fixed.

* :meth:`BedTool.random_op`, a new mechanism for implementing operations that
  you would like to apply over tens of thousands of shuffled interval files.
  This makes it easy to extend the existing :mod:`pybedtools` multiprocessing
  functionality.

* :func:`pybedtools.contrib.bigwig.bam_to_bigwig`, a helper function to create
  a libary-size-scaled bigWig file from an input BAM file.

* :class:`pybedtools.contrib.plotting.TrackCollection` class, which handles
  plotting multiple files at once, using a provided "stylesheet" configuration
  to tweak colors etc.

* :class:`pybedtools.contrib.plotting.BedToolsDemo` and
  :class:`pybedtools.contrib.plotting.ConfiguredBedToolsDemo`, useful for
  running many graphical demos of BEDTools operations using the same
  "stylesheet" configuration.  Run :file:`pybedtools/contrib/plotting.py` for
  a demo.

* chromsizes dictionaries for common assemblies now have a `default` attribute,
  which is an OrderedDict of a default set of chromosome.  For example,
  ``pybedtools.chromsizes('hg19').default`` contains only the entries for the
  autosomes and X and Y.

* :meth:`BedTool.cat` now works better with multiprocessing

* added `include_distribution` kwarg to :meth:`BedTool.randomstats`, which will
  attach the full distribution of all the randomized files to the results
  dictionary.

* New method implementing Jaccard statistic (with pvalue using randomizations):
  :meth:`BedTool.random_jaccard`

* :func:`featurefuncs.extend_fields` helper function to pad fields with `'.'`,
  useful for manipulating features with the :meth:`BedTool.each` method

* Fixed a bug where BAM files, when written to disk via :meth:`BedTool.saveas`,
  were saved as SAM files.

* Better GTF/GFF detection, and if the input had quoted attribute values, then
  the output will, too

* various minor bug fixes and improvments as documented in the github commit
  logs....


Changes in v0.6.1
-----------------
* New :class:`pybedtools.contrib.plotting.Track` class allows plotting of
  features with matplotlib.  The `Track` class subclasses
  `matplotlib.collections.PolyCollection`, making it rather fast for 1000s of
  features.

* See the `scripts/pbt_plotting_example.py` script for a way of visually showing
  the results of BEDTools operations . . . great for teaching BEDTools to new
  users.

* New :meth:`BedTool.liftover` method (needs a chain file from UCSC and the
  `liftover` program installed)

* :class:`BedTool` creation using tuples/lists of values -- everything is
  converted to string before creating an :class:`Interval` object.

* bugfix: :meth:`BedTool.window_maker` now handles the `genome` kwarg correctly

* bugfix: `pybedtools.cleanup(remove_all=True)` now works correctly when using
  the default temp dir


Changes in v0.6
---------------
* Overhaul in online documentation to hopefully make functionality easier to
  find and/or discover.  See :ref:`pybedtools reference` for summary tables of
  the different parts of :mod:`pybedtools`; each entry is linked to further
  class/method/function-specific docs.  These more detailed docs also have
  links to view the source code from within the HTML docs for more exploration.

* :func:`pybedtools.contrib.venn_maker` function that acts as an interface to
  the VennDiagram R package -- just give it some BED files and it'll do the
  rest.

* Debug mode -- :func:`pybedtools.debug_mode` -- for verbose logging messages.

* Fixed an open file leak (OSError: too many open files) that occured when
  opening thousands of streaming bed files in a single session.

* Initial support for tabix files.  Useful for extracting features from
  a single region when you don't need a full intersection.

* New :mod:`pybedtools.contrib` module (in the spirit of Django's `contrib`)
  where higher-level functionality will be built.

* :class:`pybedtools.contrib.Classifier` class for identifying the classes of
  intervals.  Useful for making pie charts of intronic/exonic/intergenic etc
  classes of peaks.  Note that this is somewhat redundant with the new `mapBed`
  program in BEDTools.

* Experimental :class:`pybedtools.contrib.IntersectionMatrix` class for
  handling pairwise intersections of a large number of interval files --
  including a local sqlite3 database to avoid re-computing already up-to-date
  results.

* :class:`Interval` objects are now hashable (it's just a hash of the string
  representation) so that you can use them as dictionary keys.

* :meth:`BedTool.split` method, which accepts a function returning an iterable
  of :class:`Interval` objects. The function is applied to each interval.
  Useful for, say, splitting each gene into TSS, TTS, upstream and downstream
  features.

* :meth:`BedTool.truncate_to_chrom` method, which truncates features to the
  chromosome sizes of the provided genome.  Useful for when you try uploading
  a MACS-generated track to the UCSC genome browser, but it complains because
  peak boundaries have been extended outside chromosome boundaries . . . this
  method fixes the problem.

* :class:`BedTool` objects now have full functionality of :class:`IntervalFile`
  objects -- that is, they have the methods :meth:`BedTool.any_hits`,
  :meth:`BedTool.all_hits`, and :meth:`BedTool.count_hits` for doing
  single-interval tests.  Sometimes this will be faster than using the tabix
  support, sometimes it won't -- it's best to try both, depending on your data.

* String representations of :class:`Interval` objects now have a newline at the
  end, just like a raw lines from a BED/GFF/VCF file.  Previously, this was
  inconsistent and sometimes led to extra blank lines in "streaming"
  :class:`BedTool` instances . . . which in turn led to problems with BEDTools
  programs using the chromsweep algorithm.

* Concatentate multiple files with one call to :meth:`BedTool.cat` (thanks Jake
  Biesinger)

* Wrapped previous BEDTools programs:
    * `unionBedGraphs` (:meth:`BedTool.union_bedgraphs`)
    * `pairToBed` (:meth:`BedTool.pair_to_bed`)
    * `pairToPair` (:meth:`BedTool.pair_to_pair`)
    * `bedpeToBam` (:meth:`BedTool.bedpe_to_bam`)

* Wrapped new BEDTools programs:
    * `mapBed` (:meth:`BedTool.map`)
    * `clusterBed` (:meth:`BedTool.cluster`)
    * `randomBed` (:meth:`BedTool.random`)
    * `multiIntersectBed` (:meth:`BedTool.multi_intersect`)
    * `expandCols` (:meth:`BedTool.expand`)
    * `windowMaker` (:meth:`BedTool.window_maker`)
    * `bamToFastq` (:meth:`BedTool.bam_to_fastq`)

* Made venn_gchart and venn_mpl tests more stable

* Automatic documenting of which args are passed implicitly for BedTool method
  calls

* More robust mechanisms for specifying custom paths for BEDTools installation
  as well as optional tabix, samtools, and R installations.  This makes it
  easier to explicitly specify which versions of the tools to use.

* Improvements to GFF attributes: handle unescaped "=" (from sim4db GFFs) and
  make Attribute class properly dict-like (thanks Libor Mořkovský)

Changes in v0.5.5
-----------------
* Use `additional_args` kwarg to pass arguments verbatim to the underlying
  BEDTools programs.  This is necessary for arguments like
  `genomeCoverageBed`'s `-5` argument, since `5=True` is not a valid Python
  expression.  For example, you can use::

     import pybedtools
     a = pybedtools.example_bedtool('a.bed')
     a.genome_coverage(bg=True, strand='+', genome='hg19', additional_args='-5')

* Brent Pedersen added support for just 2 BED files in the Venn diagram scripts

* :meth:`BedTool.all_hits` uses the underlying BEDTools C++ API to get all hits
  in a file for a particular Interval::

    a = pybedtools.example_bedtool('a.bed')
    interval = Interval('chr1', 1, 5000)
    a.all_hits(interval)

* New semantics for comparisons of Interval objects.  Visual documentation of
  this coming soon.

* More tests for latest BEDTools code

* Interval instances are now pickleable; they can now be used across processes
  for parallel code.


Changes in v0.5
---------------
* support for running random intersections in parallel.  See
  :meth:`BedTool.randomstats` and :meth:`BedTool.randomintersection` (thanks,
  Jake Biesinger)

* Cython `Interval.__copy__()` for compatibility with `copy` module

* `seek()` and `rewind()` methods for `IntervalFile` class, used for Aaron
  Quinlan's new chromsweep algorithm (https://github.com/arq5x/chrom_sweep)
  (thanks, Aaron)

* support and tests for new BEDTools programs `multiBamCov`, `tagBam`, and `nucBed`

* `output="out.bed"` kwarg for all wrapped methods for explicitly specifying
  where to save output -- no more moving tempfiles

* docs improvements:
    * direct comparison with a shell script to illustrate benefit of
      `pybedtools`; see :ref:`shell_comparison`
    * more installation details
    * 0- and 1-based coordinates discussed early on (the 3 brief examples page,
      :ref:`3examples`)
    * development history and open collaboration model (see :ref:`devmodel`)
