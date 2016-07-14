.. include:: includeme.rst

Changelog
=========
Changes in v0.7.8
-----------------
* Be more careful about BAM vs bgzipped files (#168)
* `BedTool.bgzip` now preserves the header when sorting
* In Python 3, parsed BEDTools help string is decoded properly
* Ensure integer number of processes in Python 3 (thanks Illa Shamovsky)
* Add details on IOError messages for broken pipe error
* Make converting to pandas.DataFrames easier with non-standard BED files (thanks Panos Firmpas)

Changes in v0.7.7
-----------------
* Chromsizes for dm6 and mm10 assemblies added to `genome_registry`
* Better Python 3 compatibility in the `long_range_interaction` module
* New `featurefuncs.UniqueID` class, useful for ensuring all features in a file
  have a unique ID in their name field.
* Fix error message when a specified genome file doesn't exist (thanks Saket Choudhary)

Changes in v0.7.6
-----------------
* New module `pybedtools.contrib.long_range_interaction` for working with
  HiC-like data.

Changes in v0.7.5
-----------------
* When using tabix-indexed files, `tabix` and `bgzip` are no longer required to
  be installed separately. Only `pysam` is needed.

* Recent BEDTools releases support multiple files for the `-b` argument of
  `bedtools intersect`. This version of `pybedtools` now supports multiple
  files as well. Note that it is still possible to provide a list of strings
  representing intervals as the `b` argument to `BedTool.intersect`. To
  differentiate between a list of intervals and a list of filenames, the first
  item converted into an `Interval` object; if it fails then consider the items
  to be filenames; otherwise assume strings representing intervals. This check
  only occurs if the `b` argument is a list or tuple; other iterable types are
  always assumed to be intervals.

Changes in v0.7.4
-----------------
Bugfix release.

- fix `#147 <https://github.com/daler/pybedtools/issues/147>`_ so that warnings
  are simply passed to the user without raising exceptions
- in setup.py, allow depedencies to have "loose" versions with suffixes like
  "rc1"
- fix in `BedTool.cat()` on empty files (thanks Brad Chapman (`PR #149
  <https://github.com/daler/pybedtools/pull/149>`_)

Changes in v0.7.1
-----------------
This is largely a bugfix release with the following changes:

- fix for some BAM headers (thanks Gabriel Platt)
- unified IntervalIterator to address some streaming issues (fixes #143)
- fix bug where `__add__` was not re-raising exceptions (thanks Brad Chapman
  and Dan Halligan)


Changes in v0.7.0
-----------------
This release reflects a major upgrade in the underlying code in order to
support both Python 2 and Python 3 using the same code. Aside from trivial
things like converting print statements to functions and using `next()` instead
of `.next()`, this required a substantial rewrite to support the way strings
are handled in Python 3 (in Cython and wrapped C++) and how relative modules
work.

Importantly, after converting them to Python 2- and 3-compatible syntax *all
previous tests pass* so to the end user should not notice any differences
except those noted below.

Strings from Interval fields are unicode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For consistency between Python 2 and 3, all strings from Interval objects are
now unicode. That is, in Python 2, previously we would get this::

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> a[0].name
    'feature1'

Now, we get this::

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> a[0].name
    u'feature1'


samtools no longer a dependency
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The dependency for samtools has been removed, which simplifies the installation
process. Instead, `pysam` is used for handling BAM files.

In order for existing tests to pass, `pysam.AlignedSegment` objects are
currently converted to `pybedtools.Interval` objects when iterating over a BAM
file. This will come at a performance cost if you are iterating over all reads
in a BAM file using the `pybedtools.BAM` object.

In the future, iterating over a BAM file will yield `pysam.AlignedSegment`
objects directly, but for now you can use the `pybedtools.BAM.pysam_bamfile`
attribute to access the underlying `pysam.AlignmentFile`

Cython no longer a dependency
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The Cythonized ``.cxx`` files are now shipped with the `pybedtools`
distribution, so Cython is no longer a requirement for installation.

You will however need to have Cython installed if you're developing pybedtools.

Remote BAM support clarification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Previously, `pybedtools` was able to support remote BAMs by loosely wrapping
samtools, but BAM files still needed to be fully downloaded to disk before
using with BEDTools. This was done automatically, but through an inefficient
mechanism.

Pysam does support remote BAMs, and as before, a BAM file needs to be created
on disk for use with BEDTools. But now this needs to be explicitly done by the
user, which should result in better performance.


Iterating over intervals
~~~~~~~~~~~~~~~~~~~~~~~~
Previously, when iterating over a `BedTool` object, different machinery would
be invoked depending on whether the BedTool was pointing to a file (a
cbedtools.IntervalFile would be invoked), to another iterator of Interval
objects, or to a stream like from the stdout of a BEDTools call
(cbedtools.IntervalIterator in both cases).

Everything is now an IntervalIterator, simplifying the path towards
performance optimization.

gzip support
~~~~~~~~~~~~
Thanks to Saulius Lukauskas, gzip handling is now improved, and calling
`BedTool.saveas()` with a `.gz` extension on the filename will automatically
compress the output.

Docker
~~~~~~
In the github repo there is a `docker` directory containing Dockerfiles to set
up isolated testing environments. These Dockerfiles also demonstrate how to set
up a complete environment starting from a base Ubuntu install.

Tests
~~~~~
All tests from v0.6.9 (which was Python 2 only) have been made Python 2/3
compatible and all previous tests pass.

If you have docker installed, from the top level directory, you can run the
full tests like this::

    cd docker
    ./full-tests.sh

This will build docker containers for Python 2 and Python 3 with all
depedencies, export the parent directory to the container, and run the test
suite.


Conda packages
~~~~~~~~~~~~~~
You can now install the latest versions of tabix, bedtools, pysam, and
pybedtools from conda, dramatically speeding up installation time. These
mechanisms are used for automated testing as well (see the ``condatest.sh``
script in the github repo).

To use these packages in your own environment(s), specify the `daler` conda
channel like this::

    conda install -c daler pybedtools

Note that this will not install BEDTools or tabix unless you explicitly say
so::

    conda install -c daler pybedtools bedtools tabix

.. note::

    This currently only works on Linux; contributions to Mac conda recipes (see
    the `conda` dir in the github repo) would be welcomed.

Changes in v0.6.9
-----------------
Minor bug fix release.

* improved the automatic field name handling when converting an interval file to
  a `pandas DataFrame`.
* fixed a bug in `IntervalFile` methods `all_hits`, `any_hits` and `count_hits`
  where zero-length features were being counted multiple times (thanks Brent
  Pedersen and Kyle Smith)
* bgzip and tabix paths can now be configured separately (thanks Rob Beagrie)
* fixed a bug where streaming BAM files were read fully into memory (thanks
  Alexey Sergushichev)

Changes in v0.6.8
-----------------

Bugfix: Thanks to Gabriel Pratt, `pybedtools` is no longer plagued by open filehandles
in the C code causing the notorious "Too many files open" error.

Changes in v0.6.7
-----------------
Now compatible with BEDTools v2.21.0.

The one exception is that the new `bedtools intersect` functionality that
allows multiple `-b` files is not yet implemented in `pybedtools`.

New features:

* `BedTool.fisher()` wraps the new BEDTools `fisher` tool.  The result is
  an object containing parsed results.

* `BedTool.colormap_normalize()` accepts a `percentile` argument, useful when
  applying colormaps to data with a handful of extreme outliers

* `BedTool.to_datafame()` converts a `BedTool` object into a `pandas.DataFrame`
  with columns named after the appropriate fields for the filetype (thanks
  Radhouane Aniba for the suggestion)

* `BedTool.tail()` to complement `BedTool.head()` (thanks Radhouane Aniba for
  the suggestion)

* Add hg38 and hg38.default chromsizes

Minor bug fixes:

* Ensure tuple-like args to `parallel_apply` (fixes #109)

* Temp fix for BEDTools v2.20.0 which required the `-w` arg to come before the
  `-s` arg in `bedtools makewindows` (#81)

* Better (i.e., UCSC Genome Browser-compliant) defaults for `featurefuncs.expand_fields`.

* Fix for BedTool.all_hits() and any_hits() which will now show hits for
  zero-length features intersecting with other zero-length features with the
  same coordinates.



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
