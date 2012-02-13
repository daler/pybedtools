.. include:: includeme.rst

Changelog
=========
Changes since v0.5.5
--------------------
* Overhaul in documentation to hopefully make functionality easier to find
  and/or discover.  See :ref:`pybedtools reference` for summary tables of the
  different parts of :mod:`pybedtools`; each entry is linked to further
  class/method/function-specific docs.  These more detailed docs also have
  links to view the source code from within the HTML docs for more exploration.

* :func:`pybedtools.contrib.venn_maker` function that acts as an interface to
  the VennDiagram R package -- just give it some BED files and it'll do the
  rest.

* Debug mode -- :func:`pybedtools.debug_mode` -- for verbose logging messages.

* Fixed an open file leak (OSError: too many open files)  if opening thousands
  of streaming bed files in a single session.

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

* :meth:`BedTool.split` method, which accepts a function returning an
  iterable of :class:`Interval` objects. The function is applied to each
  interval.  Useful for, say, splitting each gene into TSS, TTS, upstream and
  downstream features.

* :meth:`BedTool.truncate_to_chrom` method, which truncates features to the chromosome
  sizes of the provided genome.  Useful for when you try uploading
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

* Wrapped `mapBed`

* Wrapped `multiIntersectBed`.  Because this program uses different semantics
  than other programs (e.g., does not have an implicit file to work on; -i` is
  not a single file but a list of files), you currently need to call it with
  a list of filenames as the `i` kwarg.  Future development will allow more
  flexibility, like using other BedTool objects or streaming BedTools.

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
