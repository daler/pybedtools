.. include:: includeme.rst

Changelog
=========

Changes in v0.5
---------------
* support for running random intersections in parallel.  See :meth:`BedTool.randomstats` and
  :meth:`BedTool.randomintersection` (thanks, Jake Biesinger)
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
