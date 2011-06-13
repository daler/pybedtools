notes on BAM file semantics
---------------------------

Creating a BedTool triggers a check on the first 15 bytes of a file to see
if it's a BAM file.  If so, then `self._isbam = True`.

If `self._isbam == True`, then substitute `self.fn` for `abam` instead of
`a` (or `ibam` instead of `i`) for BEDTools programs that accept BAM input.
This is specified by the `bam` kwarg given to @_wraps.

Iterating over a file-based BedTool that points to a BAM will call
`samtools view` and yields lines which sent to `IntervalIterator`. In
`IntervalIterator`, these lines are split and passed to
`create_interval_from_list` which in turn decides on the fly whether it's
gff, bed, or sam.


Given `a = BedTool('x.bam')`:

* `c = a.intersect(b)` creates BAM output, so it returns a new BedTool with
  `c._isbam=True`.  This will allow the `__iter__` method to know to use
  `IntervalIterator(BAM(c.fn))` when iterating 

* `a.intersect(b, bed=True)` returns BED output.  `@_wraps` needs to know, if the
  input was BAM, which kwarg disables BAM output, e.g., `-bed` for
  `intersectBed`.  This is implemented with the `nonbam` kwarg for
  :func:`_wraps`.  In this case, the resulting BED file is treated like any
  other BED file.


* `c = a.intersect(b, stream=True)` returns streaming BAM output.  In this
  case, BAM needs to send the stream to stdin of the samtools call when
  iterating over `c`
