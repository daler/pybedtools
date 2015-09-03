Notes on BAM file semantics
---------------------------
These are some implementation notes about how BAM files are handled by
mod:`pybedtools` for those interested in the implementation.

The initial creation of a :class:`BedTool` that points to a file will
trigger a check on the first 15 bytes of a file to see if it's a BAM file.
If so, then the BedTool's `_isbam` attribute is set to `True`.  If the
:class:`BedTool` is a stream, then the check will not be made, and it is up
to the creator (whether it's the user on the command line or a method or
function) to set the BAM-streaming BedTool's `._isbam` attribute to `True`.
This is handled automatically for wrapped BEDTools programs (described
below).

Some BEDTools programs natively handle BAM files.  The `@_wraps` decorator
that is used to wrap each method has a `bam` kwarg that specifies what
input argument the wrapped tool will accept as BAM (for example, the
wrapper for `intersectBed` has the kwarg `bam="abam"`).

If `self._isbam == True`, then `self.fn` is passed to the `bam` input arg
instead of the default implicit input arg (so `intersectBed`, `self.fn` is
passed as `abam` instead of `-a`).

Trying to call a method that does not have a `bam` kwarg registered will
result in a ValueError, along with a message that says to use
:meth:`BedTool.bam_to_bed()` first.  For example, `subtractBed` currently
doesn't accept BAM files as input, so this doesn't work::

    >>> a = pybedtools.example_bedtool('gdc.bam')
    >>> b = pybedtools.example_bedtool('gdc.gff')

    >>> # doesn't work:
    >>> c = a.subtract(b)

However, converting to `a` to BED format first (and setting `stream=True`
to save on disk I/O) works fine::

    >>> # works:
    >>> c = a.bam_to_bed(stream=True).subtract(b)

Iterating over a file-based BedTool that points to a BAM will call
`samtools view` and yields lines which sent to `IntervalIterator`, which
splits the lines and passes them to `create_interval_from_list` which in
turn decides on the fly whether it's gff, bed, or sam.

However, we can't easily check the first 15 bytes of a streaming BedTool,
because that would consume those bytes. The `@_wraps` decorator needs to
know some information about which arguments to a wrapped program result in
BAM output and which result in non-BAM output.

Given `a = BedTool('x.bam')`:

* `c = a.intersect(b)` creates BAM output, so it returns a new BedTool with
  `c._isbam = True`.

* `a.intersect(b, bed=True)` returns BED output.  `@_wraps` needs to know, if the
  input was BAM, which kwarg[s] disable BAM output. For example, if `-bed`
  is passed to `intersectBed`, the output will NOT be BAM.  This is
  implemented with the `nonbam` kwarg for :func:`_wraps`.  In this case,
  the resulting BED file is treated like any other BED file.

* `c = a.intersect(b, stream=True)` returns streaming BAM output.  In this
  case, iterating over `c` will send the BAM stream to stdin of a samtools
  call
