.. _bedtools218:

BEDTools v2.18
==============
If you use streaming BedTool objects or use the :mod:`pybedtools` methods that
perform large numbers of operations on < 10k features it's **recommended to use
BEDTools v2.17** with :mod:`pybedtools`.

BEDTools v2.18 implements some changes that result in `substantial speedups
<http://quinlanlab.org/software-releases/bedtools-2.18.html>`_ for large files,
but those same changes make parts of :mod:`pybedtools` work suboptimally.

First, BEDTools v2.18 uses output buffering.  This causes some of the system
calls in :mod:`pybedtools` to hang while waiting for output to be written to
stdout.  This occurs when using "streaming" BedTool objects (see :ref:`BedTools
as iterators`).  This is a :mod:`pybedtools` specific issue that is being
worked on and can hopefully be resolved.

Second, a new memory pre-allocation step speeds up operations for large files,
but this same step increases the overhead of each individual call.  This means
that the parts of :mod:`pybedtools` that call `bedtools intersect` many times
-- like :meth:`BedTool.parallel_apply` -- on files with < 10k features will
have a noticable slowdown due this new overhead. This effect is independent of
:mod:`pybedtools`.

These issues are being worked on and should be resolved in future releases.
