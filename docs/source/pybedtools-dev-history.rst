.. _devmodel:

:mod:`pybedtools` development model
===================================
:mod:`pybedtools` is very much an open-source project. We do all of our
development in a public github repository
(https://github.com/daler/pybedtools). Initially Ryan Dale created `pybedtools`
as a wrapper for the BEDTools command-line that allowed whole-file operations
(e.g., intersecting two BED files).  At around the same time, Aaron Quinlan
began `bedtools-python`, a Cython wrapper to the BEDTools C++ API which allowed
per-line operations. After using both libraries, Brent Pedersen made an initial
attempt to merge the two libraries so that one could do a whole-file operation
and then iterate line-wise over the result.

All three authors -- especially Ryan -- then worked on the integration and further
improvements. We often discussed individual commits and design decisions using
the github interface. These discussions (often visible in github tickets such
as https://github.com/daler/pybedtools/issues/14) facilitated the continued
collaboration, and our daily use of the library have shaped pybedtools into
what it is today.

As three independent bioinformaticians who have not previously worked together,
using github as a place for discussing design decisions and coding standards
has been invaluable.
