.. testsetup:: *
    :options: +NORMALIZE_WHITESPACE

    import pybedtools
    from pybedtools import chromsizes


Module documentation
====================

:mod:`pybedtools` module-level functions
----------------------------------------

.. automodule:: pybedtools
    :members:

:class:`BedTool` methods that wrap BEDTools programs
----------------------------------------------------
The following methods wrap BEDTools_ programs.  This package is still in
development; the goal is to eventually support all BEDTools_ programs.

:class:`BedTool.intersect`
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.intersect

:class:`BedTool.merge`
~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.merge

:class:`BedTool.subtract`
~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.subtract

:class:`BedTool.sequence`
~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.sequence

:class:`BedTool.closest`
~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.closest

:class:`BedTool.window`
~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.window

:class:`BedTool.sort`
~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.sort

:class:`BedTool.slop`
~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.slop

:class:`BedTool.shuffle`
~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.shuffle

:class:`BedTool` methods unique to :mod:`pybedtools`
----------------------------------------------------
The following methods are currently only supported for use with BED format
files; support for other file types is under development.

:class:`BedTool.count`
~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.count

:class:`BedTool.saveas`
~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.saveas

:class:`BedTool.lengths`
~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.lengths

:class:`BedTool.features`
~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.features

:class:`BedTool.print_sequence`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.print_sequence

:class:`BedTool.save_seqs`
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.save_seqs

:class:`BedTool.cat`
~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.cat

:class:`BedTool.tostring`
~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.tostring

:class:`BedTool.sequence_coverage`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.sequence_coverage

:class:`BedTool.counts`
~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.counts

:class:`BedTool.normalized_counts`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.normalized_counts

:class:`BedTool.delete_temporary_history`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.delete_temporary_history
