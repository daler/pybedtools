Module documentation
====================

:mod:`pybedtools` module-level functions
----------------------------------------

.. automodule:: pybedtools
    :members:

:class:`bedtool` methods that wrap BEDTools programs
----------------------------------------------------
The following methods wrap BEDTools_ programs.  This package is still in
development; the goal is to eventually support all BEDTools_ programs.

:class:`bedtool.intersect`
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.intersect

:class:`bedtool.merge`
~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.merge

:class:`bedtool.subtract`
~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.subtract

:class:`bedtool.sequence`
~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.sequence

:class:`bedtool.closest`
~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.closest

:class:`bedtool.window`
~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.window

:class:`bedtool.sort`
~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.sort

:class:`bedtool.slop`
~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.slop

:class:`bedtool.shuffle`
~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.shuffle

:class:`bedtool` methods unique to :mod:`pybedtools`
----------------------------------------------------
The following methods are currently only supported for use with BED format
files; support for other file types is under development.

:class:`bedtool.count`
~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.count

:class:`bedtool.saveas`
~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.saveas

:class:`bedtool.get_chromsizes_from_ucsc`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.get_chromsizes_from_ucsc

:class:`bedtool.size_filter`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.size_filter

:class:`bedtool.lengths`
~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.lengths

:class:`bedtool.features`
~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.features

:class:`bedtool.print_sequence`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.print_sequence

:class:`bedtool.save_seqs`
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.save_seqs

:class:`bedtool.cat`
~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.cat

:class:`bedtool.tostring`
~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.tostring

:class:`bedtool.sequence_coverage`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.sequence_coverage

:class:`bedtool.counts`
~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.counts

:class:`bedtool.normalized_counts`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.normalized_counts

:class:`bedtool.delete_temporary_history`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: bedtool.delete_temporary_history
