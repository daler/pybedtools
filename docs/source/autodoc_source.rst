.. currentmodule:: pybedtools
.. contents::

.. _autodoc:

:mod:`pybedtools` Reference
+++++++++++++++++++++++++++

This section is the module reference documentation, and includes the full
docstrings for methods and functions in :mod:`pybedtools`.  It is separated
into :ref:`module-level`, :ref:`wrappers`, and :ref:`pbt-unique`.




.. autoclass:: pybedtools.BedTool
    :members: __init__

.. _module-level:

:mod:`pybedtools` module-level functions
========================================

Functions for working with example files
----------------------------------------

.. autofunction:: pybedtools.example_bedtool
.. autofunction:: pybedtools.example_filename
.. autofunction:: pybedtools.list_example_files
.. autofunction:: pybedtools.data_dir

Functions for specifying genome assemblies
------------------------------------------

.. autofunction:: pybedtools.chromsizes
.. autofunction:: pybedtools.chromsizes_to_file
.. autofunction:: pybedtools.get_chromsizes_from_ucsc

Setup
-----

.. autofunction:: pybedtools.set_tempdir
.. autofunction:: pybedtools.get_tempdir
.. autofunction:: pybedtools.set_bedtools_path

Utilities
---------
.. autofunction:: pybedtools.cleanup
.. autofunction:: pybedtools.IntervalIterator
.. autofunction:: pybedtools.find_tagged

Wrapping
--------
.. autofunction:: pybedtools.bedtool._wraps


.. _wrappers:

:class:`BedTool` methods that wrap BEDTools programs
====================================================

"Genome algebra" methods
------------------------

:meth:`intersect` (wraps "intersectBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.intersect

:meth:`merge` (wraps "mergeBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.merge

:meth:`subtract` (wraps "subtractBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.subtract

:meth:`closest` (wraps "closestBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.closest

:meth:`window` (wraps "windowBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.window

:meth:`sort` (wraps "sortBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.sort

:meth:`slop` (wraps "slopBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.slop

:meth:`complementBed` (wraps "complementBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.complement

:meth:`flank` (wraps "flankBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.flank

:meth:`shuffle` (wraps "shuffleBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.shuffle

:meth:`annotate` (wraps "annotateBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.annotate

:meth:`coverage` (wraps "coverageBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.coverage

:meth:`genome_coverage` (wraps "genomeCoverageBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.genome_coverage

:meth:`overlap` (wraps "overlap")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.overlap

:meth:`groupby` (wraps "groupBy")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.groupby

:meth:`pair_to_bed` (wraps "pairToBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.pair_to_bed

:meth:`pair_to_pair` (wraps "pairToPair")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.pair_to_pair



Methods for converting between formats
--------------------------------------

:meth:`bed6` (wraps "Bed12To6")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.bed6




Methods for working with sequences
----------------------------------

:meth:`sequence` (wraps "fastaFromBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.sequence

:meth:`mask_fasta` (wraps "maskFastaFromBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.mask_fasta




.. _pbt-unique:

:class:`BedTool` methods unique to :mod:`pybedtools`
====================================================



Introspection
-------------

:meth:`count`
~~~~~~~~~~~~~
.. automethod:: BedTool.count

:meth:`print_sequence`
~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.print_sequence

:meth:`field_count`
~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.field_count

:meth:`head`
~~~~~~~~~~~~
.. automethod:: BedTool.head



Saving
------

:meth:`saveas`
~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.saveas

:meth:`save_seqs`
~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.save_seqs





Utilities
---------

:meth:`with_attrs`
~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.with_attrs

:meth:`cat`
~~~~~~~~~~~
.. automethod:: BedTool.cat

:meth:`total_coverage`
~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.total_coverage

:meth:`delete_temporary_history`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.delete_temporary_history

:meth:`as_intervalfile`
~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.as_intervalfile





Feature-by-feature operations
-----------------------------

:meth:`each`
~~~~~~~~~~~~
.. automethod:: BedTool.each

:meth:`filter`
~~~~~~~~~~~~~~
.. automethod:: BedTool.filter

:meth:`cut`
~~~~~~~~~~~
.. automethod:: BedTool.cut

:meth:`features`
~~~~~~~~~~~~~~~~
.. automethod:: BedTool.features


Randomization helpers
---------------------

:meth:`randomintersection`
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.randomintersection

:meth:`randomstats`
~~~~~~~~~~~~~~~~~~~
.. automethod:: BedTool.randomstats


:class:`Interval` class
~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: pybedtools.Interval
    :members:

:class:`IntervalFile` class
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: pybedtools.IntervalFile
    :members:
