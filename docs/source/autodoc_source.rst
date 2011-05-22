.. include:: includeme.rst

.. testsetup:: *
    :options: +NORMALIZE_WHITESPACE

    import pybedtools
    from pybedtools import chromsizes

.. module:: pybedtools

.. _autodoc:

:mod:`pybedtools` Reference
+++++++++++++++++++++++++++

This section is the module reference documentation, and includes the full
docstrings for methods and functions in :mod:`pybedtools`.  It is separated
into :ref:`module-level`, :ref:`wrappers`, and :ref:`pbt-unique`.

.. contents::


.. _module-level:

:mod:`pybedtools` module-level functions
========================================

Functions for working with example files
----------------------------------------

.. automodule:: pybedtools
    :members: example_bedtool, example_filename, list_example_files,
              data_dir

Functions for specifying genome assemblies
------------------------------------------

.. automodule:: pybedtools
    :members: chromsizes, chromsizes_to_file, get_chromsizes_from_ucsc

Setup
-----
.. automodule:: pybedtools
    :members: set_tempdir, get_tempdir, set_bedtools_path

Utilities
---------
.. automodule:: pybedtools
    :members: cleanup, IntervalIterator, find_tagged


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

:meth:`flank` (wraps "flankBed")
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: Bedtool.flank

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



