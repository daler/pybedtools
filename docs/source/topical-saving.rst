.. include:: includeme.rst

Saving :class:`BedTool` results
===============================
In general, there are three different ways of saving results from
:class:`BedTool` operations:

Use the :meth:`BedTool.saveas` method
-------------------------------------
The :meth:`BedTool.saveas` method makes a **copy** of the results, so beware
that for large files, this can be time and/or memory-consuming.  However, when
working with a streaming or iterating :class:`BedTool`, this is a great way to
render the results to disk in the middle of a pipeline.

A good example of this is saving the results from a :meth:`BedTool.each` call:


.. doctest::

    >>> from pybedtools.featurefuncs import TSS
    >>> a = pybedtools.example_bedtool('a.bed')
    >>> result = a.each(TSS, upstream=1000, downstream=0)\
    ...     .saveas('upstream_regions.bed')

Use the :meth:`BedTool.moveto` method
-------------------------------------
The :meth:`BedTool.moveto` method does a **move** operation of the results.
This is best used when the results have been written to disk already (perhaps
to a tempfile) but you'd like to give the file a more reasonable/memorable
name.

.. doctest::

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')
    >>> c = a.intersect(b).moveto('intersection_of_a_and_b.bed')


Use the ``output`` keyword argument
-----------------------------------
If you know ahead of time that you want to save the output to a particular
file, use the ``output`` keyword argument to any wrapped :class:`BedTool`
method that returns another :class:`BedTool` object.  This will override the
default behavior of creating a tempfile.

.. doctest::


    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')
    >>> c = a.intersect(b, output='intersection_of_a_and_b.bed')

