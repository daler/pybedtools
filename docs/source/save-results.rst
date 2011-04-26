.. include:: includeme.rst

Saving the results
==================

If you want to save the results as a meaningful filename for later use, use
the :meth:`BedTool.saveas` method.  This also lets you optionally specify a
trackline for directly uploading to the UCSC Genome Browser, instead of
opening up the files afterward and manually adding a trackline:

.. doctest::
   :hide:

   >>> a = pybedtools.example_bedtool('a.bed')
   >>> b = pybedtools.example_bedtool('b.bed')
   >>> a_with_b = a.intersect(b)

.. doctest::

    >>> a_with_b.saveas('intersection-of-a-and-b.bed', trackline='track name="a and b"')
    <BedTool(intersection-of-a-and-b.bed)>

Note that the :meth:`BedTool.saveas` method returns a new :class:`BedTool`
object which points to the newly created file on disk.  This allows you to
insert a :meth:`BedTool.saveas` call in the middle of a chain of commands
(described in another section below).
