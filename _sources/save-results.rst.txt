.. include:: includeme.rst

.. _saveresults:

Saving the results
==================

If you want to save the results as a meaningful filename for later use, use
the :meth:`BedTool.saveas` method.  This does a copy operation on the file
pointed to by the :class:`BedTool`. This method also lets you optionally specify a
trackline for directly uploading to the UCSC Genome Browser, instead of
opening up the files afterward and manually adding a trackline:

.. doctest::
   :hide:

   >>> a = pybedtools.example_bedtool('a.bed')
   >>> b = pybedtools.example_bedtool('b.bed')
   >>> a_with_b = a.intersect(b)

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> c = a_with_b.saveas('intersection-of-a-and-b.bed', trackline='track name="a and b"')
    >>> print(c.fn)
    intersection-of-a-and-b.bed


    >>> # opening the underlying file shows the track line
    >>> print(open(c.fn).read())
    track name="a and b"
    chr1	155	200	feature2	0	+
    chr1	155	200	feature3	0	-
    chr1	900	901	feature4	0	+
    <BLANKLINE>

    >>> # printing file-based BedTool objects will not print the track line
    >>> print(c)
    chr1	155	200	feature2	0	+
    chr1	155	200	feature3	0	-
    chr1	900	901	feature4	0	+
    <BLANKLINE>

Note that the :meth:`BedTool.saveas` method returns a new :class:`BedTool`
object which points to the newly created file on disk.  This allows you to
insert a :meth:`BedTool.saveas` call in the middle of a chain of commands
(described in another section below).

Alternatively, if you do not want to add a track line, you can use the
:meth:`BedTool.moveto` method which can be much faster, especially on larger
files.  This does rename operation rather than a copy operation, which means
that trying to call :class:`BedTool` methods on the original will no longer
work because the underlying file no longer exists:

.. doctest::

    >>> d = a_with_b.moveto('another_location.bed')
