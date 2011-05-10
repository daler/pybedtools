.. include:: includeme.rst

Create a :class:`BedTool`
-------------------------
First, follow the :ref:`installation` instructions if you haven't already
done so to install both BEDTools_ and :mod:`pybedtools`.

Then import the :mod:`pybedtools` module:

.. doctest::

    >>> import pybedtools

This documentation will use example files that ship with :mod:`pybedtools`,
but you may find it more useful to use your own files.  To create a
:class:`BedTool` you need to specify a filename.  To get a list of example
files that ship with :mod:`pybedtools`:

.. doctest::

    >>> pybedtools.list_example_files()
    ['a.bed', 'b.bed', 'c.gff', 'd.gff', 'hg19-genes.bed.gz', 'rmsk.hg18.chr21.small.bed', 'rmsk.hg18.chr21.small.bed.gz']

Note that there are BED and GFF files, some of which are compressed.  All
files supported by BEDTools_ are supported by :mod:`pybedtools`.

Using one of these filenames, get the full path to the file (which will
depend on your operating system and how you installed the package) with:

.. doctest::

    >>> full_path = pybedtools.example_filename('a.bed')

Once you have a filename (whether it's an example file or your own file),
creating a :class:`BedTool` is easy:

.. doctest::

    >>> # create a new BedTool using that filename
    >>> a = pybedtools.BedTool(full_path)

Let's set up a second :class:`BedTool` so we can do intersections and
subtractions.  This time, we'll use a convenience function for creating
:class:`BedTool` instances out of example files (if you're using your own
files, just make another one the same way `a` was made above).

.. doctest::

    >>> # create another BedTool to play around with
    >>> b = pybedtools.example_bedtool('b.bed')

See :ref:`Creating a BedTool` for more information, including making
:class:`BedTool` objects directly from strings and iterators.
