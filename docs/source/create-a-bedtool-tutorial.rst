.. include:: includeme.rst

Create a :class:`BedTool`
-------------------------
First, follow the :ref:`installation` instructions if you haven't already
done so to install both BEDTools_ and :mod:`pybedtools`.

Then import the :mod:`pybedtools` module:

.. doctest::

    >>> import pybedtools

Then set up a :class:`BedTool` instance using a `BED format`_ file. This
can be a file that you already have, or one of the example files as shown
below.  Some methods currently only work for BED files -- see
:ref:`Limitations` for more info on this.  

For this tutorial, we'll use some example files that come with
:mod:`pybedtools`.  We can get the filename for the example files using the
:func:`pybedtools.example_files` function:

.. doctest::

    >>> # get an example filename to use
    >>> bed_filename_a = pybedtools.example_filename('a.bed')

The filename will depend on where you have installed :mod:`pybedtools`.
Once you have a filename, creating a :class:`BedTool` object is easy:

.. doctest::

    >>> # create a new BedTool using that filename
    >>> a = pybedtools.BedTool(bed_filename_a)

Set up a second one so we can do intersections and subtractions -- this
time, let's make a new :class:`BedTool` all in one line:

.. doctest::

    >>> # create another BedTool to play around with
    >>> b = pybedtools.BedTool(pybedtools.example_filename('b.bed'))

See :ref:`Creating a BedTool` for more information, including convenience
functions for working with example bed files and making :class:`BedTool`
objects directly from strings.
