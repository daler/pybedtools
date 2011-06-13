.. include:: includeme.rst

Create a :class:`BedTool`
-------------------------
First, follow the :ref:`installation` instructions if you haven't already
done so to install both BEDTools_ and :mod:`pybedtools`.

Then import the :mod:`pybedtools` module and make a new :class:`BedTool`::

    >>> import pybedtools

    >>> # use a BED file that ships with pybedtools...
    >>> a = pybedtools.example_bedtool('a.bed')

    >>> # ...or use your own by passing a filename
    >>> a = pybedtools.BedTool('peaks.bed')

This documentation uses example files that ship with :mod:`pybedtools`.  To
access these files from their installation location, we use the
:func:`example_bedtool` function.  This is convenient because if you
copy-paste the examples, they will work. If you would rather learn using
your own files, just pass the filename to a new :class:`BedTool`, like the
above example.

You can use any file that BEDTools_ supports -- this includes BED, VCF,
GFF, and gzipped versions of any of these. See :ref:`Creating a BedTool`
for more on the different ways of creating a :class:`BedTool`, including
from iterators and directly from a string.

Now, let's see how to do a common task performed on BED files: intersections.
