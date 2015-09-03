.. include:: includeme.rst

Create a :class:`BedTool`
-------------------------
First, follow the :ref:`installation` instructions if you haven't already
done so to install both BEDTools_ and :mod:`pybedtools`.

Then import the :mod:`pybedtools` module and make a new :class:`BedTool`.  A
:class:`BedTool` object encapsulates all of the available BEDTools programs and
makes them easier to use within Python.  Most of the time when working with
:mod:`pybedtools` you'll be using :class:`BedTool` objects.  In general, a
single :class:`BedTool` object points to an interval file (BED, GFF, GTF, VCF,
SAM, or BAM format).

::

    >>> import pybedtools

    >>> # use a BED file that ships with pybedtools...
    >>> a = pybedtools.example_bedtool('a.bed')

    >>> # ...or use your own by passing a filename
    >>> a = pybedtools.BedTool('peaks.bed')

This documentation uses example files that ship with :mod:`pybedtools`.  To
access these files from their installation location, we use the
:func:`example_bedtool` function.  This is convenient because if you copy-paste
the examples, they will work. When using the :func:`example_bedtool` function,
the resulting :class:`BedTool` object will point to the corresponding file in
the `test/data` directory of your :mod:`pybedtools` installation.  If you would
rather learn using your own files, just pass the filename to a new
:class:`BedTool`, like the above example.

You can use any file that BEDTools_ supports -- this includes BED, VCF,
GFF, and gzipped versions of any of these. See :ref:`Creating a BedTool`
for more on the different ways of creating a :class:`BedTool`, including
from iterators and directly from a string.

Now, let's see how to do a common task performed on BED files: intersections.
