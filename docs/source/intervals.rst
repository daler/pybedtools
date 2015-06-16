.. currentmodule:: pybedtools

.. _intervals:

Intervals
=========

An :class:`Interval` object is how :mod:`pybedtools` represents a line in a BED,
GFF, GTF, or VCF file in a uniform fashion.  This section will describe
some useful features of :class:`Interval` objects.

First, let's get a :class:`BedTool` to work with:

.. doctest::

    >>> a = pybedtools.example_bedtool('a.bed')

We can access the :class:`Intervals` of `a` several different ways.
Probably the most convenient way is by indexing a :class:`BedTool` object:

.. doctest::

    >>> feature = a[0]

:class:`BedTool` objects support slices, too:

.. doctest::

    >>> features = a[1:3]

Common :class:`Interval` attributes
-----------------------------------
Printing a feature converts it into the original line from the file:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print(feature)
    chr1	1	100	feature1	0	+

The string representation of an :class:`Interval` object is simply a valid line,
**including the newline**, for the format from which that :class:`Interval` was
created (accessible via :attr:`Interval.file_type`).

All features, no matter what the file type, have `chrom`, `start`, `stop`,
`name`, `score`, and `strand` attributes.  Note that `start` and `stop` are
integers, while everything else (including `score`) is a string.

:mod:`pybedtools` supports both Python 2 and 3. When using Python 3, all
strings are the `str` type. When using Python 2, all strings are unicode.

.. note::

    This documentation undergoes testing with Python 2 and Python 3. These
    versions handle strings differently. For example, under Python 2::

        >>> feature.chrom
        u'chr1'

    But under Python 3::

        >>> feature.chrom
        'chr1'

    Since all strings returned by Interval objects are unicode, we solve this
    by making a helper function `show_value` that converts unicode to native
    string -- but only under Python 2.


.. doctest::

    >>> import sys
    >>> def show_value(s):
    ...     """
    ...     Convert unicode to str under Python 2;
    ...     all other values pass through unchanged
    ...     """
    ...     if sys.version_info.major == 2:
    ...         if isinstance(s, unicode):
    ...             return str(s)
    ...     return s

.. doctest::

    >>> show_value(feature.chrom)
    'chr1'

    >>> show_value(feature.start)
    1

    >>> show_value(feature.stop)
    100

    >>> show_value(feature.name)
    'feature1'

    >>> show_value(feature.score)
    '0'

    >>> show_value(feature.strand)
    '+'

Let's make another feature that only has chrom, start, and stop to see how
:mod:`pybedtools` deals with missing attributes:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> feature2 = pybedtools.BedTool('chrX 500 1000', from_string=True)[0]

    >>> print(feature2)
    chrX	500	1000


    >>> show_value(feature2.chrom)
    'chrX'

    >>> show_value(feature2.start)
    500

    >>> show_value(feature2.stop)
    1000

    >>> show_value(feature2.name)
    '.'

    >>> show_value(feature2.score)
    '.'

    >>> show_value(feature2.strand)
    '.'

This illustrates that default values are the string "`.`".


Indexing into :class:`Interval` objects
---------------------------------------

:class:`Interval` objects can also be indexed by position into the original
line (like a list) or indexed by name of attribute (like a dictionary).

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print(feature)
    chr1	1	100	feature1	0	+


    >>> show_value(feature[0])
    'chr1'

    >>> show_value(feature['chrom'])
    'chr1'

    >>> show_value(feature[1])
    '1'

    >>> show_value(feature['start'])
    1


Fields
------
:class:`Interval` objects have a :attr:`Interval.fields` attribute that
contains the original line split into a list of strings.  When an integer
index is used on the :class:`Interval` (for example, `feature[3]`), it is
the `fields` attribute that is actually being indexed into.

.. doctest::

    >>> f = pybedtools.BedTool('chr1 1 100 asdf 0 + a b c d', from_string=True)[0]
    >>> [show_value(i) for i in f.fields]
    ['chr1', '1', '100', 'asdf', '0', '+', 'a', 'b', 'c', 'd']
    >>> len(f.fields)
    10

.. _zero_based_coords:

BED is 0-based, others are 1-based
----------------------------------
One troublesome part about working with multiple formats is that BED files
have a different coordinate system than GFF/GTF/VCF/ files.

* **BED files are 0-based** (the first base of the chromosome is considered
  position 0) and the **feature does not include the stop position**.

* **GFF, GTF, and VCF files are 1-based** (the first base of the chromosome
  is considered position 1) and the **feature includes the stop position**.

:mod:`pybedtools` follows the following conventions:

* The value in :attr:`Interval.start` will *always contain the
  0-based start position*, even if it came from a GFF or other 1-based
  feature.

* Getting the `len()` of an :class:`Interval` will always return
  `Interval.stop - Interval.start`, so no matter what format the
  original file was in, the length will be correct.  This greatly simplifies
  underlying code, and it means you can treat all :class:`Intervals`
  identically.

* The contents of :attr:`Interval.fields` will *always be strings*,
  which in turn always represent the original line in the file.

    * This means that for a GFF feature, :attr:`Interval.fields[3]` or
      :attr:`Interval[3]`, which is 1-based according to the file format,
      will always be one bp larger than :attr:`Interval.start`, which
      always contains the 0-based start position.  Their data types
      are different; :attr:`Interval[3]` will be a string and
      :attr:`Interval.start` will be a long.

    * Printing an :class:`Interval` object created from a GFF file will
      show the tab-delimited fields in GFF coords while printing an
      :class:`Interval` object created from a BED file will show fields in
      BED coords.

Worked example
~~~~~~~~~~~~~~
To illustrate and confirm this functionality, let's create a GFF feature and
a BED feature from scratch and compare them.

First, let's create a GFF :class:`Interval` from scratch:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> gff = ["chr1",
    ...        "fake",
    ...        "mRNA",
    ...        "51",   # <- start is 1 greater than start for the BED feature below
    ...        "300",
    ...        ".", 
    ...        "+",
    ...        ".",
    ...        "ID=mRNA1;Parent=gene1;"]
    >>> gff = pybedtools.create_interval_from_list(gff)



Then let's create a corresponding BED :class:`Interval` that represents the
same genomic coordinates of of the GFF feature, but since BED format is
zero-based we need to subtract 1 from the start:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> bed = ["chr1",
    ...        "50",
    ...        "300",
    ...        "mRNA1", 
    ...        ".",
    ...        "+"]
    >>> bed = pybedtools.create_interval_from_list(bed)

Let's confirm these new features were recognized as the right file type --
the format is auto-detected based on the position of chrom/start/stop coords in
the provided field list:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> show_value(gff.file_type)
    'gff'

    >>> show_value(bed.file_type)
    'bed'

Printing the :class:`Intervals` shows that the strings are in the appropriate
coordinates:


.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> # for testing, we make sure keys are sorted. Not needed in practice.
    >>> gff.attrs.sort_keys = True
    >>> print(gff)
    chr1	fake	mRNA	51	300	.	+	.	ID=mRNA1;Parent=gene1;

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print(bed)
    chr1	50	300	mRNA1	.	+

Since `start` attributes are always zero-based, the GFF and BED `start` values
should be identical:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> bed.start == gff.start == 50
    True

For the BED feature, the second string field (representing the start position)
and the `start` attribute should both be `50` (though one is an integer and the
other is a string) . . .

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> show_value(bed.start)
    50
    >>> show_value(bed[1])
    '50'

. . . but for the GFF feature, they differ -- the `start` attribute is
zero-based while the string representation (the fourth field of a GFF file)
remains in one-based GFF coords:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> show_value(gff.start)
    50
    >>> show_value(gff[3])
    '51'

As long as we use the integer `start` attributes, we can treat the
:class:`Interval` objects identically, without having to check for their format
every time:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> len(bed) == len(gff) == 250
    True

GFF features have access to attributes
--------------------------------------
GFF and GTF files have lots of useful information in their attributes field
(the last field in each line).  These attributes can be accessed with the
:attr:`Interval.attrs` attribute, which acts like a dictionary.  For speed,
the attributes are lazy -- they are only parsed when you ask for them.  BED
files, which do not have an attributes field, will return an empty
dictionary.

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> # original feature
    >>> print(gff)
    chr1	fake	mRNA	51	300	.	+	.	ID=mRNA1;Parent=gene1;

    >>> # original attributes
    >>> sorted(gff.attrs.items())
    [('ID', 'mRNA1'), ('Parent', 'gene1')]

    >>> # add some new attributes
    >>> gff.attrs['Awesomeness'] = "99"
    >>> gff.attrs['ID'] = 'transcript1'

    >>> # Changes in attributes are propagated to the printable feature

    >>> # for testing, we make sure keys are sorted. Not needed in practice.
    >>> gff.attrs.sort_keys = True
    >>> assert gff.attrs.sort_keys
    >>> print(gff)
    chr1	fake	mRNA	51	300	.	+	.	Awesomeness=99;ID=transcript1;Parent=gene1;


Understanding :class:`Interval` objects is important for using the powerful
filtering and mapping facilities of :class:`BedTool` objects, as described
in the next section.
