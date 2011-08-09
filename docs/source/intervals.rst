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

We can access the :class:`Intervals` of `a` several different ways. The
most common use is probably by using the :class:`BedTool` `a` as an
iterator.  For now though, let's look at a single :class:`Interval`:

.. doctest::

    >>> feature = iter(a).next()


Common :class:`Interval` attributes
-----------------------------------
Printing a feature converts it into the original line from the file:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print feature
    chr1	1	100	feature1	0	+

All features have `chrom`, `start`, `stop`, `name`, `score`, and `strand`
attributes.  Note that `start` and `stop` are long integers, while
everything else (including `score`) is a string.

.. doctest::

    >>> feature.chrom
    'chr1'

    >>> feature.start
    1L

    >>> feature.stop
    100L

    >>> feature.name
    'feature1'

    >>> feature.score
    '0'

    >>> feature.strand
    '+'

Let's make another feature that only has chrom, start, and stop to see how
:mod:`pybedtools` deals with missing attributes:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> feature2 = iter(pybedtools.BedTool('chrX 500 1000', from_string=True)).next()

    >>> print feature2
    chrX	500	1000


    >>> feature2.chrom
    'chrX'

    >>> feature2.start
    500L

    >>> feature2.stop
    1000L

    >>> feature2.name
    ''

    >>> feature2.score
    ''

    >>> feature2.strand
    ''

This illustrates that default values are empty strings.


Indexing into :class:`Interval` objects
---------------------------------------

:class:`Interval` objects can also be indexed by position into the original
line (like a list) or indexed by name of attribute (like a dictionary).

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print feature
    chr1	1	100	feature1	0	+


    >>> feature[0]
    'chr1'

    >>> feature['chrom']
    'chr1'

    >>> feature[1]
    '1'

    >>> feature['start']
    1L


Fields
------
:class:`Interval` objects have a :attr:`Interval.fields` attribute that
contains the original line split into a list of strings.  When an integer
index is used on the :class:`Interval` (for example, `feature[3]`), it is
the `fields` attribute that is actually being indexed into.

.. doctest::

    >>> f = iter(pybedtools.BedTool('chr1 1 100 asdf 0 + a b c d', from_string=True)).next()
    >>> f.fields
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

    >>> gff.file_type
    'gff'
    >>> bed.file_type
    'bed'

Printing the :class:`Intervals` shows that the strings are in the appropriate
coordinates:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print gff
    chr1	fake	mRNA	51	300	.	+	.	ID=mRNA1;Parent=gene1;

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print bed
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

    >>> bed.start
    50L
    >>> bed[1]
    '50'

. . . but for the GFF feature, they differ -- the `start` attribute is
zero-based while the string representation (the fourth field of a GFF file)
remains in one-based GFF coords:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> gff.start
    50L
    >>> gff[3]
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
    >>> print gff
    chr1	fake	mRNA	51	300	.	+	.	ID=mRNA1;Parent=gene1;

    >>> # original attributes
    >>> gff.attrs
    {'ID': 'mRNA1', 'Parent': 'gene1'}

    >>> # add some new attributes
    >>> gff.attrs['Awesomeness'] = 99
    >>> gff.attrs['ID'] = 'transcript1'

    >>> # Changes in attributes are propagated to the printable feature
    >>> print gff
    chr1	fake	mRNA	51	300	.	+	.	Awesomeness=99;ID=transcript1;Parent=gene1


Understanding :class:`Interval` objects is important for using the powerful
filtering and mapping facilities of :class:`BedTool` objects, as described
in the next section.
