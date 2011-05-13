.. _BedTools as iterators:

BedTools as iterators
=====================
The :meth:`filter` method uses a file-based format, where the new
:class:`BedTool` object refers to a new temp file.  You can use a generator
function to create a new :class:`BedTool` if you want to save disk space:

.. doctest::
    :options: +NORMALIZE_WHITESPACE


    >>> def len_filter(feature, L):
    ...     "Returns True if feature is longer than L"
    ...     return len(feature) > L

    >>> a = pybedtools.example_bedtool('a.bed')

    >>> b = pybedtools.BedTool( (i for i in a if len_filter(i, L=200)))

    >>> print b
    chr1	150	500	feature3	0	-
    <BLANKLINE>

However, keep in mind that printing `b`, which was created using a
generator expression, has now been consumed -- so printing `b` again will
do nothing:

.. doctest::

   >>> print b
   <BLANKLINE>
   <BLANKLINE>

If you create a :class:`BedTool` with a generator expression, you can
always save it as a file for later use. This is what :meth:`filter` is
doing:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.BedTool( (i for i in a if len_filter(i, L=200))).saveas('len-filtered-b.bed')
    >>> print b
    chr1	150	500	feature3	0	-
    <BLANKLINE>

    >>> print b
    chr1	150	500	feature3	0	-
    <BLANKLINE>
