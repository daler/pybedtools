.. include:: includeme.rst

Each
====
Similar to :meth:`BedTool.filter`, which applies a function to return True
or False given an :class:`Interval`, the :meth:`BedTool.each` method applies a
function to return a new, possibly modified :class:`Interval`.

The :meth:`BedTool.each` method applies a function to every feature.  Like
:meth:`BedTool.filter`, you can use your own function or some pre-defined
ones in the :mod:`featurefuncs` module.  Also like :meth:`filter`, `*args`
and `**kwargs` are sent to the function.

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')

    >>> # The results of an "intersect" with c=True will return features
    >>> # with an additional field representing the counts.
    >>> with_counts = a.intersect(b, c=True)

The :func:`featurefuncs.normalized_to_length` function divides the value at
position `N` by the length.  Here we specify `N=-1`, which refers to the
count from the previous step:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> from pybedtools import featurefuncs

    >>> normalized = with_counts.each(featurefuncs.normalized_to_length, -1)

    >>> print normalized
    chr1	1	100	feature1	0	+	0.0
    chr1	100	200	feature2	0	+	1.0000000475e-05
    chr1	150	500	feature3	0	-	2.85714299285e-06
    chr1	900	950	feature4	0	+	2.00000009499e-05
    <BLANKLINE>

