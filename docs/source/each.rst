.. include:: includeme.rst

Each
~~~~
The :meth:`BedTool.each` method applies a function to every feature.  Like
:meth:`BedTool.filter`, you can use your own function or some pre-defined ones
in the :mod:`featurefuncs` module.

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> from pybedtools import featurefuncs
    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')
    >>> with_counts = a.intersect(b, c=True)
    >>> normalized = with_counts.each(featurefuncs.normalized_to_length, -1)
    >>> print normalized
    chr1	1	100	feature1	0	+	0.0
    chr1	100	200	feature2	0	+	1.0000000475e-05
    chr1	150	500	feature3	0	-	2.85714299285e-06
    chr1	900	950	feature4	0	+	2.00000009499e-05
    <BLANKLINE>

