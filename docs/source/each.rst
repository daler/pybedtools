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

Let's define a function that will take the number of counts in each feature
as calculated above and divide by the number of bases in that feature.  We
can also supply an optional scalar, like 0.001, to get the results in
"number of intersections per kb".  We then insert that value into the score
field of the feature.  Here's the function:

.. doctest::

    >>> def normalize_count(feature, scalar=0.001):
    ...     """
    ...     assume feature's last field is the count
    ...     """
    ...     counts = float(feature[-1])
    ...     normalized = counts / len(feature) * scalar
    ...
    ...     # need to convert back to string to insert into feature
    ...     feature.score = str(normalized)
    ...     return feature

And we apply it like this:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> normalized = with_counts.each(normalize_count)
    >>> print normalized
    chr1	1	100	feature1	0.0	+	0
    chr1	100	200	feature2	1e-05	+	1
    chr1	150	500	feature3	2.85714285714e-06	-	1
    chr1	900	950	feature4	2e-05	+	1
    <BLANKLINE>
