.. currentmodule:: pybedtools


Default arguments
=================
Recall in the earlier :ref:`intersections` section that we passed the `u=True` argument to :meth:`a.intersect`:

.. doctest::

    >>> import pybedtools
    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')
    >>> a_with_b = a.intersect(b, u=True)

Let's do the same thing but use different variable names for the :class:`BedTool` objects so that
the next section is less confusing:

.. doctest::

    >>> import pybedtools
    >>> exons = pybedtools.example_bedtool('a.bed')
    >>> snps = pybedtools.example_bedtool('b.bed')
    >>> exons_with_snps = exons.intersect(snps, u=True)


While we're on the subject of arguments, note that we didn't have to specify
`-a` or `-b` arguments, like you would need if calling `intersectBed` from the
command line.  In other words, since `exons` refers to the file `a.bed` and
`snps` refers to the file `b.bed`, the following line::

    >>> exons_with_snps = exons.intersect(snps, u=True)

is equivalent to the command line usage of::

    $ intersectBed -a a.bed -b b.bed -u > tmpfile

But we didn't have to explicitly pass the argument for `-a` because
:class:`BedTool` objects make some assumptions for convenience.

We're calling a method on the :class:`BedTool` object `exons`, so
:mod:`pybedtools` assumes that the file `exons` points to (stored in the
attribute `exons.fn`) is the one we want to use as input.  So by default, we
don't need to explicitly give the keyword argument `a=exons.fn` because the
:meth:`exons.intersect` method does so automatically.

We're also calling a method that takes a second bed file as input  -- other
such methods include :meth:`BedTool.subtract` and :meth:`BedTool.closest`,
and others.  For these methods, in addition to assuming `-a` is taken care
of by the :attr:`BedTool.fn` attribute, :mod:`pybedtools` also assumes the
first unnamed argument to these methods are the second file you want to
operate on (and if you pass a :class:`BedTool`, it'll automatically use the
file in the `fn` attribute of that :class:`BedTool`).

An example may help to illustrate: these different ways of calling
:meth:`BedTool.intersect` all have the same results, with the first version
being the most compact (and probably most convenient):

.. doctest::

    >>> # these all have identical results
    >>> x1 = exons.intersect(snps)
    >>> x2 = exons.intersect(a=exons.fn, b=snps.fn)
    >>> x3 = exons.intersect(b=snps.fn)
    >>> x4 = exons.intersect(snps, a=exons.fn)
    >>> x1 == x2 == x3 == x4
    True

Note that `a.intersect(a=a.fn, b)` is not a valid Python expression, since
non-keyword arguments must come before keyword arguments, but
`a.intersect(b, a=a.fn)` works fine.

If you're ever unsure, the docstring for these methods indicates which, if
any, arguments are used as default.  For example, in the
:meth:`BedTool.intersect` help, it says::

    For convenience, the file or stream this BedTool points to is implicitly
    passed as the -a argument to intersectBed

OK, enough about arguments for now, but you can read more about them in
:ref:`similarity principle`, :ref:`default args principle` and :ref:`non
defaults principle`.
