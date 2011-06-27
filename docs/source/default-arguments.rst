.. currentmodule:: pybedtools

.. doctest::
   :hide:

   >>> import pybedtools
   >>> a = pybedtools.example_bedtool('a.bed')
   >>> b = pybedtools.example_bedtool('b.bed')

Default arguments
=================
Recall that we passed the `u=True` argument to :meth:`a.intersect`:

.. doctest::

    >>> a_with_b = a.intersect(b, u=True)

While we're on the subject of arguments, note that we didn't have to
specify `-a` or `-b` arguments, like you would need if calling
`intersectBed` from the command line.  That's because :class:`BedTool`
objects make some assumptions for convenience.

We could have supplied the arguments ``a=a.fn`` and ``b=b.fn``:


.. doctest::

    >>> another_way = a.intersect(a=a.fn, b=b.fn, u=True)
    >>> another_way == a_with_b
    True

But since we're calling a method on the :class:`BedTool` object `a`,
:mod:`pybedtools` assumes that the file `a` points to (stored in the
attribute `a.fn`) is the one we want to use as input.  So by default, we
don't need to explicitly give the keyword argument `a=a.fn` because the
:meth:`a.intersect` method does so automatically.

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
    >>> x1 = a.intersect(b)
    >>> x2 = a.intersect(a=a.fn, b=b.fn)
    >>> x3 = a.intersect(b=b.fn)
    >>> x4 = a.intersect(b, a=a.fn)
    >>> x1 == x2 == x3 == x4
    True

Note that `a.intersect(a=a.fn, b)` is not a valid Python expression, since
non-keyword arguments must come before keyword arguments, but
`a.intersect(b, a=a.fn)` works fine.

If you're ever unsure, the docstring for these methods indicates which, if
any, arguments are used as default.  For example, in the
:meth:`BedTool.intersect` help, it says::

    .. note::

    For convenience, the file this bedtool object points to is
    passed as "-a"

OK, enough about arguments for now, but you can read more about them in
:ref:`similarity principle`, :ref:`default args principle` and :ref:`non
defaults principle`.
