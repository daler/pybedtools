.. include:: includeme.rst

.. doctest::
   :hide:

   >>> import pybedtools
   >>> a = pybedtools.example_bedtool('a.bed')
   >>> b = pybedtools.example_bedtool('b.bed')

Chaining methods together (pipe)
--------------------------------

One useful thing about :class:`BedTool` methods is that they often return a
new :class:`BedTool`.  In practice, this means that we can chain together
multiple method calls all in one line, similar to piping on the command
line.

For example, this intersect and merge can be combined into one command:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> # These two lines...
    >>> x1 = a.intersect(b, u=True)
    >>> x2 = x1.merge()

    >>> # ...can be combined into one line:
    >>> x3 = a.intersect(b, u=True).merge()

    >>> x2 == x3
    True

A rule of thumb is that all methods that wrap BEDTools_ programs return
:class:`BedTool` objects, so you can chain these together. Many
:mod:`pybedtools`-unique methods return :class:`BedTool` objects too, just
check the docs (according to :ref:`good docs principle`). For example, as
we saw in one of the examples above, the :meth:`BedTool.saveas` method
returns a :class:`BedTool` object.  That means we can sprinkle those
commands within the example above to save the intermediate steps as
meaningful filenames for later use. For example:

.. doctest::

    >>> x4 = a.intersect(b, u=True).saveas('a-with-b.bed').merge().saveas('a-with-b-merged.bed')

Now we have new files in the current directory called :file:`a-with-b.bed`
and :file:`a-with-b-merged.bed`.  Since :meth:`BedTool.saveas` returns a
:class:`BedTool` object, `x4` points to the :file:`a-with-b-merged.bed`
file.

Sometimes it can be cleaner to separate consecutive calls on each line:

.. doctest::

    >>> x4 = a\
    ... .intersect(b, u=True)\
    ... .saveas('a-with-b.bed')\
    ... .merge()\
    ... .saveas('a-with-b-merged.bed')

Operator overloading
--------------------

There's an even easier way to chain together commands.

I found myself doing intersections so much that I thought it would be
useful to overload the ``+`` and ``-`` operators to do intersections.
To illustrate, these two example commands do the same thing:

.. doctest::
 
    >>> x5 = a.intersect(b, u=True)
    >>> x6 = a + b

    >>> x5 == x6
    True

Just as the `+` operator assumes `intersectBed` with the `-u` arg, the `-`
operator assumes `intersectBed` with the `-v` arg:


.. doctest::

    >>> x7 = a.intersect(b, v=True)
    >>> x8 = a - b

    >>> x7 == x8
    True


If you want to operating on the resulting :class:`BedTool` that is
returned by an addition or subtraction, you'll need to wrap the operation
in parentheses.  This is another way to do the chaining together of the
intersection and merge example from above:

.. doctest:: 

    >>> x9 = (a + b).merge()

And to double-check that all these methods return the same thing:

.. doctest::

    >>> x2 == x3 == x4 == x9
    True


You can learn more about chaining in :ref:`chaining principle`.
