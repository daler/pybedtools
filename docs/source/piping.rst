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
    >>> a_with_b = a.intersect(b, u=True)
    >>> a_with_b_merged_1 = a_with_b.merge()

    >>> # Could be combined into one line:
    >>> a_with_b_merged_2 = a.intersect(b, u=True).merge()

    >>> a_with_b_merged_1 == a_with_b_merged_2
    True


In general, methods that return :class:`BedTool` objects have the following text in
their docstring to indicate this::

        .. note::

            This method returns a new BedTool instance

A rule of thumb is that all methods that wrap BEDTools_ programs return
:class:`BedTool` objects, so you can chain these together. Other
:mod:`pybedtools`-unique methods return :class:`BedTool` objects too, just
check the docs (according to :ref:`good docs principle`). For example, as
we saw in one of the examples above, the :meth:`BedTool.saveas` method
returns a :class:`BedTool` object.  That means we can sprinkle those
commands within the example above to save the intermediate steps as
meaningful filenames for later use. For example:

.. doctest::

    >>> a_with_b_merged_3 = a.intersect(b, u=True).saveas('a-with-b.bed').merge().saveas('a-with-b-merged.bed')

Now we have new files in the current directory called :file:`a-with-b.bed`
and :file:`a-with-b-merged.bed`.  Since :meth:`BedTool.saveas` returns a
:class:`BedTool` object, ``x2`` points to the :file:`a-with-b-merged.bed`
file.

There's an even easier way to chain together commands.

I found myself doing intersections so much that I thought it would be
useful to overload the ``+`` and ``-`` operators to do intersections.
To illustrate, these two example commands do the same thing:

.. doctest::
 
    >>> a_with_b_1 = a.intersect(b, u=True)
    >>> a_with_b_2 = a+b

    >>> a_with_b_1 == a_with_b_2
    True

.. doctest::

    >>> a_without_b_1 = a.intersect(b, v=True)
    >>> a_without_b_2 = a-b

    >>> a_without_b_1 == a_without_b_2
    True

Note that the ``+`` operator assumes the ``-u`` option and the ``-``
operator assumes ``intersectBed``'s ``-v`` option:

If you want to operating on the resulting :class:`BedTool` that is
returned by an addition or subtraction, you'll need to wrap the operation
in parentheses:

.. doctest:: 

    >>> a_with_b_merged_4 = (a+b).merge()

And to double-check that all these methods return the same thing:

.. doctest::

    >>> a_with_b_merged_1 == a_with_b_merged_2 == a_with_b_merged_3 == a_with_b_merged_4
    True



You can learn more about chaining in :ref:`chaining principle`.
