.. include:: includeme.rst

Chaining methods together (pipe)
--------------------------------

One useful thing about :class:`BedTool` methods is that they often return a
new :class:`BedTool`.  In practice, this means that we can chain together
multiple method calls all in one line, similar to piping on the command
line.

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> # Intersect and then merge all on one line, displaying the first
    >>> # 10 lines of the results
    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')
    >>> a.intersect(b, u=True).merge().head()
    chr1    100 500
    chr1    900 950


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
meaningful filenames for later use:

.. doctest::

    >>> a.intersect(b, u=True).saveas('a-with-b.bed').merge().saveas('a-with-b-merged.bed')
    <BedTool(a-with-b-merged.bed)>

Now we have new files in the current directory called :file:`a-with-b.bed`
and :file:`a-with-b-merged.bed`.  


I found myself doing intersections so much that I thought it would be
useful to overload the ``+`` and ``-`` operators to do intersections.
To illustrate, these two example commands do the same thing:

.. doctest::
 
    >>> result_1 = a.intersect(b, u=True)
    >>> result_2 = a+b

    >>> # To test equality, convert to strings
    >>> str(result_1) == str(result_2)
    True

And the ``-`` operator assumes ``intersectBed``'s ``-v`` option:

.. doctest::

    >>> result_1 = a.intersect(b, v=True)
    >>> result_2 = a-b

    >>> # To test equality, convert to strings
    >>> str(result_1) == str(result_2)
    True

If you want to operating on the resulting :class:`BedTool` that is
returned by an addition or subtraction, you'll need to wrap the operation
in parentheses:

.. doctest:: 

    >>> merged = (a+b).merge()

You can learn more about chaining in :ref:`chaining principle`.
