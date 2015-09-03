.. include:: includeme.rst

.. _`working with history`:

Using the history and tags
--------------------------
`BEDTools`_ makes it very easy to do rather complex genomic algebra.  Sometimes
when you're doing some exploratory work, you'd like to rewind back to a
previous step, or clean up temporary files that have been left on disk over the
course of some experimentation.

To assist this sort of workflow, :class:`BedTool` instances keep track of
their history in the :attr:`BedTool.history` attribute.  Let's make an
example :class:`BedTool`, `c`, that has some history:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')
    >>> c = a.intersect(b, u=True)


`c` now has a history which tells you all sorts of useful things (described
in more detail below)::

    >>> print c.history
    [<HistoryStep> bedtool("/home/ryan/pybedtools/pybedtools/test/a.bed").intersect("/home/ryan/pybedtools/pybedtools/test/b.bed", u=True), parent tag: klkreuay, result tag: egzgnrvj]


There are several things to note here.  First, the history describes the full
commands, including all the names of the temp files and all the arguments that
you would need to run in order to re-create it.  Since :class:`BedTool` objects
are fundamentally file-based, the command refers to the underlying filenames
(i.e., :file:`a.bed` and :file:`b.bed`) instead of the :class:`BedTool`
instances (i.e., `a` and `b`). A simple copy-paste of the command will be
enough re-run the command. While this may be useful in some situations, be
aware that if you do run the command again you'll get *another* temp file that
has the same contents as `c`'s temp file.

To avoid such cluttering of your temp dir, the history also reports
**tags**. :class:`BedTool` objects, when created, get a random tag assigned
to them.  You can get get the :class:`BedTool` associated with tag with the
:func:`pybedtools.find_tagged` function. These tags are used to keep track
of instances during this session.

So in this case, we could get a reference to the `a` instance with::

    >>> should_be_a = pybedtools.find_tagged('klkreuay')

Here's confirmation that the parent of the first step of `c`'s history is
`a` (note that :class:`HistoryStep` objects have a
:attr:`HistoryStep.parent_tag` and :attr:`HistoryStep.result_tag`):

.. doctest::

    >>> pybedtools.find_tagged(c.history[0].parent_tag) == a
    True

Let's make something with a more complicated history:

.. doctest::

    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')
    >>> c = a.intersect(b)
    >>> d = c.slop(g=pybedtools.chromsizes('hg19'), b=1)
    >>> e = d.merge()

    >>> # this step adds complexity!
    >>> f = e.subtract(b)

Let's see what the history of `f` (the last :class:`BedTool` created) looks
like . . . note that here I'm formatting the results to make it easier to
see::

    >>> print f.history
    [
    |   [
    |   |   [
    |   |   |   [
    |   |   |   |<HistoryStep> BedTool("/usr/local/lib/python2.6/dist-packages/pybedtools/test/data/a.bed").intersect(
    |   |   |   |                      "/usr/local/lib/python2.6/dist-packages/pybedtools/test/data/b.bed", 
    |   |   |   |                      ), 
    |   |   |   |                      parent tag: rzrztxlw, 
    |   |   |   |                      result tag: ifbsanqk
    |   |   |   ],
    |   |   |
    |   |   |<HistoryStep> BedTool("/tmp/pybedtools.BgULVj.tmp").slop(
    |   |   |                      b=1,genome="hg19"
    |   |   |                      ), 
    |   |   |                      parent tag: ifbsanqk, 
    |   |   |                      result tag: omfrkwjp
    |   |   ],
    |   |
    |   |<HistoryStep> BedTool("/tmp/pybedtools.SFmbYc.tmp").merge(),
    |   |                      parent tag: omfrkwjp,
    |   |                      result tag: zlwqblvk
    |   ], 
    |
    |<HistoryStep> BedTool("/tmp/pybedtools.wlBiMo.tmp").subtract(
    |                      "/usr/local/lib/python2.6/dist-packages/pybedtools/test/data/b.bed",
    |                      ),
    |                      parent tag: zlwqblvk, 
    |                      result tag: reztxhen
    ]

Those first three history steps correspond to `c`, `d`, and `e`
respectively, as we can see by comparing the code snippet above with the
commands in each history step.  In other words, `e` can be described by the
sequence of 3 commands in the first three history steps.  In fact, if we
checked `e.history`, we'd see exactly those same 3 steps.

When `f` was created above, it operated both on `e`, which had its own
history, as well as `b` -- note the nesting of the list. You can do
arbitrarily complex "genome algebra" operations, and the history of the
:class:`BEDTools` will keep track of this.  It may not be useful in every
situtation, but the ability to backtrack and have a record of what you've
done can sometimes be helpful.

Deleting temp files specific to a single :class:`BedTool`
---------------------------------------------------------
You can delete temp files that have been created over the history of a
:class:`BedTool` with :meth:`BedTool.delete_temporary_history`.  This method
will inspect the history, figure out which items point to files in the temp dir
(which you can see with :func:`get_tempdir`), and prompt you for their
deletion::

    >>> f.delete_temporary_history()
    Delete these files?
        /tmp/pybedtools..BgULVj.tmp
        /tmp/pybedtools.SFmbYc.tmp
        /tmp/pybedtools.wlBiMo.tmp
    (y/N) y

Note that the file that `f` points to is left alone.  To clarify, the
:meth:`BedTool.delete_temporary_history` will only delete temp files that match
the pattern ``<TEMP_DIR>/pybedtools.*.tmp`` from the history of `f`, up to but
not including the file for `f` itself.  Any :class:`BedTool` instances that do
not match the pattern are left alone.  Use the kwarg `ask=False` to disable
the prompt.
