

.. currentmodule:: pybedtools

.. _tempdir: http://docs.python.org/library/tempfile.html#tempfile.tempdir

.. _filo: https://github.com/arq5x/filo

.. _R: http://www.r-project.org/

.. _BEDTools: http://github.com/arq5x/bedtools

.. _BEDTools documentation: http://code.google.com/p/bedtools/#Documentation

.. _Learn Python the Hard Way: http://learnpythonthehardway.org/static/LearnPythonTheHardWay.pdf

.. _IPython: http://ipython.scipy.org/moin/

.. _BED format: http://genome.ucsc.edu/FAQ/FAQformat#format1

.. _tutorial:
   
Tutorial
========


This tutorial assumes that 

1. You know how to use BEDTools_ (if not, check out the 
   `BEDTools documentation`_)
2. You know how to use Python (if not, check out some 
   tutorials like `Learn Python the Hard Way`_)


A brief note on conventions
---------------------------
Throughout this documentation I've tried to use consistent typography, as
follows:

* Python variables and arguments are shown in italics: *s=True*
* Files look like this: :file:`filename.bed`
* Methods, which are often linked to documentation look like this:
  :meth:`BedTool.merge`.
* BEDTools_ programs look like this: ``intersectBed``.
* Arguments that are passed to BEDTools_ programs, as if you were on the
  command line, look like this: ``-d``.
* The ">>>" in the examples below indicates a Python interpreter prompt and
  means to type the code into an interactive Python interpreter like IPython_
  (don't type the >>>) 

Onward!

Create a :class:`BedTool`
-------------------------
First, follow the :ref:`installation` instructions if you haven't already
done so to install both BEDTools_ and :mod:`pybedtools`.

Then import the :mod:`pybedtools` module:

.. doctest::

    >>> import pybedtools

Then set up a :class:`BedTool` instance using a `BED format`_ file. This
can be a file that you already have, or one of the example files as shown
below.  Some methods currently only work for BED files -- see
:ref:`Limitations` for more info on this.  

For this tutorial, we'll use some example files that come with
:mod:`pybedtools`.  We can get the filename for the example files using the
:func:`pybedtools.example_files` function:

.. doctest::
    
    >>> # get an example filename to use
    >>> bed_filename_a = pybedtools.example_filename('a.bed')

The filename will depend on where you have installed :mod:`pybedtools`.
Once you have a filename, creating a :class:`BedTool` object is easy:

.. doctest::
    
    >>> # create a new BedTool using that filename
    >>> a = pybedtools.BedTool(bed_filename_a)

Set up a second one so we can do intersections and subtractions -- this
time, let's make a new :class:`BedTool` all in one line:

.. doctest::

    >>> # create another BedTool to play around with
    >>> b = pybedtools.BedTool(pybedtools.example_filename('b.bed'))

See :ref:`Creating a BedTool` for more information, including convenience
functions for working with example bed files and making :class:`BedTool`
objects directly from strings.

Intersections
-------------

Here's how to intersect *a* with *b*:

.. doctest::

    >>> a_and_b = a.intersect(b)

*a_and_b* is a new :class:`BedTool` instance.  It now points to a temp file
on disk, which is stored in the attribute *a_and_b.fn*; this temp file contains
the intersection of *a* and *b*. 

We can either print the new :class:`BedTool` (which will show ALL features
-- use with caution if you have huge files!) or use the
:meth:`BedTool.head` method to get the first N lines (10 by default).
Here's what *a*, *b*, and *a_and_b* look like:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> a.head()
    chr1    1   100 feature1    0   +
    chr1    100 200 feature2    0   +
    chr1    150 500 feature3    0   -
    chr1    900 950 feature4    0   +
   
    >>> b.head()
    chr1    155 200 feature5    0   -
    chr1    800 901 feature6    0   +

    >>> a_and_b.head()
    chr1    155 200 feature2    0   +
    chr1    155 200 feature3    0   -
    chr1    900 901 feature4    0   +

The :meth:`BedTool.intersect` method simply wraps the BEDTools_ program
``intersectBed``.  This means that we can pass :meth:`BedTool.intersect`
any arguments that ``intersectBed`` accepts.  For example, if we want to
use the ``intersectBed`` switch ``-u`` (which acts as a True/False switch
to indicate that we want to see the features in *a* that overlapped
something in *b*), then we can use the keyword argument *u=True*, like this:


.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> # Intersection using the -u switch
    >>> a_with_b = a.intersect(b, u=True)
    >>> a_with_b.head()
    chr1    100 200 feature2    0   +
    chr1    150 500 feature3    0   -
    chr1    900 950 feature4    0   +
    
*a_with_b* is another, different temp file whose name is stored in
*a_with_b.fn*.  You can read more about the use of temp files in
:ref:`temp principle`.  More on arguments that you can pass to
:class:`BedTool` objects in a moment, but first, some info about saving
files.

Saving the results
------------------

If you want to save the results as a meaningful filename for later use, use
the :meth:`BedTool.saveas` method.  This also lets you optionally specify a
trackline for directly uploading to the UCSC Genome Browser, instead of
opening up the files afterward and manually adding a trackline:

.. doctest::

    >>> a_with_b.saveas('intersection-of-a-and-b.bed', trackline='track name="a and b"')
    <BedTool(intersection-of-a-and-b.bed)>

Note that the :meth:`BedTool.saveas` method returns a new :class:`BedTool`
object which points to the newly created file on disk.  This allows you to
insert a :meth:`BedTool.saveas` call in the middle of a chain of commands
(described in another section below).


Default arguments
-----------------
Recall that we passed the *u=True* argument to :meth:`a.intersect`::

    >>> a_with_b = a.intersect(b, u=True)

While we're on the subject of arguments, note that we didn't have to
specify *-a* or *-b* arguments, like you would need if calling
``intersectBed`` from the command line.  That's because :class:`BedTool`
objects make some assumptions for convenience.  

We could have supplied the arguments *a=a.fn* and *b=b.fn*.  But since
we're calling a method on *a*, :mod:`pybedtools` assumes that the file *a*
points to (specifically, *a.fn*) is the one we want to use as input.  So by
default, we don't need to explicitly give the keyword argument *a=a.fn*
because the :meth:`a.intersect` method does so automatically.

We're also calling a method that takes a second bed file as input  -- other
such methods include :meth:`BedTool.subtract` and :meth:`BedTool.closest`.
In these cases, :mod:`pybedtools` assumes the first unnamed argument to
these methods are the second file you want to operate on (and if you pass a
:class:`BedTool`, it'll automatically use the file in the *fn* attribute of
that :class:`BedTool`).  So ``a.intersect(b)`` is just a more convenient
form of ``a.intersect(a=a.fn, b=b.fn)``, which does the same thing.

OK, enough about arguments for now, but you can read more about them in
:ref:`similarity principle`, :ref:`default args principle` and :ref:`non
defaults principle`. 

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

Methods specific to :mod:`pybedtools`
-------------------------------------

In no particular order, here are some other useful things that
:class:`BedTool` objects can do.

Counting 
~~~~~~~~
We can easily count features in :class:`BedTool` objects:

.. doctest::

    >>> a.count()
    4

    >>> a.intersect(b, u=True).merge().count()
    2

Lengths of features
~~~~~~~~~~~~~~~~~~~
We can get the lengths of features, which is useful for things like getting
a histogram of binding site sizes after running a peak-calling algorithm on
ChIP-seq data:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> a.head()
    chr1    1   100 feature1    0   +
    chr1    100 200 feature2    0   +
    chr1    150 500 feature3    0   -
    chr1    900 950 feature4    0   +
    
    >>> [len(i) for i in a]
    [99, 100, 350, 50]


Creating genome files ('chromSizes')
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some BEDTools_ programs require a "chromSizes" file to prevent out-of-range
coordinates.  :class:`BedTool` objects have a convenience method to get
this file for you and return the resulting filename:

.. doctest:: 
    :options: +NORMALIZE_WHITESPACE
    
    >>> # Download Drosophila melanogaster chromSizes table from UCSC
    >>> chromsizes = pybedtools.get_chromsizes_from_ucsc('hg19', saveas='hg19.genome')

    >>> # print first few lines of that file
    >>> f = open('hg19.genome')
    >>> for i in range(5):
    ...     print f.readline(),
    chr1	249250621
    chr10	135534747
    chr11	135006516
    chr11_gl000202_random	40103
    chr12	133851895


Feature centers
~~~~~~~~~~~~~~~
We can get the midpoints of features, with :meth:`BedTool.feature_centers`,
which is useful for things like extracting out the center N bp of features
for de novo motif detection.

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> from pybedtools.featurefuncs import center
    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = a.each(center, width=10)
    >>> print b
    chr1	45	55	feature1	0	+
    chr1	145	155	feature2	0	+
    chr1	320	330	feature3	0	-
    chr1	920	930	feature4	0	+

Filtering
~~~~~~~~~
Only get features of a certain size:

.. todo::
   
   TODO: reimplement using the ``each()`` method.



Working with counted intersections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several methods that help you deal with counted intersections --
the files that are returned when you use the *c=True* kwarg with
:meth:`BedTool.intersect`.  First, do an intersection:

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')
    >>> print a.intersect(b, c=True)
    chr1    1   100 feature1    0   +   0
    chr1    100 200 feature2    0   +   1
    chr1    150 500 feature3    0   -   1
    chr1    900 950 feature4    0   +   1

You can retrieve these counts later using the :meth:`BedTool.counts` method:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> result = a.intersect(b, c=True)
    >>> print list(result.counts())
    [0, 1, 1, 1]


    

.. _`working with history`:

Using the history and tags
--------------------------
`BEDTools`_ makes it very easy to do rather complex genomic algebra.  Sometimes
when you're doing some exploratory work, you'd like to rewind back to a
previous step, or clean up temporary files that have been left on disk over the
course of some experimentation.

To assist this sort of workflow, :class:`BedTool` instances keep track of
their history in the :attr:`BedTool.history` attribute.  Let's make an
example :class:`BedTool`, *c*, that has some history:

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> a = pybedtools.example_bedtool('a.bed')
    >>> b = pybedtools.example_bedtool('b.bed')
    >>> c = a.intersect(b, u=True)


*c* now has a history which tells you all sorts of useful things (described
in more detail below)::

    >>> print c.history
    [<HistoryStep> bedtool("/home/ryan/pybedtools/pybedtools/test/a.bed").intersect("/home/ryan/pybedtools/pybedtools/test/b.bed", u=True), parent tag: klkreuay, result tag: egzgnrvj]


There are several things to note here.  First, the history describes the full
commands, including all the names of the temp files and all the arguments that
you would need to run in order to re-create it.  Since :class:`BedTool` objects
are fundamentally file-based, the command refers to the underlying filenames
(i.e., :file:`a.bed` and :file:`b.bed`) instead of the :class:`BedTool`
instances (i.e., *a* and *b*). A simple copy-paste of the command will be
enough re-run the command. While this may be useful in some situations, be
aware that if you do run the command again you'll get *another* temp file that
has the same contents as *c*'s temp file.

To avoid such cluttering of your temp dir, the history also reports
**tags**. :class:`BedTool` objects, when created, get a random tag assigned
to them.  You can get get the :class:`BedTool` associated with tag with the
:func:`pybedtools.find_tagged` function. These tags are used to keep track
of instances during this session.

So in this case, we could get a reference to the *a* instance with::

    >>> should_be_a = pybedtools.find_tagged('klkreuay')

Here's confirmation that the parent of the first step of *c*'s history is
*a* (note that :class:`HistoryStep` objects have a
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

Let's see what the history of *f* (the last :class:`BedTool` created) looks
like . . . note that here I'm formatting the results to make it easier to
see::

    >>> print f.history
    [
    |   [
    |   |   [
    |   |   |   [
    |   |   |   |<HistoryStep> BedTool("/home/ryan/pybedtools/pybedtools/test/a.bed").intersect(
    |   |   |   |                      "/home/ryan/pybedtools/pybedtools/test/b.bed", 
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
    |                      "/home/ryan/pybedtools/pybedtools/test/b.bed",
    |                      ),
    |                      parent tag: zlwqblvk, 
    |                      result tag: reztxhen
    ]

Those first three history steps correspond to *c*, *d*, and *e*
respectively, as we can see by comparing the code snippet above with the
commands in each history step.  In other words, *e* can be described by the
sequence of 3 commands in the first three history steps.  In fact, if we
checked *e.history*, we'd see exactly those same 3 steps.

When *f* was created above, it operated both on *e*, which had its own
history, as well as *b* -- note the nesting of the list. You can do
arbitrarily complex "genome algebra" operations, and the history of the
:class:`BEDTools` will keep track of this.  It may not be useful in every
situtation, but the ability to backtrack and have a record of what you've
done can sometimes be helpful.

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

Note that the file that *f* points to is left alone.  To clarify, the
:meth:`BedTool.delete_temporary_history` will only delete temp files that match
the pattern ``<TEMP_DIR>/pybedtools.*.tmp`` from the history of *f*, up to but
not including the file for *f* itself.  Any :class:`BedTool` instances that do
not match the pattern are left alone.

More documentation
------------------
For more info, see the :ref:`topical`.

.. doctest::
    :hide:

    Gotta clean up all the files created over the course of the tutorial...

    >>> fns_to_remove = ['a-with-b.bed', 'a-with-b-merged.bed', 'hg19.genome', 'intersection-of-a-and-b.bed','middle-100-bp.bed','shared_merged.bed']
    >>> for fn in fns_to_remove:
    ...     if os.path.exists(fn):
    ...         os.unlink(fn)
