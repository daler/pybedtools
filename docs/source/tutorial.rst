

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
  :meth:`bedtool.merge`.
* BEDTools_ programs look like this: ``intersectBed``.
* Arguments that are passed to BEDTools_ programs, as if you were on the
  command line, look like this: ``-d``.
* The ">>>" in the examples below indicates a Python interpreter prompt and
  means to type the code into an interactive Python interpreter like IPython_
  (don't type the >>>) 

Onward!

Create a :class:`bedtool`
-------------------------
First, follow the :ref:`installation` instructions if you haven't already
done so to install both BEDTools_ and :mod:`pybedtools`.

Then import the :mod:`pybedtools` module:

.. doctest::

    >>> import pybedtools

Then set up a :class:`bedtool` instance using a `BED format`_ file. This
can be a file that you already have, or one of the example files as shown
below.  Currently, only BED format files are supported -- see
:ref:`Limitations` for more info on this.  

For this tutorial, we'll use some example files that come with
:mod:`pybedtools`.  We can get the filename for the example files using the
:func:`pybedtools.example_files` function:

.. doctest::
    
    >>> # get an example filename to use
    >>> bed_filename_a = pybedtools.example_bed_fn('a.bed')

The filename will depend on where you have installed :mod:`pybedtools`.
Once you have a filename, creating a :class:`bedtool` object is easy:

.. doctest::
    
    >>> # create a new bedtool using that filename
    >>> a = pybedtools.bedtool(bed_filename_a)

Set up a second one so we can do intersections and subtractions -- this
time, let's make a new :class:`bedtool` all in one line:

.. doctest::

    >>> # create another bedtool to play around with
    >>> b = pybedtools.bedtool(pybedtools.example_bed_fn('b.bed'))

See :ref:`Creating a bedtool` for more information, including convenience
functions for working with example bed files and making :class:`bedtool`
objects directly from strings.

Intersections
-------------

Here's how to intersect *a* with *b*:

.. doctest::

    >>> a_and_b = a.intersect(b)

*a_and_b* is a new :class:`bedtool` instance.  It now points to a temp file
on disk, which is stored in the attribute *a_and_b.fn*; this temp file contains
the intersection of *a* and *b*. 

We can either print the new :class:`bedtool` (which will show ALL features
-- use with caution if you have huge files!) or use the
:meth:`bedtool.head` method to get the first N lines (10 by default).
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

The :meth:`bedtool.intersect` method simply wraps the BEDTools_ program
``intersectBed``.  This means that we can pass :meth:`bedtool.intersect`
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
:class:`bedtool` objects in a moment, but first, some info about saving
files.

Saving the results
------------------

If you want to save the results as a meaningful filename for later use, use
the :meth:`bedtool.saveas` method.  This also lets you optionally specify a
trackline for directly uploading to the UCSC Genome Browser, instead of
opening up the files afterward and manually adding a trackline:

.. doctest::

    >>> a_with_b.saveas('intersection-of-a-and-b.bed', trackline='track name="a and b"')
    <bedtool (intersection-of-a-and-b.bed)>

Note that the :meth:`bedtool.saveas` method returns a new :class:`bedtool`
object which points to the newly created file on disk.  This allows you to
insert a :meth:`bedtool.saveas` call in the middle of a chain of commands
(described in another section below).


Default arguments
-----------------
Recall that we passed the *u=True* argument to :meth:`a.intersect`::

    >>> a_with_b = a.intersect(b, u=True)

While we're on the subject of arguments, note that we didn't have to
specify *-a* or *-b* arguments, like you would need if calling
``intersectBed`` from the command line.  That's because :class:`bedtool`
objects make some assumptions for convenience.  

We could have supplied the arguments *a=a.fn* and *b=b.fn*.  But since
we're calling a method on *a*, :mod:`pybedtools` assumes that the file *a*
points to (specifically, *a.fn*) is the one we want to use as input.  So by
default, we don't need to explicitly give the keyword argument *a=a.fn*
because the :meth:`a.intersect` method does so automatically.

We're also calling a method that takes a second bed file as input  -- other
such methods include :meth:`bedtool.subtract` and :meth:`bedtool.closest`.
In these cases, :mod:`pybedtools` assumes the first unnamed argument to
these methods are the second file you want to operate on (and if you pass a
:class:`bedtool`, it'll automatically use the file in the *fn* attribute of
that :class:`bedtool`).  So ``a.intersect(b)`` is just a more convenient
form of ``a.intersect(a=a.fn, b=b.fn)``, which does the same thing.

OK, enough about arguments for now, but you can read more about them in
:ref:`similarity principle`, :ref:`default args principle` and :ref:`non
defaults principle`. 

Chaining methods together (pipe)
--------------------------------

One useful thing about :class:`bedtool` methods is that they often return a
new :class:`bedtool`.  In practice, this means that we can chain together
multiple method calls all in one line, similar to piping on the command
line.

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> # Intersect and then merge all on one line, displaying the first
    >>> # 10 lines of the results
    >>> a.intersect(b, u=True).merge().head()
    chr1    100 500
    chr1    900 950


In general, methods that return :class:`bedtool` objects have the following text in
their docstring to indicate this::

        .. note::

            This method returns a new bedtool instance

A rule of thumb is that all methods that wrap BEDTools_ programs return
:class:`bedtool` objects, so you can chain these together. Other
:mod:`pybedtools`-unique methods return :class:`bedtool` objects too, just
check the docs (according to :ref:`good docs principle`). For example, as
we saw in one of the examples above, the :meth:`bedtool.saveas` method
returns a :class:`bedtool` object.  That means we can sprinkle those
commands within the example above to save the intermediate steps as
meaningful filenames for later use:

.. doctest::

    >>> a.intersect(b, u=True).saveas('a-with-b.bed').merge().saveas('a-with-b-merged.bed')
    <bedtool (a-with-b-merged.bed)>

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

If you want to operating on the resulting :class:`bedtool` that is
returned by an addition or subtraction, you'll need to wrap the operation
in parentheses:

.. doctest:: 

    >>> merged = (a+b).merge()

You can learn more about chaining in :ref:`chaining principle`.

Methods specific to :mod:`pybedtools`
-------------------------------------

In no particular order, here are some other useful things that
:class:`bedtool` objects can do.

Counting 
~~~~~~~~
We can easily count features in :class:`bedtool` objects:

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
    
    >>> a.lengths()
    [99, 100, 350, 50]


Creating genome files ('chromSizes')
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some BEDTools_ programs require a "chromSizes" file to prevent out-of-range
coordinates.  :class:`bedtool` objects have a convenience method to get
this file for you and return the resulting filename:

.. doctest:: 
    :options: +NORMALIZE_WHITESPACE
    
    >>> # Download Drosophila melanogaster chromSizes table from UCSC
    >>> chromsizes = a.get_genome('hg19', fn='hg19.genome')

    >>> # print first few lines of that file
    >>> f = open(chromsizes)
    >>> for i in range(5):
    ...     print f.readline(),
    chr1    249250621
    chr2    243199373
    chr3    198022430
    chr4    191154276
    chr5    180915260


Feature centers
~~~~~~~~~~~~~~~
We can get the midpoints of features, with :meth:`bedtool.feature_centers`,
which is useful for things like extracting out the center N bp of features
for de novo motif detection.

.. doctest::

    >>> a.feature_centers(n=50).saveas('middle-100-bp.bed')
    <bedtool (middle-100-bp.bed)>

Filtering
~~~~~~~~~
Only get features of a certain size:

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> # only features that are smaller than 60 bp
    >>> a.size_filter(min=0, max=60).head()
    chr1    900 950 feature4    0.0 +


Working with counted intersections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several methods that help you deal with counted intersections --
the files that are returned when you use the *c=True* kwarg with
:meth:`bedtool.intersect`.  First, do an intersection:

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> print a.intersect(b, c=True)
    chr1    1   100 feature1    0   +   0
    chr1    100 200 feature2    0   +   1
    chr1    150 500 feature3    0   -   1
    chr1    900 950 feature4    0   +   1

You can retrieve these counts later using the :meth:`bedtool.counts` method:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> result = a.intersect(b, c=True)
    >>> print result.counts()
    [0, 1, 1, 1]



Topical documentation
=====================

Future work
-----------
The following is a list of items I plan to work on for future releases of
:mod:`pybedtools`.  Help and feedback are welcome!:

* Better mechanism for handling temp files.  

    * :class:`bedtool` objects could keep track of all their 'parent'
      tempfiles, and as such would retain the history of their creation.

    * indexing into the history of a :class:`bedtool` would give you access
      to these previous files
    
    * could have a :meth:`bedtool.delete_history`, which would delete all the
      intermediate tempfiles on disk, cleaning up only for this particular
      :class:`bedtool`

* Universal "interval" file support

    * parent :class:`UniversalInterval` class, GFF, GTF, VCF, etc. etc. all
      subclassed from that

    * :class:`bedtool` should auto-detect the interval file format

* Support for "derived" file types

    * currently the handling of output created by ``a.intersect(b, c=True)`` 
      is special-cased, and ``a.intersect(b, wao=True)`` is unsupported by
      custom :class:`bedtool` methods.  

    * these should probably just be implemented as subclasses of the as-yet
       hypothetical :class:`UniversalInterval` class.

* Support for ``groupBy``

    * it would be nice to have support for filo_ 

* Randomization methods are still under development


Why use ``pybedtools``?
-----------------------
:mod:`pybedtools` makes working with BEDTools_ from Python code easy.

Calling BEDTools_ from Python "by hand" gets awkward for things like geting
the number of intersecting features between two bed files, for example::

    >>> # The annoying way of calling BEDTools from Python...
    >>> p1 = subprocess.Popen(['intersectBed','-a','a.bed','-b','b.bed','-u'], stdout=subprocess.PIPE)

    >>> # then pipe it to wc -l to get a line count
    >>> p2 = subprocess.Popen(['wc','-l'], stdin=subprocess.PIPE)

    >>> # finally, parse the results
    >>> results = p2.communicate()[0]
    >>> count = int(results.split()[-1])

If we wanted to get the number of features unique to :file:`a.bed` , it
would mean another 4 lines of this, with the only difference being the
``-v`` argument instead of ``-u`` for the ``intersectBed`` call..  For me,
this got old quickly, hence the creation of :mod:`pybedtools`.

Here's how to do the same thing with :mod:`pybedtools`::
    
    >>> from pybedtools import bedtool
    >>> a = bedtool('a.bed')
    >>> count = a.intersect('b.bed', u=True).count()

Behind the scenes, the :class:`pybedtools.bedtool` class does something
very similar to the subprocess example above, but in a more Python-friendly
way.  

Furthermore, for the specific case of intersections, the ``+`` and ``-``
operators have been overloaded, making many intersections extremely easy::

    >>> a = bedtool('a.bed')
    >>> b = bedtool('b.bed')
    >>> c = bedtool('c.bed')
    
    >>> (a+b).count()   # number of features in a and b
    >>> (a-b).count()   # number of features in a not b
    >>> (a+b+c).count() # number of features in a, b and c

The other BEDTools_ programs are wrapped as well, like
:meth:`bedtool.merge`, :meth:`bedtool.slop`, and others.

In addition to wrapping the BEDtools_ programs, there are many additional
:class:`bedtool` methods provided in this module that you can use in your
Python code.  

.. _limitations:

Limitations
-----------
There are some limitations you need to be aware of.  

* :mod:`pybedtools` makes heavy use of temporary files.  This makes it
  very convenient to work with, but if you are limited by disk space,
  you'll have to pay attention to this feature (see `temp principle`_ below
  for more info).

* Second, :class:`bedtool` methods that wrap BEDTools_ programs will work on
  BAM, GFF, VCF, and everything that BEDTools_ supports.  However, many
  :mod:`pybedtools`-specific methods (for example :meth:`bedtool.lengths`
  or :meth:`bedtool.size_filter`) **currently only work on BED files**.  I
  hope to add support for all interval files soon.

.. _creating a bedtool:

Creating a :class:`bedtool`
---------------------------
To create a :class:`bedtool`, first you need to import the
:mod:`pybedtools` module.  For these examples, I'm assuming you have
already done the following:

.. doctest::

    >>> import pybedtools
    >>> from pybedtools import bedtool

Next, you need a BED file to work with. If you already have one, then great
-- move on to the next section.  If not, :mod:`pybedtools` comes with some
example bed files used for testing.  You can take a look at the list of
example files that ship with :mod:`pybedtools` with the
:func:`list_example_beds` function:

.. doctest::

   >>> # list the example bed files
   >>> pybedtools.list_example_beds()
   ['a.bed', 'b.bed']

Once you decide on a file to use, feed the your choice to the
:func:`example_bed_fn` function to get the full path:

.. doctest::

   >>> # get the full path to an example bed file
   >>> bedfn = pybedtools.example_bed_fn('a.bed') 

The full path of *bedfn* will depend on your installation (this is similar
to the ``data()`` function in R_, if you're familiar with that).

Now that you have a filename -- either one of the example files or your
own, you create a new :class:`bedtool` simply by pointing it to that
filename:

.. doctest::

    >>> # create a new bedtool from the example bed file
    >>> mybedtool = bedtool(bedfn)

Alternatively, you can construct BED files from scratch by using the
``from_string`` keyword argument.  However, all spaces will be converted to
tabs using this method, so you'll have to be careful if you add "name"
columns.  This can be useful if you want to create *de novo* BED files on
the fly:

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> # an "inline" example:
    >>> fromscratch1 = pybedtools.bedtool('chrX 1 100', from_string=True)
    >>> print fromscratch1
    chrX    1   100
    <BLANKLINE>

    >>> # using a longer string to make a bed file.  Note that
    >>> # newlines don't matter, and one or more consecutive 
    >>> # spaces will be converted to a tab character.
    >>> larger_string = """
    ... chrX 1    100   feature1  0 +
    ... chrX 50   350   feature2  0 -
    ... chr2 5000 10000 another_feature 0 +
    ... """

    >>> fromscratch2 = bedtool(larger_string, from_string=True)
    >>> print fromscratch2
    chrX    1   100 feature1    0   +
    chrX    50  350 feature2    0   -
    chr2    5000    10000   another_feature 0   +
    <BLANKLINE>

Of course, you'll usually be using your own bed files that have some
biological importance for your work that are saved in places convenient for
you, for example::

    >>> a = bedtool('/data/sample1/peaks.bed')



.. _`Design principles`:

Design principles
-----------------

.. _`temp principle`:

Principle 1: Temporary files are created automatically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Let's illustrate some of the design principles behind :mod:`pybedtools` by
merging features in :file:`a.bed` that are 100 bp or less apart (*d=100*)
in a strand-specific way (*s=True*):

.. doctest::
    
    >>> from pybedtools import bedtool
    >>> import pybedtools
    >>> a = bedtool(pybedtools.example_bed_fn('a.bed'))
    >>> merged_a = a.merge(d=100, s=True)

Now *merged_a* is a :class:`bedtool` instance that contains the results of the
merge.

:class:`bedtool` objects must always point to a file on disk.  So in the
example above, *merged_a* is a :class:`bedtool`, but what file does it
point to?  You can always check the :attr:`bedtool.fn` attribute to find
out::

    >>> # what file does *merged_a* point to?
    >>> merged_a.fn
    '/tmp/pybedtools.MPPp5f.tmp'

Note that the specific filename will be different for you since it is a
randomly chosen name (handled by Python's :mod:`tempfile` module).  This
shows one important aspect of :mod:`pybedtools`: every operation results in
a new temporary file. Temporary files are stored in :file:`/tmp` by
default, and have the form :file:`/tmp/pybedtools.*.tmp`. 

Future work on :mod:`pybedtools` will focus on streamlining the temp files,
keeping only those that are needed.  For now, when you are done using the
:mod:`pybedtools` module, make sure to clean up all the temp files created
with::
    
    >>> # Deletes all tempfiles created this session.
    >>> # Don't do this yet if you're following the tutorial!
    >>> pybedtools.cleanup()

If you forget to do this, from the command line you can always do a::

    rm /tmp/pybedtools.*.tmp

to clean everything up.

If you need to specify a different directory than that used by default by
Python's tempdir_ module, then you can set it with::

    >>> pybedtools.set_tempdir('/scratch')

You'll need write permissions to this directory, and it needs to already
exist.

.. _`similarity principle`:

Principle 2: Names and arguments are as similar as possible to BEDTools_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Returning again to this example::
    
    >>> merged_a = a.merge(d=100, s=True)

This demonstrates another :mod:`pybedtools` principle: the :class:`bedtool`
methods that wrap BEDTools_ programs do the same thing and take the exact same
arguments as the BEDTools_ program.  Here we can pass *d=100* and *s=True* only
because the underlying BEDTools_ program, ``mergeBed``, can accept these
arguments.  Need to know what arguments ``mergeBed`` can take?  See the docs
for :meth:`bedtool.merge`; for more on this see :ref:`good docs principle`.

In general, remove the "Bed" from the end of the BEDTools_ program to get
the corresponding :class:`bedtool` method.  So there's a
:meth:`bedtool.subtract` method for ``subtractBed``, a
:meth:`bedtool.intersect` method for ``intersectBed``, and so on.

Since these methods just wrap BEDTools_ programs, they are as up-to-date as
the version of BEDTools_ you have installed on disk.  If you are using a
cutting-edge version of BEDTools_ that has some hypothetical argument
``-z`` for ``intersectBed``, then you can use ``a.intersectBed(z=True)``. 

.. _`default args principle`:

Principle 3: Sensible default args
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If we were running the ``mergeBed`` program from the command line, we would
would have to specify the input file with the :option:`mergeBed -i` option.

:mod:`pybedtools` assumes that if we're calling the :meth:`merge` method on
*a*, we want to operate on the bed file that *a* points to.   

In general, BEDTools_ programs that accept a single BED file as input
(by convention typically specified with the :option:`-i` option) the
default behavior for :mod:`pybedtools` is to use the :class:`bedtool`'s
file (indicated in the :attr:`bedtool.fn` attribute) as input.  

We can still pass a file using the *i* keyword argument if we wanted to be
absolutely explicit.  In fact, the following two versions produce the same
output:

.. doctest::

    >>> # The default is to use existing file for input -- no need
    >>> # to specify "i" . . .
    >>> result1 = a.merge(d=100, s=True)
    
    >>> # . . . but you can always be explicit if you'd like
    >>> result2 = a.merge(i=a.fn, d=100, s=True)

    >>> # Confirm that the output is identical
    >>> str(result1) == str(result2)
    True

Methods that have this type of default behavior are indicated by the following text in their docstring::

    .. note::
    
        For convenience, the file this bedtool object points to is passed as "-i"

There are some BEDTools_ programs that accept two BED files as input, like
``intersectBed`` where the the first file is specified with :option:`-a`
and the second file with :option:`-b`.  The default behavior for
:mod:`pybedtools` is to consider the :mod:`bedtool`'s file as ``-a`` and
the first non-keyword argument to the method as ``-b``, like this:

.. doctest::
    
    >>> result3 = a.intersect(b)
    
This is exactly the same as passing the *a* and *b* arguments explicitly:

.. doctest::

    >>> result4 = a.intersect(a=a.fn, b=b.fn)
    >>> str(result3) == str(result4)   
    True

Furthermore, the first non-keyword argument used as ``-b`` can either be a
filename *or* another :class:`bedtool` object; that is, these commands also do the same thing:

.. doctest::

   >>> result5 = a.intersect(b=b.fn)
   >>> result6 = a.intersect(b=b)
   >>> str(result5) == str(result6)
   True

Methods that accept either a filename or another :class:`bedtool` instance as their first non-keyword argument are indicated by
the following text in their docstring::

    .. note::
        
        This method accepts either a bedtool or a file name as the first
        unnamed argument




.. _`non defaults principle`:

Principal 4: Other arguments have no defaults
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Only the BEDTools_ arguments that refer to BED (or other interval) files have
defaults.  In the current version of BEDTools_, this means only the ``-i``,
``-a``, and ``-b`` arguments have defaults.  All others have no defaults
specified by :mod:`pybedtools`; they pass the buck to BEDTools programs.  This
means if you do not specify the *d* kwarg when calling :meth:`bedtool.merge`,
then it will use whatever the installed version of BEDTools_ uses for ``-d``
(currently, ``mergeBed``'s default for ``-d`` is 0).


``-d`` is an option to BEDTools_ ``mergeBed`` that accepts a value, while
``-s`` is an option that acts as a switch.  In :mod:`pybedtools`, simply
pass a value (integer, float, whatever) for value-type options like ``-d``,
and boolean values (*True* or *False*) for the switch-type options like
``-s``.

Here's another example using both types of keyword arguments; the
:class:`bedtool` object *b* (or it could be a string filename too) is
implicitly passed to ``intersectBed`` as ``-b`` (see :ref:`default args
principle` above)::

    >>> a.intersect(b, v=True, f=0.5)

Again, any option that can be passed to a BEDTools_ program can be passed
to the corresonding :class:`bedtool` method.


.. _`chaining principle`:

Principle 5: Chaining together commands
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Most methods return new :class:`bedtool` objects, allowing you to chain
things together just like piping commands together on the command line.  To
give you a flavor of this, here is how you would get the merged regions of
features shared between :file:`a.bed` (as referred to by the
:class:`bedtool` *a* we made previously) and :file:`b.bed`: (as referred to
by the :class:`bedtool` *b*):

.. doctest::
    
    >>> a.intersect(b).merge().saveas('shared_merged.bed')
    <bedtool (shared_merged.bed)>


This is equivalent to the following BEDTools_ commands::

    intersectBed -a a.bed -b b.bed | merge -i stdin > shared_merged.bed

Methods that return a new :class:`bedtool` instance are indicated with the following text in their docstring::

    .. note::
    
        This method returns a new bedtool instance

.. _`good docs principle`:

Principle 6: Check the help
~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you're unsure of whether a method uses a default, or if you want to read
about what options an underlying BEDTools_ program accepts, check the help.
Each :class:`pybedtool` method that wraps a BEDTools_ program also wraps
the BEDTools_ program help string.  There are often examples of how to use
a method in the docstring as well.

    

Example: Flanking seqs
----------------------
The :meth:`bedtool.slop` method (which calls ``slopBed``) needs a
chromosome size file.  If you specify a genome name to the
:meth:`bedtool.slop` method, it will retrieve this file for you
automatically from the UCSC Genome Browser MySQL database.

::
    
    import pybedtools
    a = pybedtools.bedtool('in.bed')
    extended = a.slop(genome='dm3',l=100,r=100)
    flanking = extended.subtract(a).saveas('flanking.bed')
    flanking.sequence(fi='dm3.fa')
    flanking.save_seqs('flanking.fa')


Or, as a one-liner::

    pybedtools.bedtool('in.bed').slop(genome='dm3',l=100,r=100).subtract(a).sequence(fi='dm3.fa').save_seqs('flanking.fa')

Don't forget to clean up!::

    pybedtools.cleanup()

Example: Region centers that are fully intergenic
--------------------------------------------------
Useful for, e.g., motif searching::
    
    a = pybedtools.bedtool('in.bed')
    
    # Sort by score
    a = a.sorted(col=5,reverse=True)

    # Exclude some regions
    a = a.subtract('regions-to-exclude.bed')

    # Get 100 bp on either side of center
    a = a.peak_centers(100).saveas('200-bp-peak-centers.bed')


Example: Histogram of feature lengths
-------------------------------------
Note that you need matplotlib installed to plot the histogram.

::

    import pylab as p
    a = pybedtools.bedtool('in.bed')
    p.hist(a.lengths(),bins=50)
    p.show()


