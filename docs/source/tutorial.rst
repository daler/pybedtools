
.. currentmodule:: pybedtools

.. _R: http://www.r-project.org/

.. _BEDTools: http://github.com/arq5x/bedtools

.. _Tutorial:

.. _BEDTools documentation: http://code.google.com/p/bedtools/#Documentation

.. _Learn Python the Hard Way: http://learnpythonthehardway.org/static/LearnPythonTheHardWay.pdf

.. _IPython: http://ipython.scipy.org/moin/

Tutorial
========
This tutorial assumes that 

1. You know how to use BEDTools_ and
2. You know how to use Python

If this isn't the case, then 1) check out the BEDTools `BEDTools
documentation`_ and 2) a good Python tutorial like `Learn Python the Hard
Way`_.  The ">>>" means to type the code into an interactive Python interpreter like IPython_ , or optionally, a script.

First, import the :mod:`pybedtools` module:

.. doctest::

    >>> import pybedtools

Then set up a :class:`bedtool` instance.  This can be a file that you
already have, or one of the example files.  Currently, only BED format
files are supported -- see :ref:`Limitations`; see :ref:`Creating a
bedtool` for more information on creating :class:`bedtool` objects in
general.

.. doctest::
    
    >>> # get an example file that ships with pybedtools
    >>> # (or use your own BED file)
    >>> filename = pybedtools.example_files('a.bed')
    
    >>> a = pybedtools.bedtool(filename)

Set up a second one so we can do intersections and subtractions:

.. doctest::

    >>> b = pybedtools.bedtool(pybedtools.example_files('b.bed')

Intersect *a* with *b*, reporting only unique hits:

.. doctest::

    >>> in_both = a.intersect(b, u=True)

The :meth:`bedtool.intersect` method wraps the BEDTools_ program
``intersectBed``.  Note that the *u=True* kwarg has the same effect as
using the ``-u`` switch with ``intersectBed``.  The same idea holds true for
all programs wrapped by a :class:`bedtool` object -- keyword arguments
match command line switches.




Topical documentation
=====================


Why use ``pybedtools``?
-----------------------
:mod:`pybedtools` makes working with BEDTools_ from Python code easy.

Calling BEDTools_ from Python "by hand" gets awkward for things like geting
the number of intersecting features between two bed files::

    >>> # The annoying way of calling BEDTools from Python...
    >>> p1 = subprocess.Popen(['intersectBed','-a','a.bed','-b','b.bed','-u'], stdout=subprocess.PIPE)

    >>> # pipe it to wc -l to get a line count
    >>> p2 = subprocess.Popen(['wc','-l'], stdin=subprocess.PIPE)

    >>> # parse the results
    >>> results = p2.communicate()[0]
    >>> count = int(results.split()[-1])

If we wanted to get the number of features unique to :file:`a.bed` and not
:file:`b.bed` , it would mean another 4 lines of this.  For me, this got
old quickly, hence the creation of :mod:`pybedtools`.

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
  you'll have to pay attention to this feature (see `principle 1`_ below
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
:mod:`pybedtools` module.  For the rest of the tutorial, I'm assuming you
have already done the following:

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
   ['a.bed']

Once you decide on a file to use, feed the your choice to the
:func:`example_bed` function to get the full path:

.. doctest::

   >>> # get the full path to an example bed file
   >>> bedfn = pybedtools.example_bed('a.bed') 

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
    chrX	1	100
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
    chrX	1	100	feature1	0	+
    chrX	50	350	feature2	0	-
    chr2	5000	10000	another_feature	0	+
    <BLANKLINE>

Of course, you'll usually be using your own bed files that have some
biological importance for your work that are saved in places convenient for
you, for example::

    >>> a = bedtool('/data/sample1/peaks.bed')

But for the purposes of this tutorial, we'll be using two bed files that
you can get from the examples directory:

.. doctest::
    
    >>> a = bedtool(pybedtools.example_bed('a.bed'))
    >>> b = bedtool(pybedtools.example_bed('b.bed'))


.. _`Design principles`:

Design principles: an example
-----------------------------

.. _`principle 1`:

Principle 1: temporarly files are created automatically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Let's illustrate some of the design principles behind :mod:`pybedtools` by
merging features in :file:`a.bed` that are 100 bp or less apart (*d=100*)
in a strand-specific way (*s=True*):

.. doctest::
    
    >>> from pybedtools import bedtool
    >>> import pybedtools
    >>> a = bedtool(pybedtools.example_bed('a.bed')
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


Principle 2: Names and arguments are similar to BEDTools_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Returning again to this example::
    
    >>> merged_a = a.merge(d=100, s=True)

The second principle is that the :meth:`bedtool.merge` method does the same
thing and takes the same arguments as the BEDTools_ program ``mergeBed``.
In general, remove the "Bed" from the end of the BEDTools_ program to get
the corresponding :class:`bedtool` method.  So there's a
:meth:`bedtool.subtract` method, a :meth:`bedtool.intersect` method, and so
on.

Note that we used the *d=100* keyword argument.  Since :option:`mergeBed
-d` is an option for the BEDTools_ ``mergeBed`` program, it can also be
used by :meth:`bedtool.merge`.  The same goes for any other options.

Principle 3: Sensible default args
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When running BEDTools_ ``mergeBed`` program from the command line, you
would have to specify the input file with the :option:`mergeBed -i` option.

:mod:`pybedtools` assumes that if you're calling the :meth:`merge` method
on *a*, you want to operate on the bed file that *a* points to.   

In general, BEDTools_ programs that accept a single BED file as input
(typically specified with the :option:`-i` option) the default for
:mod:`pybedtools` is to use the :class:`bedtool`'s file as input.  

For BEDTools_ programs that accept two BED files as input (like
``intersectBed``, with the first file as :option:`-a` and the second file
as :option:`-b`), the default for :mod:`pybedtools` is to consider the
:mod:`bedtool`'s file as "a" and the first non-keyword argument as "b".
Furthermore, the first non-keyword argument can either be a filename *or*
another :class:`bedtool` object.

You can still pass a file in using the *i* keyword argument, if you want.
In fact, the following two versions produce the same output:

.. doctest::

    >>> # The default is to use existing file for input -- no need
    >>> # to specify "i" . . .
    >>> result1 = a.merge(d=100, s=True)
    
    >>> # . . . but you can always be explicit if you'd like
    >>> result2 = a.merge(i=a.fn, d=100, s=True)

    >>> # Confirm that the output is identical
    >>> str(result1) == str(result2)
    True
    


Principal 4: Other arguments have no defaults
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*-d* is an option to BEDTools_ ``mergeBed`` that accepts a value, while
*-s* is an option that acts as a switch.  In :mod:`pybedtools`, simply pass
boolean values (:keyword:`True` or :keyword:`False`) for the switch-type
options, and pass a value for the value-type options.

On/off switches (e.g., :option:`-c`, :option:`-u`, or :option:`v` for ``intersectBed``) are
called with a boolean kwarg; others (like :option:`-f` for ``intersectBed``) are
passed in like a normal kwarg with appropriate values.  As another
example::

    >>> a.intersect(b, v=True, f=0.5)

Other than the :option:`-i`, :option:`-a`, and :option:`-b` options for
input files mentioned above, these other options like :option:`-d` and
:option:`s` have no defaults.

Again, any option that can be passed to a BEDTools_ program can be passed
to the corresonding :class:`bedtool` method.

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

This is equivalent to the following BEDTools_ commands::

    intersectBed -a a.bed -b b.bed | merge -i stdin > shared_merged.bed

Principle 6: Check the help
~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you're unsure of whether a method uses a default, or if you want to read
about what options an underlying BEDTools_ program accepts, check the help.
Each :class:`pybedtool` method that wraps a BEDTools_ program also wraps
the BEDTools_ program help string.  There are often examples of how to use
a method in the docstring as well:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> # Checking the help!
    >>> help (a.merge)
    Help on method merge in module pybedtools.bedtool:
    <BLANKLINE>
    merge(self, **kwargs) method of pybedtools.bedtool.bedtool instance
        *pybedtools help:*
    <BLANKLINE>
                Merge overlapping features together. Returns a new bedtool object.
    <BLANKLINE>
                Example usage::
    <BLANKLINE>
                    a = bedtool('in.bed')
    <BLANKLINE>
                    # allow merging of features 100 bp apart
                    b = a.merge(d=100)
    <BLANKLINE>
    <BLANKLINE>
    <BLANKLINE>
        *Original BEDtools program help:*
    <BLANKLINE>
                Program: mergeBed (v2.10.0)
                Author:  Aaron Quinlan (aaronquinlan@gmail.com)
                Summary: Merges overlapping BED/GFF/VCF entries into a single interval.
    <BLANKLINE>
                Usage:   mergeBed [OPTIONS] -i <bed/gff/vcf>
    <BLANKLINE>
                Options: 
                        -s      Force strandedness.  That is, only merge features
                                that are the same strand.
                                - By default, merging is done without respect to strand.
    <BLANKLINE>
                        -n      Report the number of BED entries that were merged.
                                - Note: "1" is reported if no merging occurred.
    <BLANKLINE>
                        -d      Maximum distance between features allowed for features
                                to be merged.
                                - Def. 0. That is, overlapping & book-ended features are merged.
                                - (INTEGER)
    <BLANKLINE>
                        -nms    Report the names of the merged features separated by semicolons.
    <BLANKLINE>



Saving, printing, viewing
-------------------------

Let's make two files from scratch:

.. doctest::

    >>> from pybedtools import bedtool
    >>> string1 = '''
    ...           chrX 1   100
    ...           chrX 200 500
    ...           chrY 499 600
    ...           '''
    >>> string2 = '''
    ...           chrX 10  60
    ...           chrY 200 500
    ...           '''
    >>> x = bedtool(string1, from_string=True)
    >>> y = bedtool(string2, from_string=True)

We can save a new file, adding a track line as well (a newline is added for
you if you forget):

.. doctest::

    >>> y.saveas('y.bed', trackline="track name='example bed' color=135,0,85")

View the first 10 lines of the bed file:

.. doctest::

    >>> y.head()

Print the *whole* file:

.. doctest::

    >>> print y


Individual methods
------------------

:meth:`bedtool.intersect`
~~~~~~~~~~~~~~~~~~~~~~~~~

Examples of how to use :meth:`bedtool.intersect`.


First, we set up the files to use:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> # Set up the bedtools
    >>> a = bedtool('a.bed')
    >>> b = bedtool('b.bed')

Here's what they look like:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print a
    chrX	1	100
    chrX	200	500
    chrY    499 600

    >>> print b 
    chrX	10	60
    chrY	200	500

Default is to only report the part that intersects:


.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print a.intersect(b)
    chrX	10	60
    chrY    499 500

Use the `-u` flag of ``intersectBed`` to report full features in :file:`a.bed`
that intersected :file:`b.bed`:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print a.intersect(b, u=True)
    chrX   1    100
    chrY   499  600

Or the opposite -- features in :file:`a.bed` that are not in :file:`b.bed`:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print a.intersect(b, v=True)
    chrX	200	500

Report features in :file:`a.bed`, attaching an additional column indicating how many
features in :file:`b.bed` intersected:

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> print a.intersect(b, c=True)
    chrX	1	100	1
    chrX	200	500	0
    chrY	499	600	1

You can retrieve these counts later using the :meth:`bedtool.counts` method:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> result = a.intersect(b, c=True)
    >>> print result.counts()
    [1, 0, 1]

    

Example: Flanking seqs
----------------------
The ``slop()`` method (which calls ``slopBed``) needs a chromosome size
file.  If you specify a genome name to the ``slop()`` method, it will
retrieve this file for you automatically from the UCSC Genome Browser MySQL
database.

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


