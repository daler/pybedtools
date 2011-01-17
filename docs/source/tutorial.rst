
.. currentmodule:: pybedtools

.. |bt| replace:: BEDTools_
.. _BEDTools: http://github.com/arq5x/bedtools

Tutorial
========

Why use ``pybedtools``?
-----------------------
This module allows you to use |bt| directly from your Python code
without awkward system calls.  The provided :class:`pybedtools.bedtool` class
wraps the BEDtools command line programs in an intuitive and easy-to-use
interface.  As a quick illustration of the streamlining possible, here's
how to get the number of features shared between a.bed and b.bed, those
unique to :file:`a.bed`, and those unique to :file:`b.bed`::

    from pybedtools import bedtool
    a = bedtool('a.bed')
    b = bedtool('b.bed')
    (a+b).count()    # shared in a and b
    (a-b).count()    # unique to a
    (b-a).count()    # unique to b

In contrast, here's how you'd do the same from the command line:: 

    intersectBed -a a.bed -b b.bed -u | wc -l   # shared in a and b
    intersectBed -a a.bed -b b.bed -v | wc -l   # unique to a
    intersectBed -a b.bed -b a.bed -v | wc -l   # unique to b

To do the same thing in Python, *each* of these lines would have to be
wrapped in awkward, piped :func:`subprocess.Popen` calls::
    
    p1 = subprocess.Popen(['intersectBed','-a','a.bed','-b','b.bed','-u'], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['wc','-l'], stdin=subprocess.PIPE)
    results = p2.communicate()[0]
    count = results.split()[-1]

Behind the scenes, the :class:`pybedtools.bedtool` class does something very
similar to this -- but conveniently makes the functionality available as
``a+b``.

In addition to wrapping the BEDtools programs, there are many additional
:class:`bedtool` methods provided in this module that you can use in your
Python code.

Creating a :class:`pybedtools.bedtool`
--------------------------------------
First you need to import the :mod:`pybedtools` module.  For the rest of the
tutorial, I'm assuming you have already done the following:

.. doctest::

    >>> import pybedtools
    >>> from pybedtools import bedtool

Next, you need a BED file to work with.  Luckily, :mod:`pybedtools` comes
with some example bed files.  You can take a look at the list of example
files that ship with :mod:`pybedtools` with the :func:`list_example_beds`
function:

.. doctest::

   >>> # list the example bed files
   >>> pybedtools.list_example_beds()
   ['a.bed']

Once you decide on a file to use, feed the your choice to the
:func:`example_bed` function to get the full path:

.. doctest::

   >>> # get the full path to an example bed file
   >>> bedfn = pybedtools.example_bed('a.bed') 

The full path of *bedfn* will depend on your installation. 

Now that you have a filename -- either one of the example files or your
own, you create a new :class:`bedtool` simply by pointing it to that
filename:

.. doctest::

    >>> # create a new bedtool from the example bed file
    >>> mybedtool = bedtool(bedfn)

Alternatively, you can construct BED files from scratch by using the
``from_string`` keyword argument.  However, all spaces will be converted to
tabs using this method, so you'll have to be careful if you add "name"
columns:

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
biological importance for your work.  But for the purposes of this
tutorial, we'll be using two bed files that you can get from the examples
directory:

.. doctest::
    
    >>> a = bedtool(pybedtools.example_bed('a.bed'))
    >>> b = bedtool(pybedtools.example_bed('b.bed'))


Design principles: an example
-----------------------------
Let's illustrate some of
the design principles behind :mod:`pybedtools` by merging features in
:file:`a.bed` that are 100 bp or less apart (*d=100*) in a strand-specific
way (*s=True*):

.. doctest::
    
    >>> merged_a = a.merge(d=100, s=True)

Now *merged_a* is a :class:`bedtool` instance that contains the results of the
merge.

Principle 1: temporarly files are created automatically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
First, just as BEDTools_ is primarily file-based, :class:`bedtools` always
point to a file on disk.  So what file does *merged_a* point to?  You can
check the :attr:`bedtool.fn` attribute to find out::

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


Principle 2: Names and arguments are similar to BEDTools_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Second, the :meth:`bedtool.merge` method does the same thing as the
BEDTools_ program ``mergeBed``.  In general, remove the "Bed" from the end
of the BEDTools_ program to get the corresponding :class:`bedtool` method.
So there's a :meth:`bedtool.subtract` method, a :meth:`bedtool.intersect`
method, and so on.

Note that we used the *d=100* keyword argument ("kwarg").  Since
:option:`mergeBed -d` is an option for the BEDTools_ ``mergeBed`` program,
it can also be used by :meth:`bedtool.merge`.  The same goes for any other
options.

Principle 3: Sensible default args
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When running BEDTools_ ``mergeBed`` program from the command line, you
would have to specify the input file with the :option:`mergeBed -i` option.
In fact, the following two methods produce the same output:

.. doctest::

    >>> # The default is to use existing file for input -- no need
    >>> # to specify "i" . . .
    >>> result1 = a.merge(d=100, s=True)
    
    >>> # . . . but you can always be explicit if you'd like
    >>> result2 = a.merge(i=a.fn, d=100, s=True)

    >>> # Confirm that the output is identical
    >>> str(result1) == str(result2)
    True
    
For BEDTools_ programs that accept a single BED file as input (typically
specified with the :option:`-i` option), the default for :mod:`pybedtools`
is to use the :class:`bedtool`'s file as input.  

For BEDTools_ programs that accept two BED files as input (like
``intersectBed``, with the first file as :option:`-a` and the second file
as :option:`-b`), the default for :mod:`pybedtools` is to consider the
:mod:`bedtool`'s file as "a" and the first non-keyword argument as "b".

Principal 4: Other arguments have no defaults
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*-d* is an option to BEDTools_ ``mergeBed`` that accepts a value, while
*-s* is an option that acts as a switch.  In :mod:`pybedtools`, simply pass
boolean values (:keyword:`True` or :keyword:`False`) for the switch-type
options, and pass a value for the value-type options.

Other than the :option:`-i`, :option:`-a`, and :option:`-b` options
mentioned above, these other options have no defaults.

Principle 5: Check the help
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


REWRITE ENDED HERE
------------------


General usage
-------------
.. doctest::

    >>> from pybedtools import bedtool
    >>> a = '''
    ...         chrX 1   100
    ...         chrX 200 500
    ...         chrY 499 600
    ... '''
    >>> b = '''
    ...         chrX 10  60
    ...         chrY 200 500
    ... '''
    >>> a = bedtool(a, from_string=True).saveas('a.bed')
    >>> b = bedtool(b, from_string=True).saveas('b.bed')

Arguments are the same as BEDtools command line programs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The methods in the ``bedtool`` class mimic the command line usage of
``BEDtools`` as closely as possible.  Any flag that can be passed to the
``BEDtools`` programs can be passed as a keyword argument in the
corresponding ``pybedtools.bedtool`` method.

On/off switches (e.g., ``-c``, ``-u``, or ``-v`` for ``intersectBed``) are
called with a boolean kwarg; others (like ``-f`` for ``intersectBed``) are
passed in like a normal kwarg with appropriate values.  For example::

    a = bedtool('in.bed')
    a.intersect('other.bed', v=True, f=0.5)

Typically, for convenience ``-i``, ``-a``, and ``-b`` are already
passed for you although you can override this by passing these keyword
arguments explicitly.  The second line above could have equivalently been
called as::

    a.intersect(a='in.bed', b='other.bed', v=True, f=0.5)

Conveniently, the docstring for a method automatically includes the help
text of the original ``BEDtools`` program, so you can check which kwargs
you can use directly from the interpreter.

Typical workflow includes temporary files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Typical workflow is to set up a bedtool object with a bed file you already have::

    a = bedtool('in.bed')
    
`a` now references the file ``in.bed`` on disk.  

Using BEDtools from the command line, we might use the following in order
to get a new bed file of the intersection of this file with another bed
file, ``other.bed`` but only returning uniquely intersecting features from
``in.bed``::

    intersectBed -a in.bed -b other.bed -u > intersection.bed

Using ``pybedtools``::

    b = a.intersect('other.bed', u=True)

This creates a new temp file in ``/tmp`` by default, but you can change
where temp files are saved using ``pybedtools.set_tempdir()``.  To save a
file explicitly and optionally add a trackline that will label it in a
genome browser, use ``saveas()``::
    
    b.saveas('intersection.bed',trackline='track name="intersection"')

Cleaning up temporary files
~~~~~~~~~~~~~~~~~~~~~~~~~~~
When you're done, it's a good idea to clean up temporary files.  Temp files
are never deleted until you explicitly say so::

    pybedtools.cleanup()

If you forget to call ``pybedtools.cleanup()``, you can always manually
delete the files from your temp dir (typically ``/tmp``), though you can
specify this with ``pybedtools.set_tempdir()``.  They are the files that
follow the pattern ``pybedtools.*.tmp``.

Chaining together commands
~~~~~~~~~~~~~~~~~~~~~~~~~~
Most methods return new ``bedtool`` objects, allowing you to chain things
together.  This means that you can chain commands together, just like
piping things together on the command line.  To give you a flavor of this,
here's how you would get 10 random centers of features that are unique to
the file ``other.bed``::

    b = a.intersect('other.bed',v=True).feature_centers(100).random_subset(10)
 

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

