
.. _BEDTools: http://github.com/arq5x/bedtools

.. _filo: https://github.com/arq5x/filo

.. _R: http://www.r-project.org/

.. _topical:

Topical documentation
=====================

Release notes
-------------

0.2.3dev
~~~~~~~~
* Added history mechanism -- see :ref:`working with history`.
* Added tagging mechanism

Future work
-----------
The following is a list of items I plan to work on for future releases of
:mod:`pybedtools`.  Help and feedback are welcome!:

* Finish wrapping the rest of the BEDTools_ programs

* **DONE:** Better mechanism for handling temp files.  

    * **DONE:** :class:`bedtool` objects could keep track of all their 'parent'
      tempfiles, and as such would retain the history of their creation.

    * **N/A:** indexing into the history of a :class:`bedtool` would give you access
      to these previous files.  *not implemented -- complex histories would
      require a tree, which is not easily indexed.  Instead I'm using a "tag"
      system*
    
    * **DONE:** could have a :meth:`bedtool.delete_history`, which would delete all the
      intermediate tempfiles on disk, cleaning up only for this particular
      :class:`bedtool`

* Universal "interval" file support

    * Should have a parent :class:`UniversalInterval` class and then have
      GFF, GTF, VCF, etc. etc. all subclassed from that

    * :class:`bedtool` should auto-detect the interval file format

* Support for "derived" file types

    * currently the handling of output created by ``a.intersect(b, c=True)`` 
      is special-cased, and ``a.intersect(b, wao=True)`` is unsupported by
      custom :class:`bedtool` methods.  

    * these should probably just be implemented as subclasses of the as-yet
       hypothetical :class:`UniversalInterval` class.

* Support for ``groupBy``

    * it would be nice to have support for filo_ 

* Randomization methods still under development

* Have wrappers use subprocess.Popen instead of os.system calls for better
  trapping of errors

* More universal wrapper -- perhaps as decorator?

* Full unit tests to bolster the tutorial doctests

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

* :class:`bedtool` methods that wrap BEDTools_ programs
  (:meth:`bedtool.intersect`, :meth:`bedtool.merge`,
  :meth:`bedtool.subtract`, etc) will work on BAM, GFF, VCF, and everything
  that BEDTools_ supports.  However, many :mod:`pybedtools`-specific
  methods (for example :meth:`bedtool.lengths` or
  :meth:`bedtool.size_filter`) **currently only work on BED files**.  I
  hope to add support for all interval files soon.

* **Not all BEDTools programs are wrapped** -- this package is still a work
  in progress.  I wrapped the ones I use most often and still need to wrap
  the others. The following table shows what's currently wrapped:

=================   ============================ 
BEDTools program    :class:`bedtool` method name
=================   ============================ 
``intersectBed``    :meth:`bedtool.intersect`
``subtractBed``     :meth:`bedtool.subtract`
``fastaFromBed``    :meth:`bedtool.sequence`
``slopBed``         :meth:`bedtool.slop`
``windowBed``       :meth:`bedtool.window`
``closestBed``      :meth:`bedtool.closest`
``shuffleBed``      :meth:`bedtool.shuffle`
=================   ============================ 


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
:func:`list_example_files` function:

.. doctest::

   >>> # list the example bed files
   >>> pybedtools.list_example_files()
   ['a.bed', 'b.bed']

Once you decide on a file to use, feed the your choice to the
:func:`example_filename` function to get the full path:

.. doctest::

   >>> # get the full path to an example bed file
   >>> bedfn = pybedtools.example_filename('a.bed') 

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
    >>> a = bedtool(pybedtools.example_filename('a.bed'))
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
    
    >>> b = pybedtools.bedtool(pybedtools.example_filename('b.bed'))    
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


