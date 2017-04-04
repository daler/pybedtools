.. include:: includeme.rst

.. _`Design principles`:

Design principles
-----------------
Hopefully, understanding (or just being aware of) these design principles
will help in getting the most out of :mod:`pybedtools` and working
efficiently.

.. _`temp principle`:

Principle 1: Temporary files are created (and deleted) automatically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Using :class:`BedTool` instances typically has the side effect of creating
temporary files on disk.  Even when using the iterator protocol of
:class:`BedTool` objects, temporary files may be created in order to run
BEDTools programs (see :ref:`BedTools as iterators` for more on this latter topic).

Let's illustrate some of the design principles behind :mod:`pybedtools` by
merging features in :file:`a.bed` that are 100 bp or less apart (`d=100`)
in a strand-specific way (`s=True`):

.. doctest::

    >>> from pybedtools import BedTool
    >>> import pybedtools
    >>> a = BedTool(pybedtools.example_filename('a.bed'))
    >>> merged_a = a.merge(d=100, s=True)

Now `merged_a` is a :class:`BedTool` instance that contains the results of the
merge.

:class:`BedTool` objects must always point to a file on disk.  So in the
example above, `merged_a` is a :class:`BedTool`, but what file does it
point to?  You can always check the :attr:`BedTool.fn` attribute to find
out::

    >>> # what file does `merged_a` point to?
    >>> merged_a.fn
    '/tmp/pybedtools.MPPp5f.tmp'

Note that the specific filename will be different for you since it is a
randomly chosen name (handled by Python's :mod:`tempfile` module).  This
shows one important aspect of :mod:`pybedtools`: every operation results in
a new temporary file. Temporary files are stored in :file:`/tmp` by
default, and have the form :file:`/tmp/pybedtools.*.tmp`. 

By default, at exit all temp files created during the session will be deleted.
However, if Python does not exit cleanly (e.g., from a bug in client code),
then the temp files will not be deleted.

If this happens, from the command line you can always do a::

    rm /tmp/pybedtools.*.tmp

In the middle of a session, you can force a deletion of all tempfiles created thus far::

    >>> # Don't do this yet if you're following the tutorial!
    >>> pybedtools.cleanup()


Alternatively, in this session or another session you can use::

    >>> pybedtools.cleanup(remove_all=True)

to remove all files that match the pattern
:file:`<tempdir>/pybedtools.*.tmp` where `<tempdir>` is the current value
of `pybedtools.get_tempdir()`.

If you need to specify a different directory than that used by default by
Python's tempdir_ module, then you can set it with::

    >>> pybedtools.set_tempdir('/scratch')

You'll need write permissions to this directory, and it needs to already
exist.  All temp files will then be written to that directory, until the
tempdir is changed again.

.. _`similarity principle`:

Principle 2: Names and arguments are as similar as possible to BEDTools_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
As much as possible, BEDTools programs and :class:`BedTool` methods share
the same names and arguments.

Returning again to this example::

    >>> merged_a = a.merge(d=100, s=True)

This demonstrates that the :class:`BedTool` methods that wrap BEDTools_
programs do the same thing and take the exact same arguments as the
BEDTools_ program.  Here we can pass `d=100` and `s=True` only because the
underlying BEDTools_ program, `mergeBed`, can accept these arguments.
Need to know what arguments `mergeBed` can take?  See the docs for
:meth:`BedTool.merge`; for more on this see :ref:`good docs principle`.

In general, remove the "Bed" from the end of the BEDTools_ program to get
the corresponding :class:`BedTool` method.  So there's a
:meth:`BedTool.subtract` method for `subtractBed`, a
:meth:`BedTool.intersect` method for `intersectBed`, and so on.

.. _`version principle`:

Principle 3: Indifference to BEDTools version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Since :class:`BedTool` methods just wrap BEDTools_ programs, they are as up-to-date as
the version of BEDTools_ you have installed on disk.  If you are using a
cutting-edge version of BEDTools_ that has some hypothetical argument
`-z` for `intersectBed`, then you can use `a.intersectBed(z=True)`.

:mod:`pybedtools` will also raise an exception if you try to use a method
that relies on a more recent version of BEDTools than you have installed.


.. _`default args principle`:

Principle 4: Sensible default args
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If we were running the ``mergeBed`` program from the command line, we
would have to specify the input file with the :option:`mergeBed -i` option.

:mod:`pybedtools` assumes that if we're calling the :meth:`merge` method on
the :class:`BedTool`, `a`, we want to operate on the bed file that `a`
points to.

In general, BEDTools_ programs that accept a single BED file as input
(by convention typically specified with the :option:`-i` option) the
default behavior for :mod:`pybedtools` is to use the :class:`BedTool`'s
file (indicated in the :attr:`BedTool.fn` attribute) as input.  

We can still pass a file using the `i` keyword argument if we wanted to be
absolutely explicit.  In fact, the following two versions produce the same
output:

.. doctest::

    >>> # The default is to use existing file for input -- no need
    >>> # to specify "i" . . .
    >>> result1 = a.merge(d=100, s=True)

    >>> # . . . but you can always be explicit if you'd like
    >>> result2 = a.merge(i=a.fn, d=100, s=True)

    >>> # Confirm that the output is identical
    >>> result1 == result2
    True

Methods that have this type of default behavior are indicated by the following text in their docstring::

    .. note::

        For convenience, the file this BedTool object points to is passed as "-i"

There are some BEDTools_ programs that accept two BED files as input, like
``intersectBed`` where the the first file is specified with `-a` and the
second file with `-b`.  The default behavior for :mod:`pybedtools` is to
consider the :mod:`BedTool`'s file as `-a` and the first non-keyword
argument to the method as `-b`, like this:

.. doctest::

    >>> b = pybedtools.example_bedtool('b.bed')
    >>> result3 = a.intersect(b)

This is exactly the same as passing the `a` and `b` arguments explicitly:

.. doctest::

    >>> result4 = a.intersect(a=a.fn, b=b.fn)
    >>> result3 == result4
    True

Furthermore, the first non-keyword argument used as `-b` can either be a
filename *or* another :class:`BedTool` object; that is, these commands also
do the same thing:

.. doctest::

   >>> result5 = a.intersect(b=b.fn)
   >>> result6 = a.intersect(b=b)
   >>> str(result5) == str(result6)
   True

Methods that accept either a filename or another :class:`BedTool` instance
as their first non-keyword argument are indicated by the following text in
their docstring::

    .. note::

        This method accepts either a BedTool or a file name as the first
        unnamed argument

.. _`non defaults principle`:

Principal 5: Other arguments have no defaults
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Only the BEDTools_ arguments that refer to BED (or other interval) files have
defaults.  In the current version of BEDTools_, this means only the `-i`,
`-a`, and `-b` arguments have defaults.  All others have no defaults
specified by :mod:`pybedtools`; they pass the buck to BEDTools programs.  This
means if you do not specify the `d` kwarg when calling :meth:`BedTool.merge`,
then it will use whatever the installed version of BEDTools_ uses for `-d`
(currently, `mergeBed`'s default for `-d` is 0).


`-d` is an option to BEDTools_ `mergeBed` that accepts a value, while
`-s` is an option that acts as a switch.  In :mod:`pybedtools`, simply
pass a value (integer, float, whatever) for value-type options like `-d`,
and boolean values (`True` or `False`) for the switch-type options like
`-s`.

Here's another example using both types of keyword arguments; the
:class:`BedTool` object `b` (or it could be a string filename too) is
implicitly passed to `intersectBed` as `-b` (see :ref:`default args
principle` above)::

    >>> a.intersect(b, v=True, f=0.5)

Again, any option that can be passed to a BEDTools_ program can be passed
to the corresonding :class:`BedTool` method.


.. _`chaining principle`:

Principle 6: Chaining together commands
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Most methods return new :class:`BedTool` objects, allowing you to chain
things together just like piping commands together on the command line.  To
give you a flavor of this, here is how you would get the merged regions of
features shared between :file:`a.bed` (as referred to by the
:class:`BedTool` `a` we made previously) and :file:`b.bed`: (as referred to
by the :class:`BedTool` `b`):

.. doctest::
    
    >>> a.intersect(b).merge().saveas('shared_merged.bed')
    <BedTool(shared_merged.bed)>


This is equivalent to the following BEDTools_ commands::

    intersectBed -a a.bed -b b.bed | merge -i stdin > shared_merged.bed

Methods that return a new :class:`BedTool` instance are indicated with the following text in their docstring::

    .. note::
    
        This method returns a new BedTool instance

.. _`good docs principle`:

Principle 7: Check the help
~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you're unsure of whether a method uses a default, or if you want to read
about what options an underlying BEDTools_ program accepts, check the help.
Each :class:`pyBedTool` method that wraps a BEDTools_ program also wraps
the BEDTools_ program help string.  There are often examples of how to use
a method in the docstring as well. The documentation is also run through
doctests, so the code you read here is guaranteed to work and be
up-to-date.
