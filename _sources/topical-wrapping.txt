Wrapping new tools
==================
This section serves as a reference for wrapping new tools as they are added to
BEDTools.


Example program description
---------------------------
Let's assume we would like to wrap a new program, appropriately named
`newProgramBed`.  Its signature from the command line is `newProgramBed -a
<infile> -b <other file> [options]`, and it accepts `-a stdin` to indicate
data is being piped to it::

    newProgramBed -a <BED/VCF/GFF> -b <BED/VCF/GFF> [options]


Method name
-----------
Generally, I've tried to keep method names as similar as possible to
BEDTools programs while still being PEP8-compliant.  The trailing 'Bed' is
usually removed from the program name. So here the name would probably be
`new_program`.


Define a method in :class:`BedTool`
-----------------------------------
Define a method in :class:`BedTool` . . . and *don't add any content to the
function body*.  This is because the decorator we're about to add will
replace the method wholesale; anything that's in the function body will
effectively be ignored.

::

    def new_program(self):
        pass


Add the :func:`_wraps` decorator
--------------------------------
This is where most of the work happens.

Since most of the work of wrapping BEDTools programs needs to happen every
time a new program is wrapped, this work is abstracted out into the
:func:`_wraps` function.

.. note::

    The :func:`_wraps` docstring and source is the best place to learn the
    details on what it's doing; here we'll focus on using it.

Our hypothetical program, `newProgramBed`, takes `-a` as the first input.
We'd like to have `-a` implicitly be passed as whatever our
:class:`BedTool` already points to, so we use the `implicit='a'` kwarg to
:func:`_wraps` here.  `newProgramBed` also takes a second input, `-b`.  We
describe that to the wrapper with the `other='b'` kwarg.

Any other keyword args that are used when calling the method will
automatically be passed to the program.  So if `newProgramBed` has an
optional `-s` argument, we don't need to specify that here.  When the user
passes an `s=True` kwarg, it will be passed automatically to
`newProgramBed` as the argument `-s`.  If `newProgramBed` does not accept a
`-z` argument but the user passes one anyway, we rely on the BEDTools
program to do the error-checking of arguments and report any errors back to
Python.

Here's what the new method looks like so far:

::

    @_wraps(prog='newProgramBed', implicit='a', other='b')
    def new_program(self):
        pass

For wrapped programs that expect a genome file or have more complex
arguments, see the docstring and source for :func:`_wrap`.


Add doctests
------------
While the function body will be replaced wholesale by the decorator, the
docstring will be copied to the new function.  This is important because it
means we can write meaningful documentation and, even more importantly,
doctests for this method.  Writing a doctest within the method's docstring
means it will automatically be found by the test suite.

::

    @_wraps(prog='newProgramBed', implicit='a', other='b')
    def new_program(self):
        """
        Converts all features to length of 1.

        Example usage:

        >>> a = pybedtools.example_bedtool('a.bed')
        >>> b = pybedtools.example_bedtool('b.bed')
        >>> c = a.new_program(b, s=True)
        >>> print c  #+NORMALIZE_WHITESPACE
        chr1	1	2
        chr1	100	101
        chr1	150	151
        chr1	900	901
        <BLANKLINE>
        """


Add to list of known programs
-----------------------------
The last thing to do is to add the new program to the end of the tuple
`pybedtools.helpers._prog_names`.  This creates rudimentary security by only
allowing these programs to be called, and acts as sort of a central registry
for programs that have been wrapped.

Summary
-------
That's it!  We now have a method, :meth:`BedTool.new_program`, that wraps
a hypothetical `newProgramBed` BEDTools program, will accept any optional
args that `newProgramBed` does, will return a new :class:`BedTool`
containing the results, *and it's tested*.

This new method can be be chained with other :class:`BedTool` instances,
used as an iterator or generator, or anything else a normal
:class:`BedTool` can do . . . for example::

    a = pybedtools.example_bed('a.bed')
    b = pybedtools.example_bed('b.bed')
    c = a.new_program(b, s=True).filter(lambda x: x.start < 125).saveas('t.bed', trackline='track name="one-bp features"')

.. _decorator: http://www.python.org/dev/peps/pep-0318/
