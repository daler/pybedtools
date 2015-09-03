.. include:: includeme.rst

.. _creating a BedTool:

Creating a :class:`BedTool`
---------------------------
To create a :class:`BedTool`, first you need to import the
:mod:`pybedtools` module.  For these examples, I'm assuming you have
already done the following:

.. doctest::

    >>> import pybedtools
    >>> from pybedtools import BedTool

Next, you need a BED file to work with. If you already have one, then great
-- move on to the next section.  If not, :mod:`pybedtools` comes with some
example bed files used for testing.  You can take a look at the list of
example files that ship with :mod:`pybedtools` with the
:func:`list_example_files` function:

.. doctest::

   >>> # list the example bed files
   >>> files = pybedtools.list_example_files()

Once you decide on a file to use, feed the your choice to the
:func:`example_filename` function to get the full path:

.. doctest::

   >>> # get the full path to an example bed file
   >>> bedfn = pybedtools.example_filename('a.bed') 

The full path of *bedfn* will depend on your installation (this is similar
to the ``data()`` function in R_, if you're familiar with that).

Now that you have a filename -- either one of the example files or your
own, you create a new :class:`BedTool` simply by pointing it to that
filename:

.. doctest::

    >>> # create a new BedTool from the example bed file
    >>> myBedTool = BedTool(bedfn)

Alternatively, you can construct BED files from scratch by using the
``from_string`` keyword argument.  However, all spaces will be converted to
tabs using this method, so you'll have to be careful if you add "name"
columns.  This can be useful if you want to create *de novo* BED files on
the fly:

.. doctest::
    :options: +NORMALIZE_WHITESPACE
    
    >>> # an "inline" example:
    >>> fromscratch1 = pybedtools.BedTool('chrX 1 100', from_string=True)
    >>> print(fromscratch1)
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

    >>> fromscratch2 = BedTool(larger_string, from_string=True)
    >>> print(fromscratch2)
    chrX    1   100 feature1    0   +
    chrX    50  350 feature2    0   -
    chr2    5000    10000   another_feature 0   +
    <BLANKLINE>

Of course, you'll usually be using your own bed files that have some
biological importance for your work that are saved in places convenient for
you, for example::

    >>> a = BedTool('/data/sample1/peaks.bed')
