"""
Interface between pybedtools and the R package VennDiagram.

Rather than depend on the user to have rpy2 installed, this simply writes an
R script that can be edited and tweaked by the user before being run in R.
"""
from __future__ import print_function
import os
import string
import pybedtools
import six
from pybedtools import helpers
import subprocess
from collections import OrderedDict

# really just fill in x and filename...leave the rest up to the user.
#
# Note that the closing parentheses is missing -- that's so the user can add
# kwargs from the calling function
template = string.Template("""
library(VennDiagram)
venn.diagram(
    x=$x,
    filename=$filename,
    category.names = $names
""")


def _list_to_R_syntax(x):
    """
    Convert items in `x` to a string, and replace tabs with pipes in Interval
    string representations.  Put everything into an R vector and return as one
    big string.
    """
    items = []
    for i in x:
        if isinstance(i, pybedtools.Interval):
            i = str(i).replace('\t', '|')
        items.append('"%s"' % i)
    return 'c(%s)' % ','.join(items)


def _dict_to_R_named_list(d):
    """
    Calls _list_to_R_syntax for each item.  Returns one big string.
    """
    items = []
    for key, val in list(d.items()):
        items.append('"%s" = %s' % (key, _list_to_R_syntax(val)))
    return 'list(%s)' % ', '.join(items)


def truncator(feature):
    """
    Convert a feature of any format into a BED3 format.
    """
    return pybedtools.create_interval_from_list(
            [feature.chrom, str(feature.start), str(feature.stop)])


def cleaned_intersect(items):
    """
    Perform interval intersections such that the end products have identical \
    features for overlapping intervals.

    The VennDiagram package does *set* intersection, not *interval*
    intersection.  So the goal here is to represent intersecting intervals as
    intersecting sets of strings.

    Doing a simple BEDTools intersectBed call doesn't do the trick (even with
    the -u argument).  As a concrete example, what would the string be for an
    intersection of the feature "chr1:1-100" in file `x` and "chr1:50-200" in
    file `y`?

    The method used here is to substitute the intervals in `y` that overlap `x`
    with the corresponding elements in `x`.  This means that in the resulting
    sets, the overlapping features are identical.  To follow up with the
    example, both `x` and `y` would have an item "chr1:50-200" in their sets,
    simply indicating *that* one interval overlapped.

    Venn diagrams are not well suited for nested overlaps or multi-overlaps.
    To illustrate, try drawing the 2-way Venn diagram of the following two
    files. Specifically, what number goes in the middle -- the number of
    features in `x` that intersect `y` (1) or the number of features in `y`
    that intersect `x` (2)?::

        x:
            chr1  1  100
            chr1 500 6000

        y:
            chr1 50 100
            chr1 80 200
            chr9 777 888

    In this case, this function will return the following sets::

        x:
            chr1:1-100
            chr1:500-6000

        y:
            chr1:1-100
            chr9:777-888

    This means that while `x` does not change in length, `y` can.  For example,
    if there are 2 features in `x` that overlap one feature in `y`, then `y`
    will gain those two features in place of its single original feature.

    This strategy is extended for multiple intersections -- see the source for
    details.
    """
    if len(items) == 2:
        x = items[0].each(truncator).saveas()
        y = items[1].each(truncator).saveas()

        # Combine the unique-to-y intervals with the shared-with-x intervals.
        # Since x is first in x+y, resulting features are from x.
        new_y = (y - x).cat(x + y)
        return x, new_y

    if len(items) == 3:
        x = items[0].each(truncator).saveas()
        y = items[1].each(truncator).saveas()
        z = items[2].each(truncator).saveas()

        # Same as above.  Don't care about z yet; this means that y will not
        # change because of z.
        new_y = (y - x).cat(x + y)

        # Combine:
        #  unique-to-z
        #  shared-with-any-x
        #  shared-with-unique-to-y
        new_z = (z - y - x).cat(x + z).cat((y - x) + z)
        return x, new_y, new_z

    if len(items) == 4:
        x = items[0].each(truncator).saveas()
        y = items[1].each(truncator).saveas()
        z = items[2].each(truncator).saveas()
        q = items[3].each(truncator).saveas()

        # Same as 2-way
        new_y = (y - x).cat(x + y)

        # Same as 3-way
        new_z = (z - y - x).cat(x + z).cat((y - x) + z)

        # Combine:
        #  unique-to-q
        #  shared-with-any-x
        #  shared-with-unique-to-y
        #  shared-with-unique-to-z
        new_q = (q - z - y - x)\
                .cat(x + q)\
                .cat((y - x) + q)\
                .cat((z - y - x) + q)

        return x, new_y, new_z, new_q


def venn_maker(beds, names=None, figure_filename=None, script_filename=None,
        additional_args=None, run=False):
    """
    Given a list of interval files, write an R script to create a Venn \
    diagram of overlaps (and optionally run it).

    The R script calls the venn.diagram function of the R package VennDiagram
    for extremely flexible Venn and Euler diagram creation.  Uses
    `cleaned_intersect()` to create string representations of shared intervals.

    `beds` is a list of up to 4 filenames or BedTools.

    `names` is a list of names to use for the Venn diagram, in the same order
    as `beds`. Default is "abcd"[:len(beds)].

    `figure_filename` is the TIFF file to save the figure as.

    `script_filename` is the optional filename to write the R script to

    `additional_args` is list that will be inserted into the R script,
    verbatim.  For example, to use scaled Euler diagrams with different colors,
    use::

        additional_args = ['euler.d=TRUE',
                           'scaled=TRUE',
                           'cat.col=c("red","blue")']

    If `run` is True, then assume R is installed, is on the path, and has
    VennDiagram installed . . . and run the script.  The resulting filename
    will be saved as `figure_filename`.
    """

    if figure_filename is None:
        figure_filename = 'NULL'
    else:
        figure_filename = '"%s"' % figure_filename

    if names is None:
        names = "abcd"[:len(beds)]

    _beds = []
    for bed in beds:
        if not isinstance(bed, pybedtools.BedTool):
            bed = pybedtools.BedTool(bed)
        _beds.append(bed)

    cleaned = cleaned_intersect(_beds)
    results = OrderedDict(list(zip(names, cleaned)))

    s = template.substitute(
            x=_dict_to_R_named_list(results),
            filename=figure_filename,
            names=_list_to_R_syntax(names))
    if additional_args:
        s += ',' + ', '.join(additional_args)

    s += ")"

    if not script_filename:
        fn = pybedtools.BedTool._tmp()
    else:
        fn = script_filename

    fout = open(fn, 'w')
    fout.write(s)
    fout.close()

    out = fn + '.Rout'
    if run:

        if not pybedtools.settings._R_installed:
            helpers._check_for_R()

        cmds = [os.path.join(pybedtools.settings._R_path, 'R'), 'CMD', 'BATCH',
                fn, out]
        p = subprocess.Popen(
                cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        stdout, stderr = p.communicate()
        if stdout or stderr:
            print("stdout:", stdout)
            print("stderr:", stderr)

    if not script_filename:
        return s

    return None
