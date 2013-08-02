#!/usr/bin/env python
"""
Given 3 files, creates a 3-way Venn diagram of intersections using matplotlib; \
see :mod:`pybedtools.contrib.venn_maker` for more flexibility.

Numbers are placed on the diagram.  If you don't have matplotlib installed.
try venn_gchart.py to use the Google Chart API instead.

The values in the diagram assume:

    * unstranded intersections
    * no features that are nested inside larger features
"""

import argparse
import sys
import os
import pybedtools

def venn_mpl(a, b, c, colors=None, outfn='out.png', labels=None):
    """
    *a*, *b*, and *c* are filenames to BED-like files.

    *colors* is a list of matplotlib colors for the Venn diagram circles.

    *outfn* is the resulting output file.  This is passed directly to
    fig.savefig(), so you can supply extensions of .png, .pdf, or whatever your
    matplotlib installation supports.

    *labels* is a list of labels to use for each of the files; by default the
    labels are ['a','b','c']
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle
    except ImportError:
        sys.stderr.write('matplotlib is required to make a Venn diagram with %s\n' % os.path.basename(sys.argv[0]))
        sys.exit(1)

    a = pybedtools.BedTool(a)
    b = pybedtools.BedTool(b)
    c = pybedtools.BedTool(c)

    if colors is None:
        colors = ['r','b','g']

    radius = 6.0
    center = 0.0
    offset = radius / 2

    if labels is None:
        labels = ['a','b','c']

    circle_a = Circle(xy = (center-offset, center+offset), radius=radius, edgecolor=colors[0], label=labels[0])
    circle_b = Circle(xy = (center+offset, center+offset), radius=radius, edgecolor=colors[1], label=labels[1])
    circle_c = Circle(xy = (center,        center-offset), radius=radius, edgecolor=colors[2], label=labels[2])


    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111)

    for circle in (circle_a, circle_b, circle_c):
        circle.set_facecolor('none')
        circle.set_linewidth(3)
        ax.add_patch(circle)

    ax.axis('tight')
    ax.axis('equal')
    ax.set_axis_off()


    kwargs = dict(horizontalalignment='center')

    # Unique to A
    ax.text( center-2*offset, center+offset, str((a - b - c).count()), **kwargs)

    # Unique to B
    ax.text( center+2*offset, center+offset, str((b - a - c).count()), **kwargs)

    # Unique to C
    ax.text( center, center-2*offset, str((c - a - b).count()), **kwargs)

    # A and B not C
    ax.text( center, center+2*offset-0.5*offset, str((a + b - c).count()), **kwargs)

    # A and C not B
    ax.text( center-1.2*offset, center-0.5*offset, str((a + c - b).count()), **kwargs)

    # B and C not A
    ax.text( center+1.2*offset, center-0.5*offset, str((b + c - a).count()), **kwargs)

    # all
    ax.text( center, center, str((a + b + c).count()), **kwargs)

    ax.legend(loc='best')

    fig.savefig(outfn)

    plt.close(fig)

def main():
    """Create a 3-way Venn diagram, using matplotlib"""
    op = argparse.ArgumentParser(description=__doc__, prog=sys.argv[0])
    op.add_argument('-a', help='File to use for the left-most circle')
    op.add_argument('-b', help='File to use for the right-most circle')
    op.add_argument('-c', help='File to use for the bottom circle')
    op.add_argument('--labels',
                  help='Optional comma-separated list of '
                       'labels for a, b, and c', default='a,b,c')
    op.add_argument('--colors', default='r,b,g',
                  help='Comma-separated list of matplotlib-valid colors '
                       'for circles a, b, and c.  E.g., --colors=r,b,k')
    op.add_argument('-o', default='out.png', 
                  help='Output file to save as.  Extension is '
                       'meaningful, e.g., out.pdf, out.png, out.svg.  Default is "%(default)s"')
    op.add_argument('--test', action='store_true', help='run test, overriding all other options.')
    options = op.parse_args()

    reqd_args = ['a','b','c']
    if not options.test:
        for ra in reqd_args:
            if not getattr(options,ra):
                op.print_help()
                sys.stderr.write('Missing required arg "%s"\n' % ra)
                sys.exit(1)

    if options.test:
        pybedtools.bedtool.random.seed(1)
        a = pybedtools.example_bedtool('rmsk.hg18.chr21.small.bed')
        b = pybedtools.example_bedtool('venn.b.bed')
        c = pybedtools.example_bedtool('venn.c.bed')
        options.a = a.fn
        options.b = b.fn
        options.c = c.fn
        options.colors='r,b,g'
        options.o = 'out.png'
        options.labels = 'a,b,c'

    venn_mpl(a=options.a, b=options.b, c=options.c, 
             colors=options.colors.split(','),
             labels=options.labels.split(','), 
             outfn=options.o)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS).failed == 0:
        main()
