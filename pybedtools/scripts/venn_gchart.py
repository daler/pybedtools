#!/usr/bin/env python
"""
Given 3 files, creates a 3-way Venn diagram of intersections using the Google \
Chart API; see :mod:`pybedtools.contrib.venn_maker` for more flexibility.

The values in the diagram assume:

    * unstranded intersections
    * no features that are nested inside larger features
"""

import argparse
import sys
import pybedtools
import urllib
import urllib2


def venn_gchart(a, b, c=None, colors=None, labels=None, size='300x300'):
    """
    a, b, and c are filenames to BED-like files.

    *colors* is a list of 3 hex colors

    *labels* is a list of 3 labels

    *outfn* is the output PNG you want to create.

    *size* is the size in pixels for the PNG
    """
    a = pybedtools.BedTool(a)
    b = pybedtools.BedTool(b)
    if c:
        c = pybedtools.BedTool(c)

    # The order of values is meaningful to the API, see
    # http://code.google.com/apis/chart/docs/gallery/venn_charts.html
    if c:
        vals = [len(a),
                len(b),
                len(c),
                len(a + b),
                len(a + c),
                len(b + c),
                len(a + b + c)]
    else:
        # insert 0 for size of 3rd circle.
        vals = [len(a), len(b), 0, len(a + b)]
        labels = labels[:2]
    # API doesn't seem to like large numbers, so get fractions instead, then
    # join make a comma-separated list of values.
    mx = float(max(vals))
    vals = [i / mx for i in vals]
    valstr = ','.join(map(str, vals))

    data = {'cht': 'v',
            'chs': size,
            'chd': 't:' + valstr}

    # Add the optional data, if specified
    if labels:
        data['chdl'] = '|'.join(labels)
    if colors:
        data['chco'] = ','.join(colors)
    return data


def gchart(data, outfn='out.png'):
    """
    Sends data to Google Chart API
    """
    data = urllib.urlencode(data)

    url = 'https://chart.googleapis.com/chart?'

    # Request and get the PNG
    req = urllib2.Request(url, data)
    print url + data
    response = urllib2.urlopen(req)
    f = open(outfn, 'w')
    f.write(response.read())
    f.close()


def main():
    """Create a 3-way Venn diagram using Google Charts API
    """

    op = argparse.ArgumentParser(description=__doc__, prog=sys.argv[0])
    op.add_argument('-a', help='File to use for the left-most circle')
    op.add_argument('-b', help='File to use for the right-most circle')
    op.add_argument('-c', help='File to use for the bottom circle')
    op.add_argument('--colors', help='Optional comma-separated list of hex '
                    'colors for circles a, b, and c.  E.g. %(default)s',
                     default='00FF00,FF0000,0000FF')
    op.add_argument('--labels',
            help='Optional comma-separated list of labels for a, b, and c',
                    default='a,b,c')
    op.add_argument('--size', default='300x300',
                  help='Optional size of PNG, in pixels.  Default is '
                  '"%(default)s"')
    op.add_argument('-o', default='out.png',
                  help='Output file to save as, in PNG format')
    op.add_argument('--test', action='store_true',
            help='run test, overriding all other options.')
    options = op.parse_args()

    reqd_args = ['a', 'b']
    if not options.test:
        for ra in reqd_args:
            if not getattr(options, ra):
                op.print_help()
                sys.stderr.write('Missing required arg "%s"\n' % ra)
                sys.exit(1)

    if options.test:
        # Example data
        pybedtools.bedtool.random.seed(1)
        a = pybedtools.example_bedtool('rmsk.hg18.chr21.small.bed')
        b = pybedtools.example_bedtool('venn.b.bed')
        c = pybedtools.example_bedtool('venn.c.bed')
        options.a = a.fn
        options.b = b.fn
        options.c = c.fn
        options.colors = '00FF00,FF0000,0000FF'
        options.o = 'out.png'
        options.labels = 'a,b,c'

    data = venn_gchart(a=options.a, b=options.b, c=options.c,
                       colors=options.colors.split(','),
                       labels=options.labels.split(','),
                       size=options.size)
    gchart(data, outfn=options.o)

if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS).failed == 0:
        main()
