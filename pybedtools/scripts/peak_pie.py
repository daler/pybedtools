#!/usr/bin/env python
"""
Make a pie chart where peaks fall in annotations; see \
:mod:`pybedtools.contrib.Classifier` for more flexibility.

The results here are similar to CEAS (http://liulab.dfci.harvard.edu/CEAS/).

However, multi-featuretype classes are reported.  That is, if a peak falls in
an exon in one isoform and an intron in another isoform, the class is "exon,
intron".
"""

import sys
import urllib
import urllib2
import argparse
import pybedtools
from collections import defaultdict


def make_pie(bed, gff, stranded=False, out='out.png',
             include=None, exclude=None, thresh=0):

    a = pybedtools.BedTool(bed)
    b = pybedtools.BedTool(gff).remove_invalid()

    c = a.intersect(b,
                    wao=True,
                    s=stranded,
                    stream=True)

    # So we can grab just `a` features later...
    afields = a.field_count()

    # Where we can find the featuretype in the -wao output.  Assumes GFF.
    type_idx = afields + 2

    # 3 different code paths depending on include/exclude to cut down on
    # if/else checks.
    #
    # For un-included featuretypes, put them in the '.' category (unnannotated)
    if include and exclude:
        raise ValueError('Can only specify one of `include` or `exclude`.')
    d = defaultdict(set)
    if include:
        for feature in c:
            featuretype = feature[type_idx]
            key = '\t'.join(feature[:afields])
            if featuretype in include:
                d[key].update([featuretype])
            else:
                d[key].update(['.'])
    elif exclude:
        for feature in c:
            featuretype = feature[type_idx]
            key = '\t'.join(feature[:afields])
            if featuretype not in exclude:
                d[key].update([featuretype])
            else:
                d[key].update(['.'])
    else:
        for feature in c:
            featuretype = feature[type_idx]
            key = '\t'.join(feature[:afields])
            d[key].update([featuretype])

    def labelmaker(x):
        x.difference_update('.')
        label = []
        for i in list(x):
            if i == 'three_prime_UTR':
                i = "3' UTR"
            if i == 'five_prime_UTR':
                i = "5' UTR"
            label.append(i)
        return ', '.join(sorted(label))

    # Prepare results for Google Charts API
    npeaks = float(len(d))
    count_d = defaultdict(int)
    for peak, featuretypes in d.items():
        if featuretypes == set('.'):
            featuretype = 'unannotated'
        else:
            featuretype = labelmaker(featuretypes)
        count_d[featuretype] += 1

    results = count_d.items()
    results.sort(key=lambda x: x[1])
    labels, counts = zip(*results)

    labels = []
    counts_to_use = []
    for label, count in results:
        perc = count / npeaks * 100
        if perc > thresh:
            labels.append('%s: %s (%.1f%%)' % (label,
                                               count,
                                               perc))
            counts_to_use.append(perc)

    # Set up the Gchart data
    data = {'cht': 'p',
            'chs': '750x350',
            'chd': 't:' + ','.join(map(str, counts_to_use)),
            'chl': '|'.join(labels)}

    # Encode it correctly
    encoded_data = urllib.urlencode(data)

    # Send request and get data; write to file
    url = 'https://chart.googleapis.com/chart?'
    req = urllib2.Request(url, encoded_data)
    response = urllib2.urlopen(req)
    f = open(out, 'w')
    f.write(response.read())
    f.close()


def main():
    """
    Make a pie chart of features overlapping annotations (e.g., peaks in
    introns, exons, etc)
    """
    ap = argparse.ArgumentParser(description=__doc__,
                          formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bed', help='BED file of e.g. peaks')
    ap.add_argument('--gff', help='GFF file of e.g. annotations')
    ap.add_argument('--out', default='out.png', help='Output PNG file')
    ap.add_argument('--stranded', action='store_true',
                    help='Use strand-specific intersections')
    ap.add_argument('--include', nargs='*', help='Featuretypes to include')
    ap.add_argument('--exclude', nargs='*', help='Featuretypes to exclude')
    ap.add_argument('--thresh', type=float,
                    help='Threshold percentage below which output will be '
                    'suppressed')
    ap.add_argument('--test', action='store_true',
                    help='Run test, overwriting all other args. Result will '
                    'be "out.png" in current directory.')
    args = ap.parse_args()

    if not (args.bed and args.gff) and not args.test:
        ap.print_help()
        sys.exit(1)

    if not args.test:
        if args.include and args.exclude:
            raise ValueError('Cannot specify both --include and --exclude')

        make_pie(bed=args.bed,
                 gff=args.gff,
                 out=args.out,
                 thresh=args.thresh,
                 stranded=args.stranded,
                 include=args.include,
                 exclude=args.exclude)
    else:
        make_pie(bed=pybedtools.example_filename('gdc.bed'),
                 gff=pybedtools.example_filename('gdc.gff'),
                 stranded=True,
                 out='out.png',
                 include=['exon',
                          'CDS',
                          'intron',
                          'five_prime_UTR',
                          'three_prime_UTR'])


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS).failed == 0:
        main()
