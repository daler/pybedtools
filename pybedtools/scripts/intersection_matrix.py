#!/usr/bin/python
import collections
import sys
import os.path as op
import argparse
from pybedtools import BedTool, example_filename

usage = """

    Send in a list of `N` bed files, and this script will create an N by
    N matrix of their intersections, or optionally, co-localization scores.

    Run the example with::

        %s --test > matrix.txt

    You can then plot a quick heatmap in R with::

        > m = as.matrix(read.table("matrix.txt"))
        > heatmap(m)

""" % sys.argv[0]

ap = argparse.ArgumentParser(usage=usage)
ap.add_argument('beds', nargs="*", help='BED/GTF/GFF/VCF filenames, e.g., '
                'in a directory of bed files, you can use *.bed')
ap.add_argument('--frac', action='store_true',
                help='Instead of counts, report fraction overlapped')
ap.add_argument('--enrichment', action='store_true',
                help='Run randomizations (default 1000, specify otherwise '
                'with --iterations) on each pairwise comparison and compute '
                'the enrichment score as (actual intersection count + 1) / '
                '(median randomized + 1).')
ap.add_argument('--genome', help='Required argument if --enrichment is used.  '
                'Needs to be a string assembly name like "dm3" or "hg19"')
ap.add_argument('--iterations', default=1000, type=int,
                help='Number of randomizations to perform for enrichement '
                'scores')
ap.add_argument('--processes', default=None, type=int,
                help='Number of CPUs to use for randomization')
ap.add_argument('--test', action='store_true', help='Ignore any input BED '
                'files and use test BED files')
args = ap.parse_args()


if not args.beds and not args.test:
    ap.print_help()
    sys.exit(1)

if args.test:
    # insulator binding sites from ChIP-chip -- 4 proteins, 2 cell types
    # Genes Dev. 2009 23(11):1338-1350
    args.beds = [example_filename(i) for i in  [
            'Cp190_Kc_Bushey_2009.bed',
            'Cp190_Mbn2_Bushey_2009.bed',
            'CTCF_Kc_Bushey_2009.bed',
            'CTCF_Mbn2_Bushey_2009.bed',
            'SuHw_Kc_Bushey_2009.bed',
            'SuHw_Mbn2_Bushey_2009.bed',
            'BEAF_Mbn2_Bushey_2009.bed',
            'BEAF_Kc_Bushey_2009.bed'
            ]]


def get_name(fname):
    return op.splitext(op.basename(fname))[0]


def actual_intersection(a, b):
    return len(a.intersect(b, u=True))


def frac_of_a(a, b):
    len_a = float(len(a))
    return len(a.intersect(b, u=True)) / len_a

def enrichment_score(a, b):
    results = a\
            .set_chromsizes(args.genome)\
            .randomstats(b, args.iterations, processes=args.processes)
    return (results['actual'] + 1) / (results['median randomized'] + 1)

if args.enrichment:
    FUNC = enrichment_score
elif args.frac:
    FUNC = frac_of_a
else:
    FUNC = actual_intersection

matrix = collections.defaultdict(dict)
for fa in args.beds:
    a = BedTool(fa)
    for fb in args.beds:
        b = BedTool(fb)
        matrix[get_name(fa)][get_name(fb)] = FUNC(a, b)

keys = sorted(matrix.keys())

print "\t" + "\t".join(keys)
for k in keys:
    print k,
    for j in keys:
        print '\t' + str(matrix[k][j]),
    print
