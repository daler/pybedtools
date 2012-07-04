#!/usr/bin/python
"""
Create a matrix of many pairwise intersections; see \
:mod:`pybedtools.contrib.IntersectionMatrix` for more flexibility
"""

import collections
import time
import sys
import os.path as op
import argparse
import pybedtools
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


def get_name(fname):
    return op.splitext(op.basename(fname))[0]


def actual_intersection(a, b):
    return len(a.intersect(b, u=True))


def frac_of_a(a, b):
    len_a = float(len(a))
    return len(a.intersect(b, u=True)) / len_a


def enrichment_score(a, b, genome_fn, iterations=None, processes=None):
    results = a\
            .randomstats(b, new=True, genome_fn=genome_fn, iterations=iterations, processes=processes)
    return (results['actual'] + 1) / (results['median randomized'] + 1)


def create_matrix(beds, func, verbose=False, **kwargs):
    nfiles = len(beds)
    total = nfiles ** 2
    i = 0
    matrix = collections.defaultdict(dict)
    for fa in beds:
        a = BedTool(fa)
        for fb in beds:
            i += 1
            b = BedTool(fb)

            if verbose:
                sys.stderr.write(
                        '%(i)s of %(total)s: %(fa)s + %(fb)s\n' % locals())
                sys.stderr.flush()

            matrix[get_name(fa)][get_name(fb)] = func(a, b, **kwargs)

    return matrix


def main():
    """
    Creates a pairwise matrix containing overlapping feature counts for many
    BED files
    """
    ap = argparse.ArgumentParser(usage=usage)
    ap.add_argument('beds', nargs="*", help='BED/GTF/GFF/VCF filenames, e.g., '
                    'in a directory of bed files, you can use *.bed')
    ap.add_argument('--frac', action='store_true',
                    help='Instead of counts, report fraction overlapped')
    ap.add_argument('--enrichment', action='store_true',
                    help='Run randomizations (default 1000, specify otherwise '
                    'with --iterations) on each pairwise comparison and '
                    'compute the enrichment score as '
                    '(actual intersection count + 1) / (median randomized + 1)'
                    )
    ap.add_argument('--genome', help='Required argument if --enrichment is '
                    'used. Needs to be a string assembly name like "dm3" or '
                    '"hg19"')
    ap.add_argument('--iterations', default=1000, type=int,
                    help='Number of randomizations to perform for enrichement '
                    'scores')
    ap.add_argument('--processes', default=None, type=int,
                    help='Number of CPUs to use for randomization')
    ap.add_argument('--test', action='store_true', help='Ignore any input BED '
                    'files and use test BED files')
    ap.add_argument('-v', '--verbose', action='store_true',
                    help='Be verbose: print which files are '
                    'currently being intersected and timing info at the end.')
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

    if args.enrichment:
        FUNC = enrichment_score
        genome_fn = pybedtools.chromsizes_to_file(pybedtools.chromsizes(args.genome))
        kwargs = dict(genome_fn=genome_fn, iterations=args.iterations,
                processes=args.processes)

    elif args.frac:
        FUNC = frac_of_a
        kwargs = {}
    else:
        FUNC = actual_intersection
        kwargs = {}

    t0 = time.time()
    matrix = create_matrix(beds=args.beds, func=FUNC, verbose=args.verbose, **kwargs)
    t1 = time.time()

    nfiles = len(args.beds)

    if args.verbose:
        sys.stderr.write('Time to construct %s x %s matrix: %.1fs' \
                % (nfiles, nfiles, (t1 - t0)) + '\n')
    keys = sorted(matrix.keys())

    sys.stdout.write("\t" + "\t".join(keys) + '\n')
    for k in keys:
        sys.stdout.write(k)
        for j in keys:
            sys.stdout.write('\t' + str(matrix[k][j]))
        sys.stdout.write('\n')

if __name__ == "__main__":
    main()
