#!/usr/env python
"""

    Classifies reads in a BAM file by feature type as annotated in a GFF or GTF
    file.

    Reads that fall in overlapping features are placed into a class named by a
    semicolon-separated list of the features.

    You can specify a list of feature types to include (--include) or to ignore
    (--exclude); otherwise all feature types in the annotation file will be
    considered.

    The GFF file will be filtered to remove possibly problematic features; for
    example, sometimes chromosome features will have negative start coords.

    Example usage:

        classify_reads.py --gff dm3.gff --bam sample.bam --include intron exon

    Example results:

        class         count   percent of total
        total         47939   100.00
        exon          38596   80.51
        exon;intron   7072    14.75
        intron        1422    2.97
        unannotated   849     1.77

"""
import sys
import os
import time
import argparse
from collections import defaultdict
import pybedtools


def classify_reads(gff, bam, stranded=False, include=None, exclude=None,
                   output=None, pool=False, disable_cleanup=False,
                   verbose=False):
    """
    *gff* is a GFF file

    *bam* is the reads you'd like to classify

    If *stranded*, will report separate counts for sense and antisense

    *include* is a list of the feature types to include.  By default, include
    everything in the GFF

    *exclude* is a list of featuretypes to exclude.

    *output* is an optional output file to write to, otherwise stdout will be
    used.

    *disable_cleanup* will not automatically clean up the files at the end of
    script

    *verbose* is useful for debugging (esp. with *disable_cleanup*)
    """
    if include and exclude:
        raise ValueError('Both include and exclude specified')

    a = pybedtools.BedTool(gff)

    # Get a temp filename that will be cleaned up later
    tmp = a._tmp()

    if verbose:
        sys.stderr.write('Cleaning GFF file and keeping features '
                         'specified by include/exclude...')
        sys.stderr.flush()
        t0 = time.time()

    # Decide which filtering func to use.

    if include:
        def _filter(x):
            if x[2] in include:
                return True
            return False

    elif exclude:
        def _filter(x):
            if x[2] in exclude:
                return False
            return True
    else:
        def _filter(x):
            return x

    def strand_reverser(x):
        if x.strand == '+':
            x.strand = '-'
            return x
        if x.strand == '-':
            x.strand = '+'
            return x
        else:
            return x

    def renamer(x):
        x.name = x[2]
        return x

    tmp = a._tmp()
    b = a.remove_invalid().filter(_filter).each(renamer).saveas(tmp)

    if stranded:
        b_rev = b.each(strand_reverser)

    if verbose:
        sys.stderr.write('(%.1fs), file = %s\n' % ((time.time() - t0), b.fn))

    if verbose:
        sys.stderr.write('Merging features...')
        sys.stderr.flush()
        t0 = time.time()

    # Explanation for args to merge():
    #
    # nms=True
    #       creates new features like exon;intron
    # d=-1
    #       means that bookended features (which introns and exon are) will not
    #       be merged.
    # scores='sum'
    #       makes results a valid BED file (chr-start-stop-name-score-strand
    #       [OK] , vs chr-start-stop-name-strand [not valid])
    c = b.merge(nms=True, d=-1, s=stranded, scores='sum')

    sequence_space = defaultdict(int)
    for feature in c:
        nm = ';'.join(sorted(list(set(feature.name.split(';')))))
        sequence_space[nm] += len(feature)
        sequence_space['total'] += len(feature)

    if stranded:
        c_rev = b_rev.merge(nms=True, d=-1, s=stranded, scores='sum')

    if verbose:
        sys.stderr.write('(%.1fs), file=%s\n' % ((time.time() - t0), c.fn))

    if verbose:
        sys.stderr.write('Intersecting with BAM file...')
        sys.stderr.flush()
        t0 = time.time()

    # Here we're using c's filename for *b*.  Create BED output, and make sure
    # all reads are written to file.
    #
    # (note: since both abam and b are provided, it actually doesn't matter
    # which BedTool we use for this)
    d = c.intersect(abam=bam, b=c.fn, bed=True, wao=True, s=stranded)

    if stranded:
        d_rev = c_rev.intersect(abam=bam,
                                b=c_rev.fn,
                                bed=True,
                                wao=True,
                                s=stranded)

    if verbose:
        sys.stderr.write('(%.1fs), file=%s\n' % ((time.time() - t0), d.fn))

    if verbose:
        sys.stderr.write('Counting reads...')
        sys.stderr.flush()
        t0 = time.time()

    total = 0.0
    results = defaultdict(int)

    # Since strand isn't reported if no -s, we need to adjust where to look for
    # the name
    if stranded:
        feature_name_ind = -4
    else:
        feature_name_ind = -3

    for feature in d:
        total += 1

        # This is the key piece of information, and will be a
        # semicolon-separated list of intersected featuretypes
        intersected_with = feature[feature_name_ind]

        # 'pool' will count one read multiple times -- it will increment each
        # category it is found within.  Total category counts will not add to
        # the total in this case.
        if pool:
            for key in intersected_with.split(';'):
                results[key] += 1

        # Otherwise, unique and sort the components of the class, and only
        # count the read for this particular class
        else:
            key = ';'.join(sorted(list(set(intersected_with.split(';')))))
            results[key] += 1

    if stranded:
        results_rev = defaultdict(int)
        total_rev = 0.0
        for feature in d_rev:
            total_rev += 1
            intersected_with = feature[feature_name_ind]
            key = ';'.join(sorted(list(set(intersected_with.split(';')))))
            results_rev[key] += 1

    if verbose:
        sys.stderr.write('(%.1fs)\n\n' % (time.time() - t0))

    if output:
        fout = open(output, 'w')
    else:
        fout = sys.stdout

    # Format the dictionaries for output.
    results['total'] = int(total)
    if stranded:
        results_rev['total'] = int(total_rev)

    # re-name the unannotated
    try:
        results['unannotated'] = results.pop('.')
    except KeyError:
        pass
    if stranded:
        try:
            results_rev['unannotated'] = results_rev.pop('.')
        except KeyError:
            pass

    # desired output is
    # class <TAB> count <TAB> percent of total

    sorted_results = sorted(results.items(), key=lambda x: x[1], reverse=True)

    if not stranded:
        header = ['class', 'count', 'percent of total', 'sequence space']
        fout.write('\t'.join(header) + '\n')

        for key, value in sorted_results:
            values = (key,
                      value,
                      value / total * 100,
                      sequence_space[key])

            fout.write('%s\t%s\t%.2f\t%s\n' % values)

    else:
        header = ['class',
                  'count (sense)',
                  'percent of total (sense)',
                  'count (antisense)',
                  'percent of total (antisense)',
                  'sequence space',
                 ]
        fout.write('\t'.join(header) + '\n')
        for key, value in sorted_results:
            try:
                value_rev = results_rev.pop(key)
            except KeyError:
                value_rev = 0

            values = (key,
                      value,
                      value / total * 100,
                      value_rev,
                      value_rev / total_rev * 100,
                      sequence_space[key])

            fout.write('%s\t%s\t%.2f\t%s\t%.1f\t%s\n' % values)

        sorted_results_rev = sorted(results_rev.items(),
                                    key=lambda x: x[1],
                                    reverse=True)

        for key, value_rev in sorted_results_rev:
            try:
                value = results[key]
            except KeyError:
                value = 0

            values = (key,
                      value,
                      value / total * 100,
                      value_rev,
                      value_rev / total_rev * 100,
                      sequence_space[key])

            fout.write('%s\t%s\t%.2f\t%s\t%.1f\t%s\n' % values)

    if output:
        fout.close()

    if not disable_cleanup:
        pybedtools.cleanup(verbose=verbose)


def main():
    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),
                                 usage=__doc__)
    ap.add_argument('--gff', required=True,
                    help='GFF or GTF file containing annotations')
    ap.add_argument('--bam', required=True,
                    help='BAM file containing reads to be counted')
    ap.add_argument('--stranded', action='store_true',
                    help='Use strand-specific merging and overlap. '
                         'Default is to ignore strand')
    ap.add_argument('--nocleanup', action='store_true',
                    help='Disable automatic deletion of temp files '
                         'when finished. Useful for debugging.')
    ap.add_argument('-o', '--output',
                    help='Optional file to which results will be written; '
                         'default is stdout')
    ap.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose (goes to stderr)')
    ap.add_argument('--include', nargs='*', help='Feature types to include.')
    ap.add_argument('--exclude', nargs='*', help='Feature types to exclude.')
    ap.add_argument('--pool', action='store_true',
                    help='For a read that falls in multiple featuretypes '
                         '(like exon;intron), then increment each featuretype '
                         '(exons += 1, introns +=1) instead of the default '
                         'behavior of only incrementing the compound class '
                         '(exon;intron += 1)')

    args = ap.parse_args()

    classify_reads(gff=args.gff,
                   bam=args.bam,
                   stranded=args.stranded,
                   include=args.include,
                   exclude=args.exclude,
                   pool=args.pool,
                   output=args.output,
                   verbose=args.verbose,
                   disable_cleanup=args.nocleanup)

if __name__ == "__main__":
    main()
