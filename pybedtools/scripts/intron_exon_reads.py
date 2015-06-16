#!/usr/bin/env python

"""
Example from pybedtools documentation (:ref:`third example`) to count \
reads in introns and exons using multiple CPUs.
"""

from __future__ import print_function
import pybedtools
import argparse
import os
import sys
import multiprocessing


def featuretype_filter(feature, featuretype):
    """
    Only passes features with the specified *featuretype*
    """
    if feature[2] == featuretype:
        return True
    return False


def subset_featuretypes(featuretype):
    return g.filter(featuretype_filter, featuretype).saveas()


def count_reads_in_features(features):
    """
    Callback function to count reads in features
    """
    return features.intersect(abam=bam,
                             b=features.fn,
                             s=stranded,
                             bed=True,
                             stream=True).count()


def main():
    """
    Third quick example from the documentation -- count reads introns and
    exons, in parallel
    """
    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),
                                 usage=__doc__)
    ap.add_argument('--gff', required=True,
                    help='GFF or GTF file containing annotations')
    ap.add_argument('--bam', required=True,
                    help='BAM file containing reads to be counted')
    ap.add_argument('--stranded', action='store_true',
                    help='Use strand-specific merging and overlap. '
                         'Default is to ignore strand')
    ap.add_argument('--no-parallel', dest='noparallel', action='store_true',
                    help='Disables parallel computation')
    ap.add_argument('-o', '--output',
                    help='Optional file to which results will be written; '
                         'default is stdout')
    ap.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose (goes to stderr)')
    args = ap.parse_args()

    gff = args.gff
    bam = args.bam
    stranded = args.stranded
    parallel = not args.noparallel

    # Some GFF files have invalid entries -- like chromosomes with negative
    # coords or features of length = 0.  This line removes them and saves the
    # result in a tempfile
    g = pybedtools.BedTool(gff).remove_invalid().saveas()

    # Decide which version of map to use.  If parallel, we only need 3
    # processes.
    pool = multiprocessing.Pool(processes=3)

    # Get separate files for introns and exons in parallel (if specified)
    featuretypes = ('intron', 'exon')
    introns, exons = pool.map(subset_featuretypes, featuretypes)

    # Perform some genome algebra to get unique and shared regions
    exon_only = exons.subtract(introns).merge().remove_invalid().saveas()
    intron_only = introns.subtract(exons).merge().remove_invalid().saveas()
    intron_and_exon = exons\
            .intersect(introns).merge().remove_invalid().saveas()

    # Do intersections with BAM file in parallel
    features = (exon_only, intron_only, intron_and_exon)
    results = pool.map(count_reads_in_features, features)

    labels = ('      exon only:',
              '    intron only:',
              'intron and exon:')

    for label, reads in zip(labels, results):
        print('%s %s' % (label, reads))

    pybedtools.cleanup(verbose=False)

if __name__ == "__main__":
    main()
