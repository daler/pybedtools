#!/usr/bin/env python

"""
Example from pybedtools documentation: find reads in introns and exons using
multiple CPUs.
"""

from __future__ import print_function
import pybedtools
import argparse
import os
import sys
import multiprocessing

ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),
                             usage=__doc__)
ap.add_argument('--gff', required=True,
                help='GFF or GTF file containing annotations')
ap.add_argument('--bam', required=True,
                help='BAM file containing reads to be counted')
ap.add_argument('--stranded', action='store_true',
                help='Use strand-specific merging and overlap. '
                     'Default is to ignore strand')
ap.add_argument('--processes', default=1, type=int,
                help='Number of processes to use in parallel.')
ap.add_argument('-o', '--output',
                help='Optional file to which results will be written; '
                     'default is stdout')
ap.add_argument('-v', '--verbose', action='store_true',
                help='Verbose (goes to stderr)')
args = ap.parse_args()

gff = args.gff
bam = args.bam
stranded = args.stranded

if args.processes > 3:
    print(
        "Only need 3 processes (one each for exon, intron, both), so "
        "resetting processes from {0} to 3".format(args.processes)
    )
    args.processes = 3


def featuretype_filter(feature, featuretype):
    """
    Only passes features with the specified *featuretype*
    """
    if feature[2] == featuretype:
        return True
    return False


def subset_featuretypes(featuretype):
    """
    Returns the filename containing only `featuretype` features.
    """
    return g.filter(featuretype_filter, featuretype).saveas().fn


def count_reads_in_features(features):
    """
    Callback function to count reads in features
    """
    return (
        pybedtools.BedTool(bam)
        .intersect(
            features,
            s=stranded,
            bed=True,
            stream=True,
        )
    ).count()

# Some GFF files have invalid entries -- like chromosomes with negative
# coords or features of length = 0.  This line removes them and saves the
# result in a tempfile
g = pybedtools.BedTool(gff).remove_invalid().saveas()

# Set up pool of workers
pool = multiprocessing.Pool(processes=args.processes)

# Get separate files for introns and exons in parallel
featuretypes = ['intron', 'exon']
introns, exons = pool.map(subset_featuretypes, featuretypes)
introns = pybedtools.BedTool(introns)
exons = pybedtools.BedTool(exons)

# Identify unique and shared regions
exon_only = exons.subtract(introns).merge().remove_invalid().saveas()
intron_only = introns.subtract(exons).merge().remove_invalid().saveas()
intron_and_exon = (
    exons
    .intersect(introns)
    .merge()
    .remove_invalid()
    .saveas()
)

# Do intersections with BAM file in parallel. Note that we're passing filenames
# to multiprocessing.Pool rather than BedTool objects.
features = (exon_only.fn, intron_only.fn, intron_and_exon.fn)

results = pool.map(count_reads_in_features, features)

labels = ('      exon only:',
          '    intron only:',
          'intron and exon:')

for label, reads in zip(labels, results):
    print('%s %s' % (label, reads))

pybedtools.cleanup(verbose=False)
