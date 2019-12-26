#!/usr/bin/env python

"""
Example from pybedtools documentation: find reads in introns and exons using
multiple CPUs.

Prints a tab-separated file containing class (exon, intron, both) and number of
reads in each class.
"""
from __future__ import print_function
import pybedtools
import argparse
import os
import sys
import multiprocessing

if __name__ == "__main__":

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), usage=__doc__)
    ap.add_argument(
        "--gff", required=True, help="GFF or GTF file containing annotations"
    )
    ap.add_argument(
        "--bam", required=True, help="BAM file containing reads to be counted"
    )
    ap.add_argument(
        "--stranded",
        action="store_true",
        help="Use strand-specific merging and overlap. " "Default is to ignore strand",
    )
    ap.add_argument(
        "--processes",
        default=1,
        type=int,
        help="Number of processes to use in parallel.",
    )
    ap.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose (goes to stderr)"
    )
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
            pybedtools.BedTool(bam).intersect(
                features, s=stranded, bed=True, stream=True
            )
        ).count()

    # Some GFF files have invalid entries -- like chromosomes with negative coords
    # or features of length = 0.  This line removes them (the `remove_invalid`
    # method) and saves the result in a tempfile
    g = pybedtools.BedTool(gff).remove_invalid().saveas()

    # Set up pool of workers
    pool = multiprocessing.Pool(processes=args.processes)

    # Get separate files for introns and exons in parallel
    featuretypes = ["intron", "exon"]
    introns, exons = pool.map(subset_featuretypes, featuretypes)

    # Since `subset_featuretypes` returns filenames, we convert to BedTool objects
    # to do intersections below.
    introns = pybedtools.BedTool(introns)
    exons = pybedtools.BedTool(exons)

    # Identify unique and shared regions using bedtools commands subtract, merge,
    # and intersect.
    exon_only = exons.subtract(introns).merge()
    intron_only = introns.subtract(exons).merge()
    intron_and_exon = exons.intersect(introns).merge()

    # Do intersections with BAM file in parallel. Note that we're passing filenames
    # to multiprocessing.Pool rather than BedTool objects.
    features = (exon_only.fn, intron_only.fn, intron_and_exon.fn)

    # Run count_reads_in_features in parallel over features
    results = pool.map(count_reads_in_features, features)

    labels = ("exon_only", "intron_only", "intron_and_exon")

    for label, reads in zip(labels, results):
        print("{0}\t{1}".format(label, reads))
