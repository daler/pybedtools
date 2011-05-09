#!/usr/env python
"""

    Classifies reads in a BAM file as exon, intron, both (i.e., ambiguous
    because of multiple isoforms) or other (i.e., intergenic).

    Takes a GFF or GTF file as input, and assumes that introns are annotated.

    The GFF file will be filtered to remove possibly problematic features, for
    example, sometimes chromosome features will have negative start coords.

"""
import sys
import os
import time
import argparse
import pybedtools

def classify_reads(gff, bam, stranded=False, output=None, disable_cleanup=False, verbose=False):
    """
    *gff* is a GFF file with intron/exon features annotated.

    *bam* is the reads you'd like to classify.

    *stranded* will only classify a read if it falls into a feature of the same strand

    *output* is an optional output file to write to, otherwise stdout will be used.

    *disable_cleanup* will not automatically clean up the files at the end of script

    *verbose* is useful for debugging (esp. with *disable_cleanup*) 
    """

    a = pybedtools.BedTool(gff)

    # Get a temp filename that will be cleaned up later
    tmp = a._tmp()

    if verbose:
        sys.stderr.write('Cleaning GFF file and only keeping intron/exons...')
        sys.stderr.flush()
        t0 = time.time()

    # Ignore malformed lines, and make a new temp file containing only introns and
    # exons.
    def intron_exon(x):
        if x[2] in ('exon', 'intron'):
            return True
        return False

    def renamer(x):
        x.name = x[2]
        return x

    b = a.remove_invalid().filter(intron_exon).each(renamer)

    if verbose:
        sys.stderr.write('(%.1fs), file=%s\n'% ((time.time()-t0), b.fn))

    if verbose:
        sys.stderr.write('Merging intron/exons...')
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

    if verbose:
        sys.stderr.write('(%.1fs), file=%s\n'% ((time.time()-t0), c.fn))

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

    if verbose:
        sys.stderr.write('(%.1fs), file=%s\n'% ((time.time()-t0), d.fn))

    if verbose:
        sys.stderr.write('Counting reads...')
        sys.stderr.flush()
        t0 = time.time()

    exonic    = 0
    intronic  = 0
    ambiguous = 0
    other     = 0
    total     = 0.0

    # Since strand isn't reported if no -s, we need to adjust where to look for
    # the name 
    if stranded:
        feature_name_ind = -4
    else:
        feature_name_ind = -3

    for feature in d:
        total += 1
        intersected_with = feature[feature_name_ind]
        if intersected_with == 'exon': exonic += 1
        elif intersected_with == 'intron': intronic += 1
        elif 'exon' in intersected_with and 'intron' in intersected_with: ambiguous += 1
        else: other += 1

    if verbose:
        sys.stderr.write('(%.1fs)\n\n' % (time.time()-t0))

    if output:
        fout = open(output, 'w')
    else:
        fout = sys.stdout


    fout.write('   total reads: %s\n' % int(total))
    fout.write('  exonic reads: %s (%.1f%%)\n' % (exonic,    exonic    / total * 100))
    fout.write('intronic reads: %s (%.1f%%)\n' % (intronic,  intronic  / total * 100))
    fout.write('     ambiguous: %s (%.1f%%)\n' % (ambiguous, ambiguous / total * 100))
    fout.write('         other: %s (%.1f%%)\n' % (other,     other     / total * 100))

    if output:
        fout.close()

    if not disable_cleanup:
        pybedtools.cleanup(verbose=verbose)

def main():
    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), description=__doc__)
    ap.add_argument('--gff', help='GFF or GTF file containing annotations')
    ap.add_argument('--bam', help='BAM file containing reads to be counted')
    ap.add_argument('--stranded', action='store_true', help='Use strand-specific merging and overlap.')
    ap.add_argument('--nocleanup', action='store_true', help='Disable automatic deletion of temp files when finished')
    ap.add_argument('-o', help='Optional file to which results will be written; default is stdout')
    ap.add_argument('-v', action='store_true', help='Verbose (goes to stderr)') 
    args = ap.parse_args()

    classify_reads(gff=args.gff,
                   bam=args.bam,
                   stranded=args.stranded,
                   output=args.o,
                   verbose=args.v,
                   disable_cleanup=args.nocleanup)

if __name__ == "__main__":
    main()
