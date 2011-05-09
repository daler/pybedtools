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

def classify_reads(gff, bam, stranded=False, output=None, verbose=False):

    a = pybedtools.BedTool(gff)

    if verbose:
        sys.stderr.write('Cleaning GFF file and only keeping intron/exons...')
        sys.stderr.flush()
        t0 = time.time()

    # Get a temp filename that will be cleaned up later
    tmp = a._tmp()

    # Ignore malformed lines, and make a new temp file containing only introns and
    # exons.
    fout = open(tmp,'w')
    i = iter(a)
    while True:
        try:
            feature = i.next()
            if feature[2] in ('exon', 'intron'):
                feature.name = feature[2]
                fout.write(str(feature)+'\n')
        except pybedtools.MalformedBedLineError:
            continue
        except StopIteration:
            break
    fout.close()

    if verbose:
        sys.stderr.write('(%.1fs)\n'% (time.time()-t0))

    if verbose:
        sys.stderr.write('Merging intron/exons...')
        sys.stderr.flush()
        t0 = time.time()

    # nms=True creates new features like exon;intron
    # d=-1 means that bookended features (which introns and exon are) will not be
    # merged.
    b = pybedtools.BedTool(tmp).merge(nms=True, d=-1, s=stranded)

    if verbose:
        sys.stderr.write('(%.1fs)\n'% (time.time()-t0))

    if verbose:
        sys.stderr.write('Intersecting with BAM file...')
        sys.stderr.flush()
        t0 = time.time()

    # Here we're using b's filename for *b*.  Create BED output, and make sure all
    # reads are written to file.
    c = b.intersect(abam=bam, b=b.fn, bed=True, wao=True, s=stranded)

    if verbose:
        sys.stderr.write('(%.1fs)\n'% (time.time()-t0))

    if verbose:
        sys.stderr.write('Counting reads...')
        sys.stderr.flush()
        t0 = time.time()

    exonic    = 0
    intronic  = 0
    ambiguous = 0
    other     = 0
    total     = 0.0

    for feature in c:
        total += 1
        intersected_with = feature[-2]
        if intersected_with == 'exon': exonic += 1
        elif intersected_with == 'intron': intronic += 1
        elif intersected_with == 'exon;intron': ambiguous += 1
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

    pybedtools.cleanup(verbose=False)

def main():
    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), description=__doc__)
    ap.add_argument('--gff', help='GFF or GTF file containing annotations')
    ap.add_argument('--bam', help='BAM file containing reads to be counted')
    ap.add_argument('--stranded', help='Use strand-specific merging and overlap.')
    ap.add_argument('-o', help='Optional file to which results will be written; default is stdout')
    ap.add_argument('-v', action='store_true', help='Verbose (goes to stderr)') 
    args = ap.parse_args()

    classify_reads(args.gff, args.bam, args.o, args.v, args.stranded)

if __name__ == "__main__":
    main()
