#!/usr/bin/python

"""
Example pybedtools usage; please read script and comments for more info
"""
from pybedtools import BedTool
import pybedtools
import os
import time


def main():
    """
    Quick demo of some pybedtools functionality
    """

    print """

    Quick demo of some pybedtools functionality

    """

    # In order to locate files on disk, this script uses pybedtools'
    # example_filename() method to locate filenames of example data.
    #
    # To save space, example data only consists of the parts needed for this
    # script to run

    # has genes from chr1 and chr21
    gff_fn = pybedtools.example_filename('hg19.gff')

    # dbSNP 128, all of chr21 and first part of chr1
    snps_fn = pybedtools.example_filename('snps.bed.gz')

    # chr21 sequence
    fasta_fn = pybedtools.example_filename('chr21.fa')

    # tiny BAM file
    bam_fn = pybedtools.example_filename('reads.bam')

    # a couple of example exons
    exon_fn = pybedtools.example_filename('exons.gff')

    # subset and sequences will be saved as tempfiles that won't get
    # auto-deleted
    subset_fn = os.path.join(pybedtools.tempfile.gettempdir(),
                          'chr21-genes-with-snps.gff')
    seq_fn = os.path.join(pybedtools.tempfile.gettempdir(),
                          'chr21-genes-with-snps.fa')

    genes = pybedtools.example_bedtool(gff_fn)
    snps = pybedtools.example_bedtool(snps_fn)

    print '''
    Using files:
        genes:%(gff_fn)s
        snps:%(snps_fn)s
        fasta:%(fasta_fn)s
        bam:%(bam_fn)s
        exons:%(exon_fn)s
    ''' % locals()

    print 'Intersecting genes with SNPs...'
    genes_with_snps = genes.intersect(snps, u=True)
    for g in genes_with_snps[:5]:
        print g.chrom, g.start, g.end, len(g)
    print '... (only showing 5)'
    print

    def chrom_filt(g):
        return g.chrom == 'chr21'

    print 'Subsetting on chr21...',
    subset = genes_with_snps.filter(chrom_filt)
    subset = subset.saveas(subset_fn)
    print '(check saved results in "%s")' % subset_fn
    print

    print 'Extracting sequences...',
    subset.sequence(fasta_fn).\
           save_seqs(seq_fn)
    print '(check sequences in "%s")' % seq_fn
    print

    print 'Finding closest genes to intergenic SNPs...'
    intergenic_snps = (snps - genes)
    nearby = genes.closest(intergenic_snps,
                           d=True,
                           stream=True)

    # set a limit; only show 5 nearby genes.
    limit = 4
    for i, gene in enumerate(nearby):
        if i > limit:
            break
        if int(gene[-1]) < 5000:
            print gene.name
    print '... (only showing 5)'
    print

    t0 = time.time()

    reads = pybedtools.BedTool(bam_fn)
    exons = pybedtools.BedTool(exon_fn)
    on_target = reads.intersect(exons)

    print 'Coverage bedGraph of reads...'
    bedgraph = reads.genome_coverage(genome='hg19',
                                     bg=True)
    for feature in bedgraph[:5]:
        print feature
    print '... (only showing 5)'
    print

if __name__ == "__main__":
    main()
