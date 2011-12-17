"""
tests for contrib module
"""
import pybedtools
from pybedtools.contrib import Classifier

# model for gdc.
# chr2L, starts at 1 and spacing is 10bp.
# import gdc
# g = gdc.GenomeModel(chrom_start=1,
#                     chrom='chr2L',
#                     scalar=10,
#                     read_length=5)
"""
#         10       20        30
#123456789012345678901234567890
   >===||||||~~~~|||==@    #mRNA_xs2_g2_+
               <|||||||@   #tRNA_t2_t2_-
$+     -      --     +     #
$      +      + +          #
"""
#^     ^      ^^^    ^
#|     |      |      |- STRANDED: exon, UTR, mRNA, gene
#      |      |         UNSTRANDED: exon, UTR, mRNA, tRNA, gene, CDS (cause gdc considers tRNA's exon as CDS)
#|     |      |||- STRANDED: intron, mRNA, gene
#|     |      |    UNSTRANDED: intron, mRNA, gene, tRNA, CDS (because GDC considers the tRNA's exon as CDS)
#|     |      ||- STRANDED: none
#|     |      |   UNSTRANDED: intron, mRNA, gene
#|     |      |-STRANDED: +: intron, mRNA, gene; -: none
#|     |        UNSTRANDED: intron, mRNA, gene
#|     |- STRANDED: +: exon, CDS, mRNA, gene; -: none
#|        UNSTRANDED: exon, CDS, mRNA, gene
#- STRANDED: none
#  UNSTRANDED: none



def test_classifier():

    c = Classifier(
            bed=pybedtools.example_filename('gdc.bed'),
            annotations=pybedtools.example_filename('gdc.gff'))
    d = c.classify()

    print "\nclassification dict"
    for k, v in d.items():
        print k, v

    print "\nsummarized"
    for k, v in c.class_counts(d).items():
        print 'class: "%s" count: %s' % (k, v)
    raise NotImplementedError('Still working on the test...')
