"""
tests for contrib module
"""
import sys
import os
import pybedtools
from pybedtools import Interval
#from pybedtools.contrib import Classifier
from .tfuncs import setup, teardown, testdir, test_tempdir, unwriteable

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


def fix(x):
    """
    Replaces spaces with tabs, removes spurious newlines, and lstrip()s each
    line. Makes it really easy to create BED files on the fly for testing and
    checking.
    """
    s = ""
    for i in  x.splitlines():
        i = i.lstrip()
        if i.endswith('\t'):
            add_tab = '\t'
        else:
            add_tab = ''
        if len(i) == 0:
            continue
        i = i.split()
        i = '\t'.join(i) + add_tab + '\n'
        s += i
    return s

def _classifier():

    c = Classifier(
            bed=pybedtools.example_filename('gdc.bed'),
            annotations=pybedtools.example_filename('gdc.gff'))
    c.classify()

    bed = pybedtools.example_bedtool('gdc.bed')

    assert c.class_counts == {
            frozenset(['UTR', 'exon', 'mRNA', 'CDS', 'tRNA', 'gene']): 1,
            frozenset(['intron', 'gene', 'mRNA']): 3,
            frozenset([]): 1,
            frozenset(['gene', 'exon', 'mRNA', 'CDS']): 2,
            frozenset(['exon', 'mRNA', 'CDS', 'tRNA', 'intron', 'gene']): 1}

    assert c.feature_classes == {
            bed[0]: set(['.']),
            bed[1]: set(['gene', 'exon', 'mRNA', 'CDS']),
            bed[2]: set(['intron', 'gene', 'mRNA']),
            bed[3]: set(['intron', 'gene', 'mRNA']),
            bed[4]: set(['tRNA', 'UTR', 'exon', 'mRNA', 'CDS', 'gene']),
            bed[5]: set(['gene', 'exon', 'mRNA', 'CDS']),
            bed[6]: set(['intron', 'gene', 'mRNA']),
            bed[7]: set(['tRNA', 'intron', 'exon', 'mRNA', 'CDS', 'gene']),
            }

    print('use these indexes for debugging')
    for i, f in enumerate(bed):
        print(i, f)

    for k, v in list(c.class_features.items()):
        print(k)
        for i in v:
            print('\t' + str(i))

    assert c.class_features == {
            frozenset([]): [bed[0]],
            frozenset(['intron', 'gene', 'mRNA']): [bed[6], bed[2], bed[3]],
            frozenset(['gene', 'exon', 'mRNA', 'CDS']): [bed[5], bed[1]],
            frozenset(['UTR', 'exon', 'mRNA', 'CDS', 'tRNA', 'gene']): [bed[4]],
            frozenset(['exon', 'mRNA', 'CDS', 'tRNA', 'intron', 'gene']): [bed[7]],
            }

def test_cleaned_intersect():
    x = pybedtools.BedTool("""
    chr1 1 10      1
    chr1 20 30     2
    chr1 100 120   3
    """, from_string=True)
    y = pybedtools.BedTool("""
    chr1 2 7       4
    chr1 110 120   5
    chr1 200 210   6
    """, from_string=True)
    z = pybedtools.BedTool("""
    chr1 25 40     7
    chr1 190 205   8
    chr1 1000 1001 9
    """, from_string=True)

    # Two-way test
    #
    x2, y2 = pybedtools.contrib.venn_maker.cleaned_intersect([x, y])

    # x should be the same -- 1, 2, 3
    # y should have 1, 3, 6

    assert x2 == fix("""
    chr1 1 10
    chr1 20 30
    chr1 100 120
    """)

    assert y2 == fix("""
    chr1 1 10
    chr1 100 120
    chr1 200 210""")

    # Three-way test
    #
    x3, y3, z3 = pybedtools.contrib.venn_maker.cleaned_intersect([x, y, z])

    # x should be the same -- 1, 2, 3
    # y should have 1, 3, 6
    # z should have 2, 6

    assert x3 == fix("""
    chr1 1 10
    chr1 20 30
    chr1 100 120
    """)

    assert y3 == fix("""
    chr1 1 10
    chr1 100 120
    chr1 200 210""")

    assert z3 == fix("""
    chr1 20 30
    chr1 200 210
    chr1 1000 1001""")

    try:
        pybedtools.helpers._check_for_R()
        print(pybedtools.contrib.venn_maker.venn_maker(
                beds=[x, y, z],
                names=['x','y','z'],
                figure_filename='out.tiff',
                additional_args = ['euler.d=TRUE', 'scaled=TRUE', 'fill=c("red","blue", "orange")'],
                run=True))
    except ValueError:
        sys.stderr.write('R installation not found; skipping test')

    if os.path.exists('out.tiff'):
        os.unlink('out.tiff')


def test_venn_maker_3way_1empty():
    # Fix issue #95. The problem was that BedTool.cat() failed when checking
    # field counts on an empty file.  The fix was to make cat() aware of empty
    # files when checking field num and field type.
    a = pybedtools.BedTool("""
    chr1 10 100
    chr2 10 100""", from_string=True)
    b = pybedtools.BedTool("""
    chr1 12 80""", from_string=True)
    c = pybedtools.BedTool("""
    chr2 20 40""", from_string=True)
    try:
        pybedtools.contrib.venn_maker.venn_maker([a, b, c], run=True, figure_filename='t.tiff')
    except ValueError:
        print("R not installed, skipping test")
        #os.unlink('t.tiff')
