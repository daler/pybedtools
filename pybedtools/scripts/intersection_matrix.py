#!/usr/bin/python
import collections
import sys
import os.path as op
from pybedtools import BedTool, example_filename

usage = """
Send in a list of `N` bed files, and this script will create an N by N matrix
of their intersections.

Usage:

    %s *.bed > matrix.txt

To run with test data, use:

    %s --test

""" % (op.basename(sys.argv[0]),
        op.basename(sys.argv[0]))

try:
    bed_files = sys.argv[1:]

    if bed_files[0] == '--test':
        # insulator binding sites from ChIP-chip -- 4 proteins, 2 cell types
        # Genes Dev. 2009 23(11):1338-1350
        bed_files = ['Cp190_Kc_Bushey_2009.bed',
                     'Cp190_Mbn2_Bushey_2009.bed',
                     'CTCF_Kc_Bushey_2009.bed',
                     'CTCF_Mbn2_Bushey_2009.bed',
                     'SuHw_Kc_Bushey_2009.bed',
                     'SuHw_Mbn2_Bushey_2009.bed',
                     'BEAF_Mbn2_Bushey_2009.bed',
                     'BEAF_Kc_Bushey_2009.bed'
                     ]
        bed_files = [example_filename(i) for i in bed_files]

except IndexError:
    sys.stderr.write(usage)
    sys.exit(1)

def get_name(fname):
    return op.splitext(op.basename(fname))[0]

matrix = collections.defaultdict(dict)
for fa in bed_files:
   a = BedTool(fa)
   for fb in bed_files:
       matrix[get_name(fa)][get_name(fb)] = len(a.intersect(fb, u=True))

keys = sorted(matrix.keys())

print "\t".join(keys)
for k in keys:
    print k + '\t',
    for j in keys:
        print '\t' + str(matrix[k][j]),
    print
