#!/usr/bin/python
"""
Example from the manuscript to print the names of genes that are <5000 bp away
from intergenic SNPs.  See sh_ms_example.sh for the shell script equivalent.
"""
from pybedtools import BedTool

snps = BedTool('../test/data/snps.bed.gz')
genes = BedTool('../test/data/hg19.gff')

intergenic_snps = (snps - genes)

nearby = genes.closest(intergenic_snps, d=True, stream=True)

for gene in nearby:
    if int(gene[-1]) < 5000:
        print gene.name
