#!/usr/bin/python
"""
Example from the manuscript; see sh_ms_example.sh for the shell script \
equivalent.

Prints the names of genes that are <5000 bp away from intergenic SNPs.
"""
from os import path
from pybedtools import BedTool


def main():
    """
    Runs Python example from the manuscript
    """
    bedtools_dir = path.split(__file__)[0]
    snps = BedTool(path.join(bedtools_dir, '../test/data/snps.bed.gz'))
    genes = BedTool(path.join(bedtools_dir, '../test/data/hg19.gff'))

    intergenic_snps = (snps - genes)

    nearby = genes.closest(intergenic_snps, d=True, stream=True)

    for gene in nearby:
        if int(gene[-1]) < 5000:
            print(gene.name)

if __name__ == "__main__":
    main()
