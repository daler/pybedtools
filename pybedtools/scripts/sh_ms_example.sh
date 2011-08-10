#!/bin/bash

# Shell script to do the same analysis as py_ms_example.py -- namely,
# print the names of genes that are <5000 bp away from intergenic SNPs
#
# See below for timing comparisons.

snps=../test/data/snps.bed.gz
genes=../test/data/hg19.gff
intergenic_snps=/tmp/intergenic_snps

# see note #1 below
snp_fields=`zcat $snps | awk '(NR == 2){print NF; exit;}'`
gene_fields=9
distance_field=$(($gene_fields + $snp_fields + 1))

intersectBed -a $snps -b $genes -v > $intergenic_snps

# see note #2 below
closestBed -a $genes -b $intergenic_snps -d \
| awk -v dist_field="$distance_field" '{if ($dist_field < 5000) print $9;}' \
| perl -ne 'm/[ID|Name|gene_id]=(.*?);/; print "$1\n"'

# see note #3 below
rm $intergenic_snps







# -----corresponding pybedtools script (see py_ms_example.py in this same dir)----
#
#  from pybedtools import BedTool
#  snps = BedTool('../test/data/snps.bed.gz')
#  genes = BedTool('../test/data/hg19.gff')
#
#  intergenic_snps = (snps - genes)
#  nearby = genes.closest(intergenic_snps,
#                         d=True,
#                         stream=True)
#  for gene in nearby:
#      if int(gene[-1]) < 5000:
#         print gene.name

#------------------------------------------------------------------------------
# Note 1:
#   pybedtools allows the user to index into the fields of a feature using
#   Python indexing.  The "gene[-1]" in the python example above gets last item
#   in the line, no matter what the line sizes of the input files.
#
#   In bash, we need to explicitly check how many fields the original files
#   have so we can know the total field count after the the closestBed
#   operation.  This will change depending on the input BED file, since BED
#   files can have 3-6, 9, or 12 columns.  GFF is defined to have 9 fields.
#   Pre-calculating these numbers saves a little time by not having to get the
#   number of fields for each line in the awk script
#------------------------------------------------------------------------------
# Note 2:
#   First, we use the distance field calculated previously.  The closestBed
#   operation was driven by a GFF file, so that file's fields will come first.
#   GFF files are defined to have 9 fields, and the 9th contains the
#   attributes . . . so we can rely on the 9th field being GFF attributes.
#
#   pybedtools stores GFF attributes as a dictionary, and can pull out the name
#   of a gene with the `name` attribute, as in the Python line "print
#   gene.name".  There is no standard for GFF formats; the gene name could be
#   listed as "ID", "Name", or "gene_id" attributes. Here, the Perl regular
#   expression performs the same function
#------------------------------------------------------------------------------
# Note 3:
#   pybedtools automatically deletes temporary files from disk. The Python line
#   "intergenic_snps = (snps - genes)" saves a temporary file for the length of
#   the script, and that file is deleted at exit.  Here, we need to explicitly
#   remove the tempfile we created.
#------------------------------------------------------------------------------

# =============================================================================
# Timing
# =============================================================================
# The running times of the Python version and the bash script version are
# equivalent:
#
# time ./sh_ms_example.sh > sh.out
# (best of 3)
# --------------------------------
# real	0m42.059s
# user	0m41.970s
# sys	0m0.190s

# time ./py_ms_example.py > py.out
# (best of 3)
# --------------------------------
# real	0m42.000s
# user	0m42.050s
# sys	0m0.230s


# Confirm identical results
# -------------------------
# diff sh.out py.out | wc -l
# 0
