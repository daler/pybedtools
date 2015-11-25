Overview
--------

.. image:: https://travis-ci.org/daler/pybedtools.png?branch=master
    :target: https://travis-ci.org/daler/pybedtools

.. image:: https://badge.fury.io/py/pybedtools.svg?style=flat
    :target: http://badge.fury.io/py/pybedtools

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
    :target: http://bioconda.github.io

The `BEDTools suite of programs <http://bedtools.readthedocs.org/>`_ is widely
used for genomic interval manipulation or "genome algebra".  `pybedtools` wraps
and extends BEDTools and offers feature-level manipulations from within
Python.

See full online documentation, including installation instructions, at
http://daler.github.io/pybedtools/.

Why `pybedtools`?
-----------------

Here is an example to get the names of genes that are <5 kb away from
intergenic SNPs:

.. code-block:: python

    from pybedtools import BedTool

    snps = BedTool('snps.bed.gz')  # [1]
    genes = BedTool('hg19.gff')    # [1]

    intergenic_snps = snps.subtract(genes)                       # [2]
    nearby = genes.closest(intergenic_snps, d=True, stream=True) # [2, 3]

    for gene in nearby:             # [4]
        if int(gene[-1]) < 5000:    # [4]
            print gene.name         # [4]

Useful features shown here include:

* `[1]` support for all BEDTools-supported formats (here gzipped BED and GFF)
* `[2]` wrapping of all BEDTools programs and arguments (here, `subtract` and `closest` and passing
  the `-d` flag to `closest`);
* `[3]` streaming results (like Unix pipes, here specified by `stream=True`)
* `[4]` iterating over results while accessing feature data by index or by attribute
  access (here `[-1]` and `.name`).

In contrast, here is the same analysis using shell scripting.  Note that this
requires knowledge in Perl, bash, and awk.  The run time is identical to the
`pybedtools` version above:

.. code-block:: bash

    snps=snps.bed.gz
    genes=hg19.gff
    intergenic_snps=/tmp/intergenic_snps

    snp_fields=`zcat $snps | awk '(NR == 2){print NF; exit;}'`
    gene_fields=9
    distance_field=$(($gene_fields + $snp_fields + 1))

    intersectBed -a $snps -b $genes -v > $intergenic_snps

    closestBed -a $genes -b $intergenic_snps -d \
    | awk '($'$distance_field' < 5000){print $9;}' \
    | perl -ne 'm/[ID|Name|gene_id]=(.*?);/; print "$1\n"'

    rm $intergenic_snps

See the `Shell script comparison <http://daler.github.io/pybedtools/sh-comparison.html>`_ in the docs
for more details on this comparison, or keep reading the full documentation at
http://daler.github.io/pybedtools.


