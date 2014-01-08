Overview
--------

.. image:: https://travis-ci.org/daler/pybedtools.png?branch=master
    :target: https://travis-ci.org/daler/pybedtools

`pybedtools` is a Python wrapper for Aaron Quinlan's `BEDtools` programs
(https://github.com/arq5x/bedtools), which are widely used for genomic interval
manipulation or "genome algebra".  `pybedtools` extends `BEDTools` by offering
feature-level manipulations from with Python. See full online documentation,
including installation instructions, at http://pythonhosted.org/pybedtools/.

Why `pybedtools`?
-----------------

Here is an example to get the names of genes that are <5 kb away from
intergenic SNPs:

.. code-block:: python

    from pybedtools import BedTool

    snps = BedTool('snps.bed.gz')
    genes = BedTool('hg19.gff')

    intergenic_snps = snps.subtract(genes)
    nearby = genes.closest(intergenic_snps, d=True, stream=True)

    for gene in nearby:
        if int(gene[-1]) < 5000:
            print gene.name

Useful features shown here include:

* support for all BEDTools-supported formats (here gzipped BED and GFF)
* wrapping of all BEDTools programs and arguments (here, `subtract` and `closest` and passing
  the `-d` flag to `closest`);
* streaming results (like Unix pipes, here specified by `stream=True`)
* iterating over results while accessing feature data by index or by attribute
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

See the `Shell script comparison <http://pythonhosted.org/pybedtools/sh-comparison.html>`_ in the docs
for more details on this comparison, or keep reading the full documentation at
http://pythonhosted.org/pybedtools/index.html.


