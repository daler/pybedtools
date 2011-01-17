Quick examples for the impatient
================================

Example 1
---------
Get the sequences of the 100 bp on either side of features.

The :meth:`bedtool.slop()` method automatically downloads chromSizes from
UCSC for dm3 genome, though you could pass your own file using the standard
BEDTools ``slop`` argument of ``g``.  This examples also assumes you have a
local copy of the entire dm3 genome saved as ``dm3.fa``.

::

    from pybedtools import bedtool
    bedtool('in.bed').slop(genome='dm3',l=100,r=100).subtract('in.bed')
    flanking_features.sequence(fi='dm3.fa').save_seqs('flanking.fa')

Example 2
---------
Get values for a 3-way Venn diagram of overlaps::

    import pybedtools
    a = pybedtools.bedtool('a.bed')
    b = pybedtools.bedtool('b.bed')
    c = pybedtools.bedtool('c.bed')

    (a-b-c).count()  # unique to a
    (a+b-c).count()  # in a and b, not c
    (a+b+c).count()  # common to all 
    # ... and so on, for all the combinations.
    
Example 3
---------
Intersections, plus adding track names::


    import pybedtools
    a = pybedtools.bedtool('a.bed')
    a.intersect('b.bed').saveas('a-and-b.bed', trackline="track name='a and b'")
   
