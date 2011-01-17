``pybedtools`` overview and examples
====================================

``pybedtools`` is a Python wrapper for Aaron Quinlan's ``BEDtools``
(http://code.google.com/p/bedtools/), designed to leverage the "genome
algebra" power of ``BEDtools`` from within Python scripts.

See full online documentation at http://daler.github.com/pybedtools (this
README contains just installation instructions and a few examples).


Installation
------------
For the latest version:

    * go to http://github.com/daler/pybedtools 
    * click the Downloads link
    * choose either a .tar.gz or a .zip file, whatever you're 
      comfortable with
    * unzip into a temporary directory
    * from the command line, run::
        
        python setup.py install

You will also need to install ``BEDTools`` -- see
https://github.com/arq5x/bedtools





Two examples
------------

Example: Save a new BED file of intersections between ``a.bed`` and
``b.bed``, adding a track line to the output::

    from pybedtools import bedtool
    a = bedtool('a.bed')
    a.intersect('b.bed').saveas('a-and-b.bed', trackline="track name='a and b' color=128,0,0")

Example: Get values for a 3-way Venn diagram of overlaps.  This
demonstrates operator overloading of ``bedtool`` objects::

    from pybedtools import bedtool
    
    # set up 3 different bedtools
    a = bedtool('a.bed')
    b = bedtool('b.bed')
    c = bedtool('c.bed')
    
    # 
    (a-b-c).count()  # unique to a
    (a+b-c).count()  # in a and b, not c
    (a+b+c).count()  # common to all 
    # ... and so on, for all the combinations.

See the full documentation at http://daler.github.com/pybedtools for more
examples and more information.
