.. include:: includeme.rst

Intro
=====


This tutorial assumes that 

1. You know how to use BEDTools_ (if not, check out the 
   `BEDTools documentation`_)
2. You know how to use Python (if not, check out some 
   tutorials like `Learn Python the Hard Way`_)


A brief note on conventions
---------------------------
Throughout this documentation I've tried to use consistent typography, as
follows:

* Python variables and arguments are shown in italics: *s=True*
* Files look like this: :file:`filename.bed`
* Methods, which are often linked to documentation look like this:
  :meth:`BedTool.merge`.
* BEDTools_ programs look like this: ``intersectBed``.
* Arguments that are passed to BEDTools_ programs, as if you were on the
  command line, look like this: ``-d``.
* The ">>>" in the examples below indicates a Python interpreter prompt and
  means to type the code into an interactive Python interpreter like IPython_
  (don't type the >>>) 

Onward!
