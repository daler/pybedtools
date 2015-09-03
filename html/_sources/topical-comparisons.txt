.. include:: includeme.rst

Comparisons
===========
Sometimes it is useful to be able to do quick comparisons between features to
see if they overlap or if they are to the left or to the right.  Comparsion
operators (`<`, `<=`, `==`, `=>`, `>`) are defined for intervals.  Note that
these comparsions **ignore strand**; if you need more control then it's
probably better to write a quick one-off function to do the comparisons.

In general, `>` and `<` are True if the features are completely
separate from each other; if they overlap then `>=` and `<=` are True as well.
Nested features are not comparable, so a NotImplementedError will be raised.

It's probably easiest to describe these operators "ASCII-graphically"::

    # a == b, a >= b, a <= b
    a ---------
    b ---------

    # a < b, a <= b
    a ----
    b       -----

    # a <= b
    a ----
    b     -----  (book-ended)

    # a >= b
    a     -----
    b ----      (book-ended)

    # a > b, a >= b
    a       ------
    b ----

    # a >= b
    a  ------------
    b  ---------

    # a >= b
    a   -----------
    b -------------

    # a <= b
    a -------------
    b   -----------

    # a <= b
    a  ---------
    b  ------------

    # a <= b
    a -----------
    b        -----------

    # a >= b
    a        -----------
    b -----------

    # undefined!
    a    ----
    b -----------

    # undefined!
    a -----------
    b    ----

    # a <= b
    a -----------
    b           -

    # a >= b
    a           -
    b -----------

    # a == b, a <= b, a >= b
    a -
    b -  (starts and stops are identical for all features)
