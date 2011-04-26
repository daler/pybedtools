.. include:: includeme.rst

More documentation
------------------
For more info, see the :ref:`topical`.

.. doctest::
    :hide:

    Gotta clean up all the files created over the course of the tutorial...

    >>> fns_to_remove = ['a-with-b.bed', 'a-with-b-merged.bed', 'hg19.genome', 'intersection-of-a-and-b.bed','middle-100-bp.bed','shared_merged.bed']
    >>> for fn in fns_to_remove:
    ...     if os.path.exists(fn):
    ...         os.unlink(fn)
