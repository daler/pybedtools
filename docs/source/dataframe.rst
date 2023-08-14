.. include:: includeme.rst

.. _saveresults:

Exporting to a dataframe
==================

If you want to export the results as a dataframe for more analysis, use
the :meth:`BedTool.to_dataframe` method to export to a pandas dataframe or the :meth:`BedTool.to_polars_dataframe` method to export to a polars dataframe. This method also lets you optionally specify column names for the dataframes instead of the default columns names that pybedtools uses. You can use the same arguments you would normally use while reading a file into a pandas (`names=`) or polars (`new_columns=`) dataframe. By default, pybedtools assumes that there is no header line in the bed file. If your bed file already has names in the first row, you can set the `disable_auto_names` argument to `False`.

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> import pandas
    >>> import polars
    >>> a = pybedtools.example_bedtool('a.bed')
    <BLANKLINE>
    
    >>> pandas_df = a.to_dataframe()
    >>> print(pandas_df)
        chrom   start    end        name   score  strand
    0    chr1       1    100    feature1       0       +
    1    chr1     100    200    feature2       0       +
    2    chr1     150    500    feature3       0       -
    3    chr1     900    950    feature4       0       +
    <BLANKLINE>
    
    >>> polars_df = a.to_polars_dataframe()
    >>> print(polars_df)
    ——————————————————————————————————————————————————————
    │ chrom ┆ start ┆ end   ┆ name      ┆ score ┆ strand │
    │ ---   ┆ ---   ┆ ---   ┆ ---       ┆ ---   ┆ ---    │
    │ str   ┆ i64   ┆ i64   ┆ str       ┆ i64   ┆ str    │
    ══════════════════════════════════════════════════════
    │ chr1  ┆ 1     ┆ 100   ┆ feature1  ┆ 0     ┆ +      │
    │ chr1  ┆ 100   ┆ 200   ┆ feature2  ┆ 0     ┆ +      │
    │ chr1  ┆ 150   ┆ 500   ┆ feature3  ┆ 0     ┆ -      │
    │ chr1  ┆ 900   ┆ 950   ┆ feature4  ┆ 0     ┆ +      │
    ——————————————————————————————————————————————————————
    <BLANKLINE>

You can also generate a :class:`BedTool` object from a pandas or polars dataframe using the  :meth:`BedTool.from_dataframe` or :meth:`BedTool.from_polars_dataframe` method respectively.
