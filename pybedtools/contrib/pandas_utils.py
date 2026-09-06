import pybedtools


def intersect_dfs(df1, df2, intersect_kwargs=None, read_table_names=None, other_read_table_kwargs=None):
    """
    Intersect two pandas DataFrames using pybedtools.intersect.

    @param df1: file A
    @type df1: pandas.DataFrame
    @param df2: file B
    @type df2: pandas.DataFrame
    @param intersect_kwargs: kwargs passed to pybedtools.intersect
    @type intersect_kwargs: dict
    @param read_table_names: list of column names passed to pandas.read_table (instead of the default ones given by
    pybedtools)
    @type read_table_names: list[str]
    @param other_read_table_kwargs: kwargs passed to pandas.read_table other than `header` and `names`
    @type other_read_table_kwargs: dict
    @return: intersected_df
    @rtype: pandas.DataFrame
    """
    bed1 = pybedtools.BedTool.from_dataframe(df1)
    bed2 = pybedtools.BedTool.from_dataframe(df2)

    intersect_kwargs = {} if intersect_kwargs is None else intersect_kwargs

    intersected_bed = bed1.intersect(bed2, **intersect_kwargs)

    read_table_kwargs = {}
    if read_table_names is not None:
        read_table_kwargs["header"] = None
        read_table_kwargs["names"] = read_table_names
    if other_read_table_kwargs is not None:
        read_table_kwargs.update(other_read_table_kwargs)

    intersected_df = intersected_bed.to_dataframe(**read_table_kwargs)

    return intersected_df
