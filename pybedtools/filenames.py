"""
Provides access to example files and keeps track of all temp files created
during a Python session.
"""
import os

TEMPFILES = []


def data_dir():
    """
    Returns the data directory that contains example files for tests and
    documentation.
    """
    return os.path.join(os.path.dirname(__file__), 'test', 'data')


def example_filename(fn):
    """
    Return a bed file from the pybedtools examples directory.  Use
    func:`list_example_files` to see a list of files that are included.
    """
    fn = os.path.join(data_dir(), fn)
    if not os.path.exists(fn):
        raise ValueError("%s does not exist" % fn)
    return fn


def list_example_files():
    """
    Returns a list of files in the examples dir.  Choose one and pass it to
    :func:`example_filename` to get the full path to an example file.

    Example usage:

        >>> choices = list_example_files()
        >>> assert 'a.bed' in choices
        >>> bedfn = example_filename('a.bed')
        >>> mybedtool = BedTool(bedfn)

    """
    candidate_fns = os.listdir(data_dir())
    exts = ('.bed', '.gff', '.gtf', '.bed.gz', '.bam', '.gff.gz')
    valid_fns = [f for f in candidate_fns if f.endswith(exts)]
    return sorted(valid_fns)
