from bedtool import bedtool, get_tempdir, set_tempdir, cleanup
import pybedtools
import os

__version__ = '0.2.2dev'

def data_dir():
    """
    Returns the data directory that contains example files for tests and
    documentation.
    """
    return os.path.join(os.path.split(pybedtools.__file__)[0], 'test')

def example_bed_fn(bed):
    """
    Return a bed file from the pybedtools examples directory.  Use
    :func:`list_example_beds` to see a list of files that are included.
    """
    fn = os.path.join(data_dir(), bed)
    if not os.path.exists(fn):
        raise ValueError, "%s does not exist" % fn
    return fn

def example_bedtool(fn):
    """
    Return a bedtool using a bed file from the pybedtools examples directory.
    Use :func:`list_example_beds` to see a list of files that are included.
    """
    fn = os.path.join(data_dir(), fn)
    if not os.path.exists(fn):
        raise ValueError, "%s does not exist" % fn
    return bedtool(fn)

def list_example_beds():
    """
    Returns a list of bed files in the examples dir.  Choose one and pass
    it to :func:`example_bed` to get the full path to an example BED file.

    Example usage::

        >>> choices = list_example_beds()
        >>> bedfn = example_bed(choices[0])
        >>> mybedtool = bedtool(bedfn)


    """
    return sorted([i for i in os.listdir(data_dir()) if os.path.splitext(i)[-1] == '.bed'])
