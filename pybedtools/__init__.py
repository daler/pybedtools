from bedtool import bedtool, get_tempdir, set_tempdir, cleanup
import pybedtools
import os

def data_dir():
    """
    Returns the data directory that contains example files for tests and
    documentation
    """
    return os.path.join(os.path.split(pybedtools.__file__)[0], 'test')

def example_bed(bed):
    """
    Return a bed file from the pybedtools examples directory
    """
    fn = os.path.join(data_dir(), bed)
    if not os.path.exists(fn):
        raise ValueError, "%s does not exist" % fn
    return fn

def list_example_beds():
    """
    returns a list of bed files in the examples dir
    """
    return [i for i in os.listdir(data_dir()) if os.path.splitext(i)[-1] == '.bed']
