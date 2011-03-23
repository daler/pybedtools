import os
import sys
import subprocess

from bedtool import bedtool, get_tempdir, set_tempdir, cleanup, find_tagged

__version__ = '0.2.2dev'

def check_for_bedtools():
    try:
        p = subprocess.Popen(['intersectBed'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except OSError as err:
        if err.errno == 2:
            print "Please make sure you have installed BEDTools (https://github.com/arq5x/bedtools) and that it's on the path."
            sys.exit(1)

check_for_bedtools()

def data_dir():
    """
    Returns the data directory that contains example files for tests and
    documentation.
    """
    return os.path.join(os.path.dirname(__file__), 'test')

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
