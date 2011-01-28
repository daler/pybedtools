#from distutils.core import setup
from setuptools import setup, find_packages

long_description = """
``pybedtools`` is a wrapper around Aaron Quinlan's BEDtools suite
(http://code.google.com/p/bedtools/), used for comparing genomic features.

``pybedtools`` allows you to intuitively call BEDtools programs from within
Python without writing awkward system calls.

Development version, as well as documentation, can be found on github:

    http://github.com/daler/pybedtools

"""

setup( 
        name="pybedtools",
        version="0.2.2dev",
        packages=['pybedtools','pybedtools.test'],
        author="Ryan Dale",
        description='Wrapper around BEDTools for bioinformatics work',
        long_description=long_description,
        url="none",
        author_email="dalerr@niddk.nih.gov"
    )
