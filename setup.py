import glob
#from setuptools import setup
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

exts = [ Extension("pybedtools.cbedtools",
              sources=["pybedtools/cbedtools.pyx"] + glob.glob("src/*.cpp"),
              libraries=["stdc++", 'z'],
              include_dirs=["src/"],
              depends = glob.glob("src/*.h"),
              language="c++"), ]

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
        version="0.2.3dev",
        ext_modules=exts,
        packages=['pybedtools','pybedtools.test'],
        author="Ryan Dale",
        description='Wrapper around BEDTools for bioinformatics work',
        long_description=long_description,
        url="none",
        package_data = {'src': ['*.pyx', "*.c", "*.cpp", "*.h", "README.rst"]},
        package_dir = {"pybedtools": "pybedtools"},
        cmdclass = {'build_ext': build_ext},
        author_email="dalerr@niddk.nih.gov",
        test_suite='nose.collector',
    )
