import glob
#from setuptools import setup
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os

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
def get_data_files():
    here = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(here, 'pybedtools', 'test','data')
    data_files = []
    for df in os.listdir(data_dir):
        df = os.path.join(data_dir, df)
        df = os.path.relpath(df, start=here)
        data_files.append(df)
    return data_files

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
        data_files = [('pybedtools/pybedtools/test/data', get_data_files())],
        package_dir = {"pybedtools": "pybedtools"},
        scripts = ['scripts/venn_gchart.py', 'scripts/venn_mpl.py'],
        cmdclass = {'build_ext': build_ext},
        author_email="dalerr@niddk.nih.gov",
        test_suite='nose.collector',
    )
