import ez_setup
ez_setup.use_setuptools()

import glob
import os
import sys

try:
    from Cython.Distutils import build_ext

    # Work around setuptools bug by providing a fake Pyrex, as is done here:
    #    http://codespeak.net/svn/lxml/trunk/
    #
    # This seems to solve errors where the compiler is looking for .c files
    # (which don't exist for these C++ wrappers!)
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "fake_pyrex"))

except ImportError:
    sys.stderr.write('''
-------------------------------------------------------------------------------
Please install Cython 0.14.1 or greater (http://docs.cython.org) before running
setup.py. You can use:

    $ easy_install cython

-------------------------------------------------------------------------------
''')
    sys.exit(1)

# Strangely, this needs to come AFTER the Cython import
from setuptools import setup
from setuptools.extension import Extension

version_py = os.path.join(os.path.dirname(__file__), 'pybedtools', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','')

exts = [ Extension("pybedtools.cbedtools",
                   sources=["pybedtools/cbedtools.pyx",
                            "pybedtools/cbedtools.pxi",
                            "pybedtools/cbedtools.pxd"] \
                            + glob.glob("src/*.cpp"),
                   libraries=["stdc++", 'z'],
                   include_dirs=["src/"],
                   depends = glob.glob("src/*.h"),
                   language="c++"),

         Extension('pybedtools.featurefuncs',
                   sources=["pybedtools/featurefuncs.pyx",
                            "pybedtools/cbedtools.pyx",
                            "pybedtools/cbedtools.pxi",
                            "pybedtools/cbedtools.pxd"] \
                            + glob.glob("src/*.cpp"),
                   libraries=["stdc++", 'z'],
                   include_dirs=["src/"],
                   depends = glob.glob("src/*.h"),
                   language="c++"),

         Extension('pybedtools._Window',
                    sources=['pybedtools/_Window.pyx'],),
        ]

long_description = """
``pybedtools`` is a Python extension of Aaron Quinlan's BEDtools suite
(http://code.google.com/p/bedtools/), used for comparing genomic features.

``pybedtools`` allows you to intuitively call BEDtools programs from within
Python without writing awkward system calls, and allows you to manipulate data
on the file level as well as on the individual feature level.

Development version can be found on github:

    http://github.com/daler/pybedtools

and see full documentation and tutorial at:

    http://packages.python.org/pybedtools

"""
tests_require = ['nose']
setup(
        name="pybedtools",
        version=version,
        ext_modules=exts,
        install_requires=['argparse','cython'],
        tests_require=tests_require,
        extras_require={'test': tests_require},
        packages=['pybedtools',
                  'pybedtools.test',
                  'pybedtools.scripts',
                  'pybedtools.test.data'],
        author="Ryan Dale",
        description='Wrapper around BEDTools for bioinformatics work',
        long_description=long_description,
        url="none",
        package_data = {'pybedtools':["test/data/*",
                                      "*.pyx",
                                      "*.pxi",
                                      "*.pxd",
                                      "*.cpp"]
                        },
        package_dir = {"pybedtools": "pybedtools"},
        scripts = ['pybedtools/scripts/venn_gchart.py',
                   'pybedtools/scripts/venn_mpl.py',
                   'pybedtools/scripts/annotate.py',
                   'pybedtools/scripts/intron_exon_reads.py',
                   'pybedtools/scripts/pybedtools_demo.py'],
        cmdclass = {'build_ext': build_ext},
        author_email="dalerr@niddk.nih.gov",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']

    )
