import ez_setup
ez_setup.use_setuptools()

import glob
import os
import sys
from setuptools import setup
from distutils.extension import Extension
try:
    from Cython.Distutils import build_ext

except ImportError:
    sys.stderr.write("""
    ==================================================

    Please install Cython (http://cython.org/),
    which is required to build pybedtools. Usually
    you can do:

        pip install -U cython

    or

        easy_install -U cython

    ==================================================
    """)
    sys.exit(1)

if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__

version_py = os.path.join(os.path.dirname(__file__), 'pybedtools', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','')
sources=["src/bedFile.pyx",
         "src/fileType.pxi",
         "src/gzstream.pxd",
         "pybedtools/cbedtools.cpp"]

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

install_requires = []
if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        install_requires = ['argparse', 'ordereddict']

tests_require = ['nose>=0.11', 'pyyaml']
setup(
        name="pybedtools",
        version=version,
        ext_modules=exts,
        cmdclass = {'build_ext': build_ext},
        tests_require=tests_require,
        install_requires=install_requires,
        extras_require={'test': tests_require},
        packages=['pybedtools',
                  'pybedtools.test',
                  'pybedtools.contrib',
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
                   'pybedtools/scripts/peak_pie.py',
                   'pybedtools/scripts/intersection_matrix.py',
                   'pybedtools/scripts/intron_exon_reads.py',
                   'pybedtools/scripts/pybedtools_demo.py',
                   'pybedtools/scripts/examples/pbt_plotting_example.py',
                   'pybedtools/scripts/pybedtools'],
        author_email="dalerr@niddk.nih.gov",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']

    )
