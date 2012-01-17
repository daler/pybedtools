import ez_setup
ez_setup.use_setuptools()

import glob
import os
import sys
from setuptools import setup
from distutils.extension import Extension

if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__

version_py = os.path.join(os.path.dirname(__file__), 'pybedtools', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','')
sources=["src/bedFile.cpp",
         "src/fileType.cpp",
         "src/gzstream.cpp",
         "pybedtools/cbedtools.cpp"]
exts = [ Extension("pybedtools.cbedtools",
                   sources=sources,
                   libraries=["stdc++", 'z'],
                   library_dirs=["src/"],
                   include_dirs=["src/"],
                   language="c++"),

         Extension('pybedtools.featurefuncs',
                   sources=sources + ["pybedtools/featurefuncs.cpp"],
                   libraries=["stdc++", 'z'],
                   include_dirs=["src/"],
                   library_dirs=["src/"],
                   language="c++"),

         Extension('pybedtools._Window',
                    sources=['pybedtools/_Window.c'],),
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

tests_require = ['nose>=0.11', 'pyyaml']
setup(
        name="pybedtools",
        version=version,
        ext_modules=exts,
        install_requires=['argparse'],
        tests_require=tests_require,
        setup_requires=tests_require,
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
                   'pybedtools/scripts/pybedtools'],
        author_email="dalerr@niddk.nih.gov",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']

    )
