import glob
import os
import sys
from setuptools import setup
#from distutils.core import setup
from distutils.extension import Extension

import os.path

def no_cythonize(extensions, **_ignore):
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in ('.pyx', '.py'):
                if extension.language == 'c++':
                    ext = '.cpp'
                else:
                    ext = '.c'
                sfile = path + ext
            sources.append(sfile)
        extension.sources[:] = sources
    return extensions


try:
    from Cython.Build import cythonize
    USE_CYTHON = True

except ImportError:
    USE_CYTHON = False

ext = '.pyx' if USE_CYTHON else '.cpp'


version_py = os.path.join(os.path.dirname(__file__), 'pybedtools', 'version.py')
version = open(version_py).read().strip().split(' = ')[-1].replace('"','')

extensions = [ Extension("pybedtools.cbedtools",
                   sources=["pybedtools/cbedtools.pyx",
                            ] \
                            + glob.glob("src/*.cpp"),
                   libraries=["stdc++", 'z'],
                   include_dirs=["src/"],
                   depends = glob.glob("src/*.h"),
                   language="c++"),

         Extension('pybedtools.featurefuncs',
                   sources=["pybedtools/featurefuncs.pyx",
                            "pybedtools/cbedtools.pyx",
                            ] \
                            + glob.glob("src/*.cpp"),
                   libraries=["stdc++", 'z'],
                   include_dirs=["src/"],
                   depends = glob.glob("src/*.h"),
                   language="c++"),

         Extension('pybedtools._Window',
                    sources=['pybedtools/_Window.pyx'],),
        ]


if USE_CYTHON:
    extensions = cythonize(extensions)
else:
    extensions = no_cythonize(extensions)

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

install_requires = ['pysam']
if sys.version_info[0] == 2:
    if sys.version_info[1] < 7:
        install_requires = ['argparse', 'ordereddict']

tests_require = ['nose>=0.11', 'pyyaml']
setup(
        name="pybedtools",
        version=version,
        ext_modules=extensions,
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
