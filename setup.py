import ez_setup
ez_setup.use_setuptools()

import glob
import os
import sys
from setuptools import setup
from distutils.extension import Extension

# optional cython
try:
  from Cython.Distutils import build_ext
except ImportError:
  from distutils.command import build_ext as _build_ext
  class build_ext(_build_ext.build_ext):

      description = "change pyx files to corresponding .c/.cpp (fallback when cython is not installed)"

      def build_extensions(self):
          # First, sanity-check the 'extensions' list
          self.check_extensions_list(self.extensions)
          
          for extension in self.extensions:
              iscpp = extension.language and extension.language.lower() == 'c++'
              target_ext = '.cpp' if iscpp else '.c'

              patchedsrc = []
              for source in extension.sources:
                (root, ext) = os.path.splitext(source)
                if ext == '.pyx':
                  patchedsrc.append(root + target_ext)
                else:
                  patchedsrc.append(source)

              extension.sources = patchedsrc
              self.build_extension(extension)
  

if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__

version_py = os.path.join(os.path.dirname(__file__), 'pybedtools', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','')
sources=["src/bedFile.cpp",
         "src/fileType.cpp",
         "src/gzstream.cpp",
         "pybedtools/cbedtools.pyx"]


exts = [ Extension("pybedtools.cbedtools",
                   sources=sources,
                   libraries=["stdc++", 'z'],
                   library_dirs=["src/"],
                   include_dirs=["src/"],
                   language="c++"),

         Extension('pybedtools.featurefuncs',
                   sources=sources + ["pybedtools/featurefuncs.pyx"],
                   libraries=["stdc++", 'z'],
                   include_dirs=["src/"],
                   library_dirs=["src/"],
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

tests_require = ['nose>=0.11', 'pyyaml']
setup(
        cmdclass= {'build_ext': build_ext},
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
                   'pybedtools/scripts/pybedtools',
                   'pybedtools/scripts/examples/pbt_plotting_example.py'],
        author_email="dalerr@niddk.nih.gov",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']

    )
