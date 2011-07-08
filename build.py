from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import glob
import os
import sys
sys.argv.append('build_ext')

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

setup(
        ext_modules=exts,
        cmdclass = {'build_ext': build_ext},
    )
