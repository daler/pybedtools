import os
import sys
import glob

usage = """

Cython .pyx files as well as the created .cpp files are included in the source
distribution. The following information is useful for developers working on the
Cython source code.

    Install in development mode. Will cythonize .pyx files first if needed.

        python setup.py cythonize develop

    Rebuild .cpp files from .pyx, and then stop:

        python setup.py cythonize

    Build extensions from existing .cpp

        python setup.py build_ext

    Source distribution:

        python setup.py clean cythonize sdist
"""

try:
    from Cython.Build import build_ext as cython_build_ext
    from Cython.Build import cythonize
    HAVE_CYTHON = True
except ImportError:
    HAVE_CYTHON = False

if '--usage' in sys.argv:
    print(usage)
    sys.exit(0)

elif 'cythonize' in sys.argv:
    USE_CYTHON = True

else:
    USE_CYTHON = False

if USE_CYTHON and not HAVE_CYTHON:
    raise ValueError(
        '''
        Cython could not be found. Please install Cython and try again.
        ''')


# Try bootstrapping setuptools if it doesn't exist. This is for using the
# `develop` command, which is very useful for in-place development work.
try:
    import pkg_resources
    try:
        pkg_resources.require("setuptools>=0.6c5")
    except pkg_resources.VersionConflict:
        from ez_setup import use_setuptools
        use_setuptools(version="0.6c5")
    from setuptools import setup, Command
except ImportError:
    sys.exit(
        'pybedtools uses setuptools '
        '(https://packaging.python.org/installing/) '
        'for installation but setuptools was not found')

curdir = os.path.abspath(os.path.dirname(__file__))

# These imports need to be here; setuptools needs to be imported first.
from distutils.extension import Extension  # noqa: E402
from distutils.command.build import build  # noqa: E402
from distutils.command.build_ext import build_ext  # noqa: E402
from distutils.command.sdist import sdist  # noqa: E402
import distutils.log


MAJ = 0
MIN = 8
REV = 2
VERSION = '%d.%d.%d' % (MAJ, MIN, REV)


class CleanCommand(Command):
    """
    Custom distutils command to clean the various files created by cython.

    E.g.,

        pybedtools/featurefuncs.cpp
        pybedtools/cbedtools.cpp
        pybedtools/cbedtools.cpython-36m-x86_64-linux-gnu.so
        pybedtools/featurefuncs.cpython-36m-x86_64-linux-gnu.so

    """
    user_options = []

    def initialize_options(self):
        self._clean_me = []
        self._clean_trees = []

        # Add files to be protected here
        self._clean_exclude = ['bedFile.cpp', 'fileType.cpp', 'gzstream.cpp']

        for root, dirs, files in list(os.walk('pybedtools')):
            for f in files:
                if f in self._clean_exclude:
                    continue
                if os.path.splitext(f)[-1] in ('.pyc', '.so', '.o', '.pyo',
                                               '.pyd', '.c', '.cpp', '.cxx',
                                               '.orig'):
                    self._clean_me.append(os.path.join(root, f))
            for d in dirs:
                if d == '__pycache__':
                    self._clean_trees.append(os.path.join(root, d))

        for d in ('build',):
            if os.path.exists(d):
                self._clean_trees.append(d)

    def finalize_options(self):
        pass

    def run(self):
        for clean_me in self._clean_me:
            try:
                print('removing', clean_me)
                os.unlink(clean_me)
            except Exception:
                pass
        for clean_tree in self._clean_trees:
            try:
                import shutil
                print('removing directory', clean_tree)
                shutil.rmtree(clean_tree)
            except Exception:
                pass


class CythonBuildExt(build_ext):
    """
    Subclass build_ext to get clearer report if Cython is necessary.
    """
    def build_extensions(self):
        for ext in self.extensions:
            cythonize(ext)
        build_ext.build_extensions(self)


class Cythonize(Command):
    """
    Generate .cpp files and then stop
    """
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        cythonize(extensions)


def find_missing_src(ext):
    """
    Check whether sources needed to build extension *ext* are present.

    Will return a list of missing source files that might be obtained by running cythonize.

    Will raise a FileNotFoundError if some missing .cpp source files do not have a corresponding .pyx.
    """
    missing_src = []
    for src in ext.sources:
        if not os.path.exists(src):
            # raise Exception(
            #     """
            #     Cython-generated file '%s' not found.

            #     Please install Cython and run

            #         python setup.py cythonize

            #     """ % src)
            (root, extn) = os.path.splitext(src)
            if extn == ".cpp":
                alt_src = root + ".pyx"
                if not os.path.exists(alt_src):
                    raise FileNotFoundError(
                        "Files %s and %s not found." % (src, alt_src))
                missing_src.append(src)
            else:
                raise FileNotFoundError(
                    "File %s not found." % src)
    return missing_src


class InformativeBuildExt(build_ext):
    def build_extensions(self):
        for ext in self.extensions:
            missing_src = find_missing_src(ext)
            if missing_src:
                if not HAVE_CYTHON:
                    raise ValueError(
                        '''
                        Cython could not be found.
                        Please install Cython and try again.
                        ''')
                self.announce(
                    "Trying to generate the following missing files:\n%s" % "\n".join(missing_src),
                    level=distutils.log.INFO)
                for src in missing_src:
                    assert src in ext.sources
                    (root, extn) = os.path.splitext(src)
                    assert extn == ".cpp"
                    cythonize(root + ".pyx")
                still_missing = find_missing_src(ext)
                if still_missing:
                    raise ValueError(
                        '''
                        Some source files are missing to build an extension.

                        %s''' % "\n".join(still_missing))
        build_ext.build_extensions(self)


class SDist(sdist):

    def run(self):
        cythonize(extensions)
        for ext in extensions:
            for src in ext.sources:
                if not os.path.exists(src):
                    raise Exception(
                        "Cython-generated file '{0}' not found. "
                        "Run 'python setup.py --usage' for details.".format(src, usage))
        sdist.run(self)


EXT = '.pyx' if USE_CYTHON else '.cpp'

extensions = [
    Extension(
        'pybedtools.cbedtools',
        depends=glob.glob('pybedtools/include/*h'),
        libraries=['stdc++', 'z'],
        include_dirs=['pybedtools/include/'],
        sources=['pybedtools/cbedtools' + EXT] + glob.glob('pybedtools/include/*.cpp'),
        language='c++'),

    Extension(
        'pybedtools.featurefuncs',
        depends=glob.glob('pybedtools/include/*h'),
        libraries=['stdc++', 'z'],
        include_dirs=['pybedtools/include/'],
        sources=['pybedtools/featurefuncs' + EXT] + glob.glob('pybedtools/include/*.cpp'),
        language='c++'),
]


cmdclass = {
    'clean': CleanCommand,
    'build': build,
    'sdist': SDist,
}

if USE_CYTHON:
    cmdclass['build_ext'] = cython_build_ext
    cmdclass['cythonize'] = Cythonize
else:
    cmdclass['build_ext'] = InformativeBuildExt

if __name__ == "__main__":

    with open(os.path.join(curdir, 'pybedtools/version.py'), 'w') as fout:
        fout.write(
            "\n".join(["",
                       "# THIS FILE IS GENERATED FROM SETUP.PY",
                       "version = '{version}'",
                       "__version__ = version"]).format(version=VERSION)
        )

    README = open(os.path.join(curdir, "README.rst")).read()

    setup(
        name='pybedtools',
        maintainer='Ryan Dale',
        version=VERSION,
        ext_modules=extensions,
        maintainer_email='ryan.dale@nih.gov',
        description='Wrapper around BEDTools for bioinformatics work',
        license='GPLv2',
        url='https://github.com/daler/pybedtools',
        download_url='',
        long_description=README,
        zip_safe=False,
        setup_requires=[],
        install_requires=['six', 'pysam'],
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Programming Language :: Python',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.6',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.3',
            'Topic :: Software Development :: Libraries :: Python Modules',
        ],
        cmdclass=cmdclass,
        packages=['pybedtools',
                  'pybedtools.test',
                  'pybedtools.contrib',
                  'pybedtools.test.data'],
        package_data={'pybedtools': ["test/data/*",
                                     "*.pyx",
                                     "*.pxi",
                                     "*.pxd",
                                     "*.cxx",
                                     "*.c",
                                     "*.cpp",
                                     "*.h"],
                      'src': ['src/*'],
                      },
        include_package_data=True,
    )
