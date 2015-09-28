"""
Much help from statsmodels, numpy, scipy, pandas, pyzmq, lxml
http://pages.uoregon.edu/cfulton/posts/building_python_modules.html
https://github.com/statsmodels/statsmodels/blob/master/setup.py

The idea is to provide just the .cpp files in the distribution, and only run
Cython when creating an sdist.
"""

import os
from os.path import relpath, join as pjoin
import sys
import subprocess
import re
from distutils.version import StrictVersion
import glob

no_frills = (len(sys.argv) >= 2 and ('--help' in sys.argv[1:] or
                                     sys.argv[1] in ('--help-commands',
                                                     'egg_info', '--version',
                                                     'clean')))

# try bootstrapping setuptools if it doesn't exist
try:
    import pkg_resources
    try:
        pkg_resources.require("setuptools>=0.6c5")
    except pkg_resources.VersionConflict:
        from ez_setup import use_setuptools
        use_setuptools(version="0.6c5")
    from setuptools import setup, Command, find_packages
    _have_setuptools = True
except ImportError:
    # no setuptools installed
    from distutils.core import setup, Command
    _have_setuptools = False

if _have_setuptools:
    setuptools_kwargs = {"zip_safe": False,
                         "test_suite": "nose.collector"}
else:
    setuptools_kwargs = {}
    if sys.version_info[0] >= 3:
        sys.exit("Need setuptools to install pybedtools for Python 3.x")


curdir = os.path.abspath(os.path.dirname(__file__))
README = open(pjoin(curdir, "README.rst")).read()


DISTNAME = 'pybedtools'
DESCRIPTION = 'Wrapper around BEDTools for bioinformatics work'
LONG_DESCRIPTION = README
MAINTAINER = 'Ryan Dale'
MAINTAINER_EMAIL = 'dalerr@niddk.nih.gov'
URL = 'https://github.com/daler/pybedtools'
LICENSE = 'GPLv2'
DOWNLOAD_URL = ''


# These imports need to be here; setuptools needs to be imported first.
from distutils.extension import Extension
from distutils.command.build import build
from distutils.command.build_ext import build_ext as _build_ext


class build_ext(_build_ext):
    def build_extensions(self):

        # Pybedtools doesn't need NumPy. Comment out adding the include files,
        # but keep this wrapped function so that later calls still work.
        #
        # numpy_incl = pkg_resources.resource_filename('numpy', 'core/include')
        # for ext in self.extensions:
        #     if (hasattr(ext, 'include_dirs') and
        #             not numpy_incl in ext.include_dirs):
        #         ext.include_dirs.append(numpy_incl)
        #
        _build_ext.build_extensions(self)


def generate_cython():
    cwd = os.path.abspath(os.path.dirname(__file__))
    print("Cythonizing sources")
    p = subprocess.call([sys.executable,
                         os.path.join(cwd, 'tools', 'cythonize.py'),
                         'pybedtools'],
                        cwd=cwd)
    if p != 0:
        raise RuntimeError("Running cythonize failed!")


def strip_rc(version):
    return re.sub(r"rc\d+$", "", version)


def check_dependency_versions(min_versions):
    """
    Don't let pip/setuptools do this all by itself.  It's rude.
    For all dependencies, try to import them and check if the versions of
    installed dependencies match the minimum version requirements.  If
    installed but version too low, raise an error.  If not installed at all,
    return the correct ``setup_requires`` and ``install_requires`` arguments to
    be added to the setuptools kwargs.  This prevents upgrading installed
    dependencies like numpy (that should be an explicit choice by the user and
    never happen automatically), but make things work when installing into an
    empty virtualenv for example.
    """
    setup_requires = []
    install_requires = ['six']

    if 'pysam' in min_versions:
        try:
            from pysam import __version__ as pysam_version
        except ImportError:
            install_requires.append('pysam')
        else:
            if not (StrictVersion(pysam_version) >= min_versions['pysam']):
                raise ImportError("Pysam version is %s. Requires >= %s" %
                                  (pysam_version, min_versions['pysam']))

    if 'pandas' in min_versions:
        try:
            from pandas import __version__ as pandas_version
        except ImportError:
            install_requires.append('pandas')
        else:
            if not (StrictVersion(pandas_version) >= min_versions['pandas']):
                raise ImportError("pandas version is %s. Requires >= %s" %
                                  (pandas_version, min_versions['pandas']))

    return setup_requires, install_requires


MAJ = 0
MIN = 7
REV = 1
ISRELEASED = True
VERSION = '%d.%d.%d' % (MAJ, MIN, REV)


def git_version():
    """Return the git revision as a string"""
    def _minimal_ext_cmd(cmd):
        """construct minimal environment"""
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v

        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(" ".join(cmd), stdout=subprocess.PIPE, env=env,
                               shell=True).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"
    return GIT_REVISION


def write_version_py(filename=pjoin(curdir, 'pybedtools/version.py')):
    cnt = "\n".join(["",
                     "# THIS FILE IS GENERATED FROM SETUP.PY",
                     "short_version = '%(version)s'",
                     "version = '%(version)s'",
                     "full_version = '%(full_version)s'",
                     "git_revision = '%(git_revision)s'",
                     "release = %(isrelease)s", "",
                     "__version__ = version",
                     "if not release:",
                     "    version = full_version"])

    # Adding the git rev number needs to be done inside write_version_py(),
    # otherwise the import of numpy.version messes up the build under Python 3.
    FULLVERSION = VERSION
    dowrite = True
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists(filename):
        # must be a source distribution, use existing version file
        try:
            from pybedtools.version import git_revision as GIT_REVISION
        except ImportError:
            dowrite = False
            GIT_REVISION = "Unknown"
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev0+' + GIT_REVISION[:7]

    if dowrite:
        with open(filename, 'w') as a:
            a.write(cnt % {'version': VERSION,
                           'full_version': FULLVERSION,
                           'git_revision': GIT_REVISION,
                           'isrelease': str(ISRELEASED)})


class CleanCommand(Command):
    """Custom distutils command to clean the .so and .pyc files."""

    user_options = [("all", "a", "")]

    def initialize_options(self):
        self.all = True
        self._clean_me = []
        self._clean_trees = []
        self._clean_exclude = []

        for root, dirs, files in list(os.walk('pybedtools')):
            for f in files:
                if f in self._clean_exclude:
                    continue
                if os.path.splitext(f)[-1] in ('.pyc', '.so', '.o', '.pyo',
                                               '.pyd', '.c', '.cpp', '.cxx',
                                               '.orig'):
                    self._clean_me.append(pjoin(root, f))
            for d in dirs:
                if d == '__pycache__':
                    self._clean_trees.append(pjoin(root, d))

        for d in ('build',):
            if os.path.exists(d):
                self._clean_trees.append(d)

    def finalize_options(self):
        pass

    def run(self):
        for clean_me in self._clean_me:
            try:
                os.unlink(clean_me)
            except Exception:
                pass
        for clean_tree in self._clean_trees:
            try:
                import shutil
                shutil.rmtree(clean_tree)
            except Exception:
                pass


class CheckingBuildExt(build_ext):
    """Subclass build_ext to get clearer report if Cython is necessary."""

    def check_cython_extensions(self, extensions):
        for ext in extensions:
            for src in ext.sources:
                if not os.path.exists(src):
                    raise Exception("""Cython-generated file '%s' not found.
        Cython is required to compile pybedtools from a development branch.
        Please install Cython or download a source release of pybedtools.
                """ % src)

    def build_extensions(self):
        self.check_cython_extensions(self.extensions)
        build_ext.build_extensions(self)


class DummyBuildSrc(Command):
    """ numpy's build_src command interferes with Cython's build_ext.
    """
    user_options = []

    def initialize_options(self):
        self.py_modules_dict = {}

    def finalize_options(self):
        pass

    def run(self):
        pass


cmdclass = {'clean': CleanCommand,
            'build': build}

cmdclass["build_src"] = DummyBuildSrc
cmdclass["build_ext"] = CheckingBuildExt


ext_data = dict(
    cbedtools={
        'name': 'pybedtools/cbedtools.cxx',
        'depends': glob.glob('src/*.h'),
        'libraries': ["stdc++", 'z'],
        'include_dirs': ['src/'],
        'sources': glob.glob("src/*.cpp"),
        'language': 'c++'},

    featurefuncs={
        'name': 'pybedtools/featurefuncs.cxx',
        'depends': glob.glob('src/*.h'),
        'libraries': ["stdc++", 'z'],
        'include_dirs': ['src/'],
        'sources': glob.glob("src/*.cpp"),
        'language': 'c++'},
)

extensions = []
for name, data in ext_data.items():
    data['sources'] = data.get('sources', []) + [data['name']]
    destdir = '.'.join(os.path.dirname(data['name']).split('/'))
    data.pop('name')
    obj = Extension('%s.%s' % (destdir, name), **data)
    extensions.append(obj)


if __name__ == "__main__":
    min_versions = {
        'pysam': '0.8.1',
    }
    (setup_requires,
     install_requires) = check_dependency_versions(min_versions)
    if _have_setuptools:
        setuptools_kwargs['setup_requires'] = setup_requires
        setuptools_kwargs['install_requires'] = install_requires
        write_version_py()

    cwd = os.path.abspath(os.path.dirname(__file__))
    if not os.path.exists(os.path.join(cwd, 'PKG-INFO')) and not no_frills:
        # Generate Cython sources, unless building from source release
        generate_cython()

    setup(
        name=DISTNAME,
        maintainer=MAINTAINER,
        version=VERSION,
        ext_modules=extensions,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        license=LICENSE,
        url=URL,
        download_url=DOWNLOAD_URL,
        long_description=LONG_DESCRIPTION,
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
                  'pybedtools.scripts',
                  'pybedtools.test.data'],
        package_data={'pybedtools': ["test/data/*",
                                     "*.pyx",
                                     "*.pxi",
                                     "*.pxd",
                                     "*.cxx",
                                     "*.c",
                                     "*.cpp"]
                      },
        include_package_data=False,
        scripts=['pybedtools/scripts/venn_gchart.py',
                 'pybedtools/scripts/venn_mpl.py',
                 'pybedtools/scripts/annotate.py',
                 'pybedtools/scripts/peak_pie.py',
                 'pybedtools/scripts/intersection_matrix.py',
                 'pybedtools/scripts/intron_exon_reads.py',
                 'pybedtools/scripts/examples/pbt_plotting_example.py',
                 'pybedtools/scripts/pybedtools'],
        **setuptools_kwargs)
