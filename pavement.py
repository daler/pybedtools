# -*- Import: -*-
from paver.easy import *
from paver.setuputils import setup
from setuptools import find_packages

try:
    # Optional tasks, only needed for development
    # -*- Optional import: -*-
    from github.tools.task import *
    import paver.doctools
    import paver.virtual
    import paver.misctasks
    ALL_TASKS_LOADED = True
except ImportError, e:
    info("some tasks could not not be imported.")
    debug(str(e))
    ALL_TASKS_LOADED = False



version = '0.1'

classifiers = [
    # Get more strings from http://www.python.org/pypi?%3Aaction=list_classifiers
    "Development Status :: 1 - Planning",
    ]

install_requires = [
    # -*- Install requires: -*-
    'setuptools',
    ]

entry_points="""
    # -*- Entry points: -*-
    """

# compatible with distutils of python 2.3+ or later
setup(
    name='pybedtools',
    version=version,
    description="Wrapper for Aaron Quinlan's BEDtools, plus other useful methods for working with BED files",
    long_description=open('README.rst', 'r').read(),
    classifiers=classifiers,
    keywords='bioinformatics genome',
    author='Ryan Dale',
    author_email='dalerr@niddk.nih.gov',
    url='',
    license='BSD',
    packages = find_packages(exclude=['bootstrap', 'pavement',]),
    include_package_data=True,
    test_suite='nose.collector',
    zip_safe=False,
    install_requires=install_requires,
    entry_points=entry_points,
    )

options(
    # -*- Paver options: -*-
    minilib=Bunch(
        extra_files=[
            # -*- Minilib extra files: -*-
            ]
        ),
    sphinx=Bunch(
        docroot='docs',
        builddir="_build",
        sourcedir=""
        ),
    virtualenv=Bunch(
        packages_to_install=[
            # -*- Virtualenv packages to install: -*-
            'github-tools',
            "nose",
            "Sphinx>=0.6b1",
            "pkginfo", 
            "virtualenv"],
        dest_dir='./virtual-env/',
        install_paver=True,
        script_name='bootstrap.py',
        paver_command_line=None
        ),
    )

options.setup.package_data=paver.setuputils.find_package_data(
    'pybedtools', package='pybedtools', only_in_packages=False)

if ALL_TASKS_LOADED:
    @task
    @needs('generate_setup', 'minilib', 'setuptools.command.sdist')
    def sdist():
        """Overrides sdist to make sure that our setup.py is generated."""