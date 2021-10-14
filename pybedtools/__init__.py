import os
import sys
import subprocess
import tempfile
import logging
from six.moves import urllib
from six.moves import copyreg
from .cbedtools import (
    Interval,
    IntervalFile,
    overlap,
    Attributes,
    MalformedBedLineError,
    IntervalIterator,
)
from . import contrib
from .helpers import (
    get_tempdir,
    set_tempdir,
    cleanup,
    find_tagged,
    set_bedtools_path,
    chromsizes,
    get_chromsizes_from_ucsc,
    chromsizes_to_file,
    create_interval_from_list,
)
from . import helpers
from .bedtool import BedTool
from . import genome_registry
from . import stats
from .version import __version__
from .filenames import data_dir, example_filename, list_example_files
from .bedtool import example_bedtool

from . import settings
from .logger import logger, ch

example_files = ["a.bed.", "b.bed", "test.fa", "a.bam"]


def debug_mode(x):
    """
    Enable debug mode.

    Use debug_mode(True) to show debug log events in the console and to save
    calling info in BedTool objects, and turn it off again with
    debug_mode(False).

    Note that `pybedtools.KEEP_TEMPFILES` will be set as well, so you will need
    to clean up the tempfile directory manually after using debug mode.
    """
    if x:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
        _DEBUG = True
        KEEP_TEMPFILES = True
        logger.info(
            "Debug mode enabled.  You may also want to set "
            "pybedtools.KEEP_TEMPFILES=True to prevent automatic deletion "
            "of files upon exit."
        )
    else:
        logger.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)
        _DEBUG = False
        KEEP_TEMPFILES = False
        logger.info("Debug mode disabled")


def check_for_bedtools(*args, **kwargs):
    """
    For backwards compatibility; please use helpers._check_for_bedtools()
    """
    return helpers._check_for_bedtools(*args, **kwargs)


# Allow Interval objects to be pickled -- required if you want to pass them
# across process boundaries
def interval_constructor(fields):
    return create_interval_from_list(list(fields))


def interval_reducer(interval):
    return interval_constructor, (tuple(interval.fields),)


copyreg.pickle(Interval, interval_reducer, interval_constructor)


def load_path_config(fn):
    """
    You can use a config file to specify installation paths of various programs
    used by pybedtools.  This can be useful for testing, or using different
    versions of programs.

    `fn` is a config file with the following format.  If an entry is blank,
    then assume it's already on the path. All items must be lowercase::

        [paths]
        bedtools=/tools/BEDTools/bin
        r=
        tabix=
        bgzip=

    You only need to specify paths you need to change, so this is a valid file
    that will only specify the path to use for R::

        [paths]
        r=/usr/bin/R-dev

    If `fn` is not a string, then assume it is a dictionary of (program,
    paths). This is used primarily for testing.
    """
    setters = dict(
        bedtools=helpers.set_bedtools_path,
        r=helpers.set_R_path,
        tabix=helpers.set_tabix_path,
        bgzip=helpers.set_bgzip_path,
    )

    if isinstance(fn, dict):
        for prog, setter in list(setters.items()):
            try:
                path = fn[prog]
                setter(path)
            except KeyError:
                pass

    if isinstance(fn, str):
        from six.moves import configparser

        c = configparser.SafeConfigParser()
        c.read(fn)
        if c.sections() != ["paths"]:
            raise ValueError(
                "Invalid path config -- must have " "only one section, [paths]."
            )
        for prog, setter in list(setters.items()):
            try:
                path = c.get("paths", prog)
                setter(path)

            except configparser.NoOptionError:
                pass
