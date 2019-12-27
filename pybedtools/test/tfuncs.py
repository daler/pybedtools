import pytest
import pybedtools
import os

testdir = os.path.dirname(__file__)
test_tempdir = os.path.join(os.path.abspath(testdir), "tmp")
unwriteable = os.path.join(os.path.abspath(testdir), "unwriteable")


def setup_module():
    if not os.path.exists(test_tempdir):
        os.system("mkdir -p %s" % test_tempdir)
    pybedtools.set_tempdir(test_tempdir)


def teardown_module():
    if os.path.exists(test_tempdir):
        os.system("rm -r %s" % test_tempdir)
    pybedtools.cleanup()
