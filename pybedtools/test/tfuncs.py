import pybedtools
import os

testdir = os.path.dirname(__file__)
unwriteable = os.path.join(os.path.abspath(testdir), "unwriteable")


def teardown_module():
    pybedtools.cleanup()
