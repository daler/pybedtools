#!/usr/bin/env python
"""
pybedtools utility scripts:

"""
from __future__ import print_function
import sys
import textwrap


def import_module(name):
    __import__("%s" % (name,), globals(), locals(), [], -1)
    return sys.modules[name]


def script_names(module):
    return sorted(
            [script for script in  module.__all__ if not script[:2] == "__"])


def main():
    m = import_module("pybedtools.scripts")
    scripts = script_names(m)
    mods = [import_module("pybedtools.scripts.%s" % s) for s in scripts]

    if (len(sys.argv) != 1 and not sys.argv[1] in scripts) \
       or len(sys.argv) == 1:

        print(__doc__.strip() + "\n")
        for name, mod in zip(scripts, mods):
            scriptname = " %-22s:" % name
            padding = ' ' * (len(scriptname) + 1)
            doclines = textwrap.wrap(textwrap.dedent(mod.main.__doc__), 50)
            print(scriptname, doclines[0])
            for line in doclines[1:]:
                print(padding, line)
            print()

    else:
        mname = sys.argv.pop(1)
        i = scripts.index(mname)
        module = mods[i]
        module.main()

if __name__ == "__main__":
    main()
