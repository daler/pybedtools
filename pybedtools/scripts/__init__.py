import os.path as op
import glob

patt = op.join(op.abspath(op.dirname(__file__)), "*.py")
__all__ = [op.splitext(op.basename(p))[0] for p in glob.glob(patt)]
