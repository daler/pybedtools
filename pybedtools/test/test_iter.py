import itertools
import yaml
import os
import pybedtools
# The functools.partial trick to get descriptions to be valid is from:
#
#   http://code.google.com/p/python-nose/issues/detail?id=244#c1
from functools import partial

this_dir = os.path.dirname(__file__)
config_fn = os.path.join(this_dir, 'test_cases.yaml')


def fix(x):
    """
    Replaces spaces with tabs, removes spurious newlines, and lstrip()s each
    line. Makes it really easy to create BED files on the fly for testing and
    checking.
    """
    s = ""
    for i in  x.splitlines():
        i = i.strip()
        if len(i) == 0:
            continue
        i = i.split()
        i = '\t'.join(i)+'\n'
        s += i
    return s


def parse_yaml(infile):
    x = yaml.load(open(infile).read())
    for test_case in x:
        method = test_case['method']
        send_kwargs = test_case['kwargs']
        expected = test_case['expected']
        yield method, send_kwargs, expected


def run(method, a, b, expected, **kwargs):
    result = getattr(a, method)(b, **kwargs)
    print 'Got:'
    print result
    print 'Expected:'
    print expected
    assert str(result) == expected

def test_a_b_methods():
    """
    Generator that yields tests, inserting different versions of `a` and `b` as
    needed
    """
    for method, send_kwargs, expected in parse_yaml(config_fn):

        if not ('a' in send_kwargs) and ('b' in send_kwargs):
            continue

        orig_a = pybedtools.example_bedtool(send_kwargs['a'])
        orig_b = pybedtools.example_bedtool(send_kwargs['b'])

        del send_kwargs['a']
        del send_kwargs['b']

        converter = {'file': lambda x: pybedtools.BedTool(x.fn),
                     'generator': lambda x: pybedtools.BedTool(i for i in x),
                     'stream': lambda x: pybedtools.BedTool(open(x.fn))
                    }
        done = []
        for kind_a, kind_b in itertools.permutations(('file', 'generator', 'stream'), 2):
                a = converter[kind_a](orig_a)
                b = converter[kind_b](orig_b)
                kind = 'a=%(kind_a)s, b=%(kind_b)s' % locals()
                f = partial(run, method, a, b, expected, **send_kwargs)
                f.description = '%(method)s, %(kind)s, %(send_kwargs)s' % locals()
                yield (f, )

