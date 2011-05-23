
import pybedtools
# The functools.partial trick to get descriptions to be valid is from:
#
#   http://code.google.com/p/python-nose/issues/detail?id=244#c1
from functools import partial



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


# Actual test looks like this.
def run(method, a, b, expected, **kwargs):
    result = getattr(a, method)(b, **kwargs)
    print result
    assert str(result) == expected

# Generator that yields tests, inserting `a` and `b` as needed
def test():

    a = pybedtools.example_bedtool('a.bed')
    orig_b = pybedtools.example_bedtool('b.bed')

    # Contains (kwargs, expected) tuples.
    # TODO: refactor this into YAML
    expecteds = (
                   (dict(s=True), fix("""
                             chr1	155	200	feature3	0	-
                             chr1	900	901	feature4	0	+
                             """)),

                   (dict(s=False), fix("""
                             chr1	155	200	feature2	0	+
                             chr1	155	200	feature3	0	-
                             chr1	900	901	feature4	0	+
                             """)),
                )

    for send_kwargs, expected in expecteds:
        b = orig_b
        kind = 'file'
        method = 'intersect'
        f = partial(run, method, a, b, expected, **send_kwargs)
        f.description = '%(method)s, %(kind)s, %(send_kwargs)s' % locals()
        yield (f, )

        b = (i for i in orig_b)
        kind = 'generator'
        method = 'intersect'
        f = partial(run, method, a, b, expected, **send_kwargs)
        f.description = '%(method)s, %(kind)s, %(send_kwargs)s' % locals()
        yield (f, )

        b = open(orig_b.fn)
        kind = 'stream'
        method = 'intersect'
        f = partial(run, method, a, b, expected, **send_kwargs)
        f.description = '%(method)s, %(kind)s, %(send_kwargs)s' % locals()
        yield (f, )
