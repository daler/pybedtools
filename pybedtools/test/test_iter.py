from __future__ import print_function 
import difflib
import itertools
import yaml
import os
import gzip
import pybedtools

# The functools.partial trick to get descriptions to be valid is from:
#
#   http://code.google.com/p/python-nose/issues/detail?id=244#c1
from functools import partial

this_dir = os.path.dirname(__file__)
config_fn = os.path.join(this_dir, 'test_cases.yaml')

def gz(x):
    """
    Gzips a file to a tempfile, and returns a new BedTool using the gzipped
    version.
    """
    fin = open(x.fn)
    gzfn = pybedtools.BedTool._tmp()
    with gzip.open(gzfn, 'wb') as out_:
        with open(x.fn, 'rb') as in_:
            out_.writelines(in_)
    return pybedtools.BedTool(gzfn)

def fix(x):
    """
    Replaces spaces with tabs, removes spurious newlines, and lstrip()s each
    line. Makes it really easy to create BED files on the fly for testing and
    checking.
    """
    s = ""
    for i in x.splitlines():
        i = i.strip('\n\r')
        if len(i) == 0:
            continue

        # If the expected output contains tabs, then use those to split,
        # otherwise space.  This allows you to have expected output with blank
        # fields (e.g., "\t\t")
        if '\t' in i:
            i = i.split('\t')
        else:
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


def run(method, bedtool, expected, **kwargs):
    result = getattr(bedtool, method)(**kwargs)
    res = str(result)
    expected = fix(expected)
    try:
        assert res == expected

    except AssertionError:
        print(result.fn)
        print('Method call:')
        args = []
        for key, val in list(kwargs.items()):
            args.append(('%s=%s' % (key, val)).strip())

        args = ', '.join(args)
        print('BedTool.%(method)s(%(args)s)' % locals())
        print('Got:')
        print(res)
        print('Expected:')
        print(expected)
        print('Diff:')
        for i in difflib.unified_diff(res.splitlines(1), expected.splitlines(1)):
            print(i, end=' ')

        # Make tabs and newlines visible
        spec_res = res.replace('\t', '\\t').replace('\n', '\\n\n')
        spec_expected = expected.replace('\t', '\\t').replace('\n', '\\n\n')

        print('Showing special characters:')
        print('Got:')
        print(spec_res)
        print('Expected:')
        print(spec_expected)
        print('Diff:')
        for i in difflib.unified_diff(spec_res.splitlines(1), spec_expected.splitlines(1)):
            print(i, end=' ')
        raise


# List of methods that *only* take BAM as input
bam_methods = ('bam_to_bed',)

# List of supported BedTool construction from BAM files.  Currently only
# file-based.
supported_bam = ('filename', )

converter = {'filename': lambda x: pybedtools.BedTool(x.fn),
             'generator': lambda x: pybedtools.BedTool(i for i in x),
             'stream': lambda x: pybedtools.BedTool(open(x.fn)),
             'gzip': gz,
            }

def test_a_b_methods():
    """
    Generator that yields tests, inserting different versions of `a` and `b` as
    needed
    """
    for method, send_kwargs, expected in parse_yaml(config_fn):
        a_isbam = False
        b_isbam = False

        if 'abam' in send_kwargs:
            send_kwargs['abam'] = pybedtools.example_filename(send_kwargs['abam'])
            send_kwargs['a'] = send_kwargs['abam']
            a_isbam = True

        if not (('a' in send_kwargs) and ('b' in send_kwargs)):
            continue

        # If abam, makes a BedTool out of it anyway.
        orig_a = pybedtools.example_bedtool(send_kwargs['a'])
        orig_b = pybedtools.example_bedtool(send_kwargs['b'])

        del send_kwargs['a']
        del send_kwargs['b']

        if orig_a._isbam:
            a_isbam = True
        if orig_b._isbam:
            b_isbam = True

        for kind_a, kind_b in itertools.permutations(('filename', 'generator', 'stream', 'gzip'), 2):

            if a_isbam and (kind_a not in supported_bam):
                continue

            if b_isbam and (kind_b not in supported_bam):
                continue

            # Convert to file/generator/stream
            bedtool = converter[kind_a](orig_a)
            b = converter[kind_b](orig_b)

            kind = 'a=%(kind_a)s, b=%(kind_b)s abam=%(a_isbam)s bbam=%(b_isbam)s' % locals()

            send_kwargs['b'] = b

            f = partial(run, method, bedtool, expected, **send_kwargs)

            # Meaningful description
            f.description = '%(method)s, %(kind)s, %(send_kwargs)s' % locals()
            yield (f, )

def test_i_methods():
    """
    Generator that yields tests, inserting different versions of `i` as needed
    """
    for method, send_kwargs, expected in parse_yaml(config_fn):
        i_isbam = False
        if 'ibam' in send_kwargs:
            i_isbam = True
            send_kwargs['ibam'] = pybedtools.example_filename(send_kwargs['ibam'])
            send_kwargs['i'] = send_kwargs['ibam']

        if ('a' in send_kwargs) and ('b' in send_kwargs):
            continue

        if ('i' not in send_kwargs) and ('ibam' not in send_kwargs):
            continue

        if 'files' in send_kwargs:
            send_kwargs['files'] = [pybedtools.example_filename(i) for i in send_kwargs['files']]

        orig_i = pybedtools.example_bedtool(send_kwargs['i'])
        if orig_i._isbam:
            i_isbam = True

        del send_kwargs['i']

        done = []
        for kind_i in ('filename', 'generator', 'stream', 'gzip'):
            if i_isbam:
                if (kind_i not in supported_bam):
                    continue
            i = converter[kind_i](orig_i)
            kind = 'i=%(kind_i)s ibam=%(i_isbam)s' % locals()
            f = partial(run, method, i, expected, **send_kwargs)
            f.description = '%(method)s, %(kind)s, %(send_kwargs)s' % locals()
            yield (f, )

def test_bed_methods():
    """
    Generator that yields tests, inserting different versions of `bed` as needed
    """
    for method, send_kwargs, expected in parse_yaml(config_fn):
        ignore = ['a', 'b','abam','i']
        skip_test = False
        for i in ignore:
            if i in send_kwargs:
                skip_test = True
        if skip_test:
            continue
        if 'bed' not in send_kwargs:
            continue

        if 'files' in send_kwargs:
            send_kwargs['files'] = [pybedtools.example_filename(i) for i in send_kwargs['files']]

        if 'bams' in send_kwargs:
            send_kwargs['bams'] = [pybedtools.example_filename(i) for i in send_kwargs['bams']]

        if 'fi' in send_kwargs:
            send_kwargs['fi'] = pybedtools.example_filename(send_kwargs['fi'])

        orig_bed = pybedtools.example_bedtool(send_kwargs['bed'])

        del send_kwargs['bed']

        done = []
        for kind_bed in ('filename', 'generator', 'stream', 'gzip'):
            bed = converter[kind_bed](orig_bed)
            kind = 'i=%(kind_bed)s' % locals()
            f = partial(run, method, bed, expected, **send_kwargs)
            f.description = '%(method)s, %(kind)s, %(send_kwargs)s' % locals()
            yield (f, )

def teardown():
    pybedtools.cleanup(remove_all=True)
