import pybedtools

fn = pybedtools.example_filename('a.bed')

def show_open_fds(func):
    doc = func.__doc__
    print
    print doc
    print "." * len(doc)
    orig_fds = pybedtools.helpers.n_open_fds()
    print '\t', max(func(fn)) - orig_fds, 'file-based'


def func1(src):
    "create bedtool in loop"
    for i in range(10):
        x = pybedtools.BedTool(src)
        yield pybedtools.helpers.n_open_fds()


def func2(src):
    'create bedtool in loop and check length'
    for i in range(10):
        x = pybedtools.BedTool(src)
        len(x)
        yield pybedtools.helpers.n_open_fds()


def func3(src):
    'create bedtool outside of loop; check length inside'
    x = pybedtools.BedTool(src)
    for i in range(10):
        len(x)
        yield pybedtools.helpers.n_open_fds()


def func4(src):
    'create and len in loop; don\'t assign to var'
    for i in range(10):
        len(pybedtools.BedTool(src))
        yield pybedtools.helpers.n_open_fds()

if __name__ == "__main__":
    for k, v in sorted(locals().items()):
        if k.startswith('func'):
            show_open_fds(v)
