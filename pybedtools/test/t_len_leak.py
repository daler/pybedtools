import pybedtools
fn = pybedtools.example_filename('a.bed')
orig_fds = pybedtools.helpers.n_open_fds()

n1 = []
for i in range(10):
    x = pybedtools.BedTool(fn)
    n1.append(pybedtools.helpers.n_open_fds())
print max(n1) - orig_fds
# 2

n2 = []
for i in range(10):
    x = pybedtools.BedTool(fn)
    len(x)
    n2.append(pybedtools.helpers.n_open_fds())
print max(n2) - orig_fds
# 10

n3 = []
for i in range(10):
    x = pybedtools.BedTool(fn)
    sum(1 for _ in x)
    n3.append(pybedtools.helpers.n_open_fds())
print max(n3) - orig_fds
