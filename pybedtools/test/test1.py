import pybedtools

pybedtools.set_tempdir('.')
a = pybedtools.bedtool('a.bed')
a.intersect(a)
pybedtools.cleanup(verbose=True)
