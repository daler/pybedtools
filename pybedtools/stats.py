import multiprocessing
from . import helpers
import pybedtools


def random_jaccard(x, y, genome_fn, shuffle_kwargs, jaccard_kwargs):
    z = x.shuffle(g=genome_fn, **shuffle_kwargs).sort()
    result = z.jaccard(y, **jaccard_kwargs)
    helpers.close_or_delete(z)
    return result


def random_intersection(x, y, genome_fn, shuffle_kwargs, intersect_kwargs):
    z = x.shuffle(g=genome_fn, **shuffle_kwargs)
    zz = z.intersect(y, stream=True, **intersect_kwargs)
    result = len(zz)
    helpers.close_or_delete(z, zz)
    return result


def random_intersection_bp(x, y, genome_fn, shuffle_kwargs, intersect_kwargs):
    z = x.shuffle(g=genome_fn, **shuffle_kwargs)
    zz = z.intersect(y, stream=True, **intersect_kwargs)
    result = sum(len(i) for i in zz)
    helpers.close_or_delete(z, zz)
    return result
