import time
import os
import pybedtools
from pybedtools.contrib import plotting
from matplotlib import pyplot as plt

colors = ['r', 'b', 'g']


def plot_a_b_tool(a, b, method, **kwargs):
    """
    Use for BEDTools programs that use -a and -b input arguments.  Filenames
    `a` and `b` are used for `method` of the BedTool class, and `kwargs` are
    sent to that method.

    The result is a plot of `a`, `b`, and the result, with the commandline
    argument as the plot title.
    """
    a = pybedtools.BedTool(a)
    b = pybedtools.BedTool(b)
    kwargs['b'] = b
    result = getattr(a, method)(**kwargs)

    fig = plt.figure(figsize=(8, 2))
    ax = fig.add_subplot(111)
    ybase = 0
    yheight = 1
    ylabels = []
    yticks = []
    for color, bt, label in zip(
            colors,
            [result, b, a],
            ['result', os.path.basename(b.fn), os.path.basename(a.fn)]):
        ylabels.append(label)
        track = plotting.Track(
                bt, visibility='squish', alpha=0.5, ybase=ybase, color=color)
        yticks.append(track.midpoint)
        ybase = track.ymax + 0.1
        ax.add_collection(track)
    ax.set_yticklabels(ylabels)
    ax.set_yticks(yticks)
    ax.set_title(' '.join([os.path.basename(i) for i in result._cmds]))
    ax.axis('tight')
    fig.subplots_adjust(top=0.8, bottom=0.15)
    return ax


def plot_i_tool(i, method, **kwargs):
    a = pybedtools.BedTool(i)
    result = getattr(a, method)(**kwargs)
    fig = plt.figure(figsize=(8, 2))
    ax = fig.add_subplot(111)
    res_track = plotting.Track(result, color='r', alpha=0.5, ybase=0, visibility='squish')
    a_track = plotting.Track(a, color='b', alpha=0.5, ybase=1.1, visibility='squish')
    ax.add_collection(res_track)
    ax.add_collection(a_track)
    ax.set_yticks([res_track.midpoint, a_track.midpoint])
    ax.set_yticklabels(['result', os.path.basename(a.fn)])
    ax.axis('tight')
    ax.set_title(' '.join([os.path.basename(k) for k in result._cmds]))
    fig.subplots_adjust(top=0.8, bottom=0.15)
    return ax




if __name__ == "__main__":
    a = pybedtools.example_filename('a.bed')
    b = pybedtools.example_filename('b.bed')

    plot_a_b_tool(a, b, 'intersect', u=True)
    plot_a_b_tool(a, b, 'intersect')
    plot_a_b_tool(a, b, 'subtract')

    plot_i_tool(a, 'merge')



    # Check performance -- should be <2s for 15k features
    t0 = time.time()
    fig = plt.figure(figsize=(8,2))
    ax = fig.add_subplot(111)
    big = pybedtools.example_bedtool('dm3-chr2L-5M.gff.gz')
    gene_track = plotting.Track(
            big.filter(lambda x: x[2] == 'gene'),
            color='k', visibility='squish', alpha=0.5, label='genes')
    exon_track = plotting.Track(
            big.filter(lambda x: x[2] == 'exon'),
            color='r', visibility='squish', alpha=0.5, ybase=gene_track.ymax,
            label='exons')
    ax.add_collection(gene_track)
    ax.add_collection(exon_track)

    ax.legend(loc='best')
    ax.axis('tight')

    print'%.2fs' % (time.time() - t0)


    plt.show()
