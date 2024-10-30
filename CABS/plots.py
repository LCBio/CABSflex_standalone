import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import EngFormatter

import numpy as np
from itertools import chain

try:
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['#666666', '#ff4000'])
except AttributeError:
    pass


def set_fixed_ar(plt, ratio):
    """
    Arguments:
    plt -- matplotlib subplot instance.
    ratio -- aspect ratio to be set.
    """
    xvs = list(map(float, plt.get_xlim()))
    yvs = list(map(float, plt.get_ylim()))
    plt.set_aspect(ratio * ((xvs[1] - xvs[0]) / (yvs[1] - yvs[0])), adjustable='box')

def mk_discrete_plot(splot, xvals, series, xlim=None, ylim=None, joined=False):
    """
    Arguments:
    splot -- plt.Axes instance; figure subplot to plot on.
    xvals -- sequence of x-axis values.
    series -- nested sequence of data series (floats or ints) of x-axis values len.
    """
    fmt = "o" + (":" if joined else " ")
    for sser, xsvals in zip(series, xvals):
        splot.plot(xsvals, sser, fmt, markersize=1, linewidth=1)
    if xlim:
        splot.set_xlim(xlim)
    if ylim:
        splot.set_ylim(ylim)
    return splot

def mk_histo(sfig, vls, lbls, ylim=(0., 1.)):
    """
    Arguments:
    sfig -- matplotlib subplot or list of subplots.
    vls -- list of subsequent histogram values (or list of corresponding number of lists if sfig is a list).
    lbls -- list of subsequent labels (or list of lists if sfig is a list).
    """
    if type(sfig) is not list:
        sfig = [sfig]
        vls = [vls]
        lbls = [lbls]
    for n, (vl, lb, sf) in enumerate(zip(vls, lbls, sfig)):
        xloc = [.5 * i for i in range(len(lb))]
        sf.set_ylim(*ylim)
        sf.bar(xloc, vl, width=.51)
        sf.set_xticks(xloc)
        sf.set_xticklabels(lb)
        sf.tick_params(labelsize=6)
    return sfig

def mk_histos_series(series, labels, fname, titles=None, fmt='svg', n_y_ticks=5):
    """
    Arguments:
    series -- list of sequences of data to be plotted.
    labels -- corresponding ticks labels.
    fname -- file name to be created.
    titles -- dict int: str; keys are indexes of histos, values are title to be set.
    fmt -- format of file to be created; 'svg' byt default.
    n_y_ticks -- int; maximal number of tickes on y axis.
    See plt.savefig for more formats.
    """
    fig, sfigarr = plt.subplots(len(series), squeeze=False)

    try:
        ylim = max(chain((n_y_ticks,), *series)) + 1
    except ValueError:  # for empty series
        ylim = n_y_ticks + 1
    get_xloc = lambda x: [.5 * i for i in range(len(x))]

    fig.set_figheight(len(series))

    for n, (vls, ticks) in enumerate(zip(series, labels)):
        xloc = get_xloc(ticks)
        sfigarr[n, 0].set_ylim((0, ylim))
        sfigarr[n, 0].bar(xloc, vls, width=.51)
        #~ sfigarr[n, 0].yaxis.set_major_locator(MaxNLocator(n_y_ticks))
        #~ sfigarr[n, 0].yaxis.set_major_formatter(FuncFormatter(lambda x, pos: "%i" % x))
        sfigarr[n, 0].set_xticks(xloc)
        sfigarr[n, 0].set_xticklabels(ticks)
        sfigarr[n, 0].tick_params(labelsize=6)

    try:
        for k, title in titles.items():
            sfigarr[k, 0].set_title(title)
    except TypeError:
        pass

    plt.tight_layout()
    plt.savefig(fname + '.' + fmt, format=fmt)
    plt.close(fig)

def drop_csv_file(fname, columns, fmts="%s"):
    """
    Creates *fname* csv file and writes given columns to this file.

    Arguments:
    fname -- name of file to be created.
    columns -- sequences of data sequences. Lists will be truncated to len of the shortest one.
    fmts -- str or sequence of str. C-style string formats to be used (e.g. "%s", "%.3f", ...).
    If only one string is given -- same format is used for all data.
    Otheriwse subsequent fmts are used for corresponding solumns.
    """
    if type(fmts) is str:
        fmts = [fmts for dummy in columns]
    with open(fname + '.csv', 'w') as f:
        for vals in zip(*columns):
            f.write("\t".join([fmt % val for fmt, val in zip(fmts, vals)]))
            f.write('\n')

def plot_E_RMSD(trajectories, rmsds, labels, fname, fmt='svg', interaction=True):
    """
    Creates energy(RMSD) plots.

    Arguments:
    trajectories -- sequence of trajectories to be used.
    rmsds -- nested sequence of RMDSs.
    fname -- file name to be created.
    fmt -- format of file to be created; 'svg' byt default.
    interaction -- bool; if True plots both, total and interaction plots, otherwise plots only total energy.
    See plt.savefig for more formats.

    Plots figure with three subplots: total and internal energy vs RMSD
    and histogram of RMSDs. All three will be plotted for given data series,
    so nested arrays or lists are expected. Plots will be written in given format
    to file of given name.
    """
    max_data = 5
    if interaction:
        sets, labels = (0, 1), ('total','interaction')
    else:
        sets, labels = (0,), ('total',)
    for ind, etp in zip(sets, labels):
        fig = plt.figure(figsize=(9, 12))
        grid = plt.GridSpec(2, 1)
        plot = plt.subplot(grid[0, 0])
        histo = plt.subplot(grid[1, 0])

        data = [[h.get_energy(mode=etp, number_of_peptides=traj.number_of_peptides) for h in traj.headers] for traj in trajectories]
        xlim = (0, max(chain((max_data,), *rmsds)))
        ylim = (min(chain((-max_data,), *data)), max(chain(*data)))
        mk_discrete_plot(plot, rmsds, data, xlim, ylim)
        drop_csv_file(fname + "_%s" % etp, (rmsds[0], data[0]), fmts="%.3f")

        for traj, rmsd_list, lab in zip(trajectories, rmsds, labels):
            n_bins = np.arange(0, max(max_data, np.max(rmsd_list)), 1)
            histo.hist(rmsd_list, n_bins, label=lab)

        for sfig in (plot, histo):
            sfig.xaxis.set_major_locator(MaxNLocator(10))
            sfig.xaxis.set_minor_locator(MaxNLocator(20))

        plot.set_xlabel('RMSD')
        plot.set_ylabel('%s energy' % etp.capitalize())
        plot.set_title('CABS %s energy vs. RMSD' % etp)
        set_fixed_ar(plot, .75)

        histo.set_xlim(xlim)
        histo.yaxis.set_major_formatter(EngFormatter(range(int(histo.get_ylim()[1] + 1))))
        histo.set_xlabel('RMSD')
        histo.set_ylabel('Number of frames')
        set_fixed_ar(histo, .75)

        histo.legend(bbox_to_anchor=(0., -.202, 1., -.102), loc=3, ncol=len(labels), mode="expand", borderaxespad=0.)

        plt.savefig(fname + '_%s.' % etp + fmt, format=fmt)
        plt.close()

def plot_RMSD_N(rmsds, fname, fmt='svg'):
    """Plots and saves to a file RMSD(Nframe) plot.

    Arguments:
    rmsds -- nested sequence of RMSDs.
    fname -- file name to be created.
    fmt -- format of file to be created; 'svg' byt default.
    See plt.savefig for more formats.
    """
    for n, rmsd_lst in enumerate(rmsds):
        tfname = fname + '_replica_%i' % n
        nfs = range(len(rmsd_lst))

        fig, sfig = plt.subplots(1)
        mk_discrete_plot(sfig, [nfs], [rmsd_lst], joined=True)
        sfig.set_ylabel('RMSD')
        sfig.set_xlabel('Frame index')
        plt.savefig(tfname + '.' + fmt, format=fmt)
        plt.close(fig)
        drop_csv_file(tfname, (map(str, nfs), rmsd_lst), fmts=("%s", "%.3f"))

def graph_RMSF(trajectory, chains, fname, fmt='svg'):
    #~ if hist: # old version producting histo instead of plot; self is class from CABS.job
        #~ rmsf_vals = _chunk_lst(trajectory.rmsf(self.initial_complex.receptor_chains), 15, 0)
        #~ lbls = _chunk_lst([i.fmt() for i in trajectory.template.atoms if i.chid in self.initial_complex.receptor_chains], 15, "")
        #~ mk_histos_series(rmsf_vals, lbls, fname + '_hist')
    rmsf_vals = [trajectory.rmsf(chains)]
    lbls = [i.fmt() for i in trajectory.template.atoms if i.chid in chains]
    plot_RMSF_seq(rmsf_vals, lbls, fname + '_seq',fmt)
    drop_csv_file(fname, (lbls, tuple(chain(*rmsf_vals))), fmts=("%s", "%.3f"))

def plot_RMSF_seq(series, labels, fname, fmt='svg'):
    """
    Arguments:
    series -- list of sequences of data to be plotted.
    labels -- corresponding ticks labels.
    fname -- file name to be created.
    fmt -- format of file to be created; 'svg' byt default.
    See plt.savefig for more formats.
    """
    fig, sfig = plt.subplots(1)
    xvals = [range(len(series[0]))]
    mk_discrete_plot(sfig, xvals, series, xlim=(min(xvals[0]), max(xvals[0])), joined=True)
    sfig.set_ylabel('RMSF')
    sfig.set_xlabel('Residue index')
    sfig.xaxis.set_major_locator(MaxNLocator(25, integer=True))
    sfig.xaxis.set_major_formatter(FuncFormatter(lambda x, p: labels[min(int(x), len(labels) - 1)]))
    for tick in sfig.get_xticklabels():
        tick.set_rotation(90)
    fig.tight_layout()
    plt.savefig(fname + '.' + fmt, format=fmt)
    plt.close(fig)
