"""
Plotting utilities for CABS with comprehensive type annotations.
"""

import matplotlib.pyplot as plt

plt.switch_backend("agg")
from itertools import chain
import json
from typing import Dict, List, Optional, Sequence, Tuple, Union

from matplotlib.axes import Axes
from matplotlib.patches import Patch
from matplotlib.ticker import EngFormatter, FuncFormatter, MaxNLocator
import numpy as np
import numpy.typing as npt

try:
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=["#666666", "#ff4000"])
except AttributeError:
    pass

# CABS SS code → (hex color, display label), matching webserver palette
# 1=Coil, 2=Helix, 3=Turn, 4=Sheet
_SS_COLORS: Dict[int, Tuple[str, str]] = {
    1: ("#AAAAAA", "Coil"),
    2: ("#AA00AA", "Helix"),
    3: ("#00AA00", "Turn"),
    4: ("#FF8800", "Sheet"),
}


def _draw_ss_band(ax: Axes, ss_vals: List[int], show_x_axis: bool = False) -> None:
    """Fill *ax* with a secondary-structure colour strip (one colour block per run)."""
    n = len(ss_vals)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, 1)
    if show_x_axis:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.get_yaxis().set_visible(False)
    else:
        ax.axis("off")
    i = 0
    while i < n:
        ss = ss_vals[i]
        j = i
        while j < n and ss_vals[j] == ss:
            j += 1
        color, _ = _SS_COLORS.get(ss, ("#AAAAAA", "Coil"))
        ax.broken_barh([(i - 0.5, j - i)], (0, 1), facecolors=color, edgecolors="none")
        i = j


def set_fixed_ar(plt_axes: Axes, ratio: float) -> None:
    """
    Set fixed aspect ratio for matplotlib subplot.

    Arguments:
        plt_axes: matplotlib subplot instance.
        ratio: aspect ratio to be set.
    """
    xvs = list(map(float, plt_axes.get_xlim()))
    yvs = list(map(float, plt_axes.get_ylim()))
    plt_axes.set_aspect(
        ratio * ((xvs[1] - xvs[0]) / (yvs[1] - yvs[0])), adjustable="box"
    )


def mk_discrete_plot(
    splot: Axes,
    xvals: Sequence[npt.NDArray[np.float64]],
    series: Sequence[npt.NDArray[np.float64]],
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    joined: bool = False,
) -> Axes:
    """
    Create discrete plot.

    Arguments:
        splot: plt.Axes instance; figure subplot to plot on.
        xvals: sequence of x-axis values.
        series: nested sequence of data series (floats or ints) of x-axis values len.
        xlim: x-axis limits.
        ylim: y-axis limits.
        joined: whether to join points with lines.
    """
    fmt = "o" + (":" if joined else " ")
    for sser, xsvals in zip(series, xvals):
        splot.plot(xsvals, sser, fmt, markersize=1, linewidth=1)
    if xlim:
        splot.set_xlim(xlim)
    if ylim:
        splot.set_ylim(ylim)
    return splot


def mk_histo(
    sfig: Union[Axes, List[Axes]],
    vls: Union[List[float], List[List[float]]],
    lbls: Union[List[str], List[List[str]]],
    ylim: Tuple[float, float] = (0.0, 1.0),
) -> None:
    """
    Create histogram plot.

    Arguments:
        sfig: matplotlib subplot or list of subplots.
        vls: list of subsequent histogram values (or list of corresponding number of lists if sfig is a list).
        lbls: list of subsequent labels (or list of lists if sfig is a list).
        ylim: y-axis limits.
    """
    if type(sfig) is not list:
        sfig = [sfig]
        vls = [vls]
        lbls = [lbls]
    for n, (vl, lb, sf) in enumerate(zip(vls, lbls, sfig)):
        xloc = [0.5 * i for i in range(len(lb))]
        sf.set_ylim(*ylim)
        sf.bar(xloc, vl, width=0.51)
        sf.set_xticks(xloc)
        sf.set_xticklabels(lb)
        sf.tick_params(labelsize=6)
    return sfig


def mk_histos_series(series, labels, fname, titles=None, fmt="svg", n_y_ticks=5):
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
    get_xloc = lambda x: [0.5 * i for i in range(len(x))]

    fig.set_figheight(len(series))

    for n, (vls, ticks) in enumerate(zip(series, labels)):
        xloc = get_xloc(ticks)
        sfigarr[n, 0].set_ylim((0, ylim))
        sfigarr[n, 0].bar(xloc, vls, width=0.51)

        sfigarr[n, 0].set_xticks(xloc)
        sfigarr[n, 0].set_xticklabels(ticks)
        sfigarr[n, 0].tick_params(labelsize=6)

    try:
        for k, title in titles.items():
            sfigarr[k, 0].set_title(title)
    except TypeError:
        pass

    plt.tight_layout()
    plt.savefig(fname + "." + fmt, format=fmt)
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
    with open(fname + ".csv", "w") as f:
        for vals in zip(*columns):
            f.write("\t".join([fmt % val for fmt, val in zip(fmts, vals)]))
            f.write("\n")


def plot_E_RMSD(trajectories, rmsds, labels, fname, fmt="svg", interaction=True):
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
        sets, labels = (0, 1), ("total", "interaction")
    else:
        sets, labels = (0,), ("total",)
    for ind, etp in zip(sets, labels):
        fig = plt.figure(figsize=(9, 12))
        grid = plt.GridSpec(2, 1)
        plot = plt.subplot(grid[0, 0])
        histo = plt.subplot(grid[1, 0])

        data = [
            [
                h.get_energy(mode=etp, number_of_peptides=traj.number_of_peptides)
                for h in traj.headers
            ]
            for traj in trajectories
        ]
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

        plot.set_xlabel("RMSD")
        plot.set_ylabel("%s energy" % etp.capitalize())
        plot.set_title("CABS %s energy vs. RMSD" % etp)
        set_fixed_ar(plot, 0.75)

        histo.set_xlim(xlim)
        histo.yaxis.set_major_formatter(
            EngFormatter(range(int(histo.get_ylim()[1] + 1)))
        )
        histo.set_xlabel("RMSD")
        histo.set_ylabel("Number of frames")
        set_fixed_ar(histo, 0.75)

        histo.legend(
            bbox_to_anchor=(0.0, -0.202, 1.0, -0.102),
            loc=3,
            ncol=len(labels),
            mode="expand",
            borderaxespad=0.0,
        )

        plt.savefig(fname + "_%s." % etp + fmt, format=fmt)
        plt.close()


def plot_RMSD_N(rmsds, fname, fmt="svg"):
    """Plots and saves to a file RMSD(Nframe) plot.

    Arguments:
    rmsds -- nested sequence of RMSDs.
    fname -- file name to be created.
    fmt -- format of file to be created; 'svg' byt default.
    See plt.savefig for more formats.
    """
    for n, rmsd_lst in enumerate(rmsds):
        tfname = fname + "_replica_%i" % n
        nfs = range(len(rmsd_lst))

        fig, sfig = plt.subplots(1)
        mk_discrete_plot(sfig, [nfs], [rmsd_lst], joined=True)
        sfig.set_ylabel("RMSD")
        sfig.set_xlabel("Frame index")
        plt.savefig(tfname + "." + fmt, format=fmt)
        plt.close(fig)
        drop_csv_file(tfname, (map(str, nfs), rmsd_lst), fmts=("%s", "%.3f"))


def graph_RMSF(trajectory, chains, fname, fmt="svg"):
    atoms = [i for i in trajectory.template.atoms if i.chid in chains]
    rmsf_vals = [trajectory.rmsf(chains)]
    lbls = [i.fmt() for i in atoms]
    ss_vals = [int(i.occ) for i in atoms]
    plot_RMSF_seq(rmsf_vals, lbls, fname + "_seq", fmt, ss_vals=ss_vals)
    drop_csv_file(fname, (lbls, tuple(chain(*rmsf_vals))), fmts=("%s", "%.3f"))
    rmsf_min = np.min(rmsf_vals[0])
    rmsf_max = np.max(rmsf_vals[0])
    rmsf_med = np.median(rmsf_vals[0])
    rmsf_med2 = 2 * rmsf_med - rmsf_min
    stats_dict = {
        "min": f"{rmsf_min:.3f}",
        "med": f"{rmsf_med:.3f}",
        "med2": f"{rmsf_med2:.3f}",
        "max": f"{rmsf_max:.3f}",
    }
    with open(fname + "_stats.json", "w") as f:
        json.dump(stats_dict, f)


def plot_RMSF_seq(series, labels, fname, fmt="svg", ss_vals=None):
    """
    Arguments:
    series   -- list of sequences of data to be plotted.
    labels   -- corresponding tick labels (strings like "A:1:MET").
    fname    -- file name to be created (without extension).
    fmt      -- image format; 'svg' by default.
    ss_vals  -- optional list of int CABS SS codes (1=Coil, 2=Helix, 3=Turn,
                4=Sheet) per residue. When supplied, coloured secondary-structure
                bands are drawn above and below the RMSF trace, matching the
                webserver style.
    """
    n = len(series[0])
    x = np.arange(n)

    unique_chains = set()
    for lbl in labels:
        if ":" in lbl:
            unique_chains.add(lbl.split(":")[0])
        else:
            for i, c in enumerate(lbl):
                if c.isdigit():
                    unique_chains.add(lbl[:i])
                    break
    show_chain = len(unique_chains) > 1

    if ss_vals is not None and len(ss_vals) == n:
        fig = plt.figure(figsize=(12, 5), constrained_layout=True)
        gs = fig.add_gridspec(
            2, 1, height_ratios=[1, 0.06], hspace=0.02
        )
        sfig = fig.add_subplot(gs[0])
        ax_bot = fig.add_subplot(gs[1], sharex=sfig)

        _draw_ss_band(ax_bot, ss_vals, show_x_axis=True)

        # sfig shares X-axis with ax_bot, we keep sfig xticklabels hidden and let ax_bot show them.
        plt.setp(sfig.get_xticklabels(), visible=False)

        sfig.set_title("RMSF with Secondary Structure")

        legend_handles = [
            Patch(facecolor=c, label=name)
            for _, (c, name) in sorted(_SS_COLORS.items())
        ]
        ax_bot.legend(
            handles=legend_handles,
            loc="upper left",
            bbox_to_anchor=(0.0, -1.5),
            ncol=4,
            fontsize=8,
            title="Secondary Structure",
            title_fontsize=8,
            frameon=False,
        )
        ax_x = ax_bot
    else:
        fig, sfig = plt.subplots(1, figsize=(10, 4))
        ax_x = sfig

    # RMSF trace
    sfig.plot(x, series[0], "-o", color="#1f77b4", markersize=2, linewidth=1)
    sfig.set_ylabel("RMSF")
    sfig.set_xlim(x[0] - 0.5, x[-1] + 0.5)

    # X-axis ticks: show just the residue number part of "A:1:MET"
    ax_x.xaxis.set_major_locator(MaxNLocator(25, integer=True))
    def format_label(v, p):
        idx = int(round(v))
        if 0 <= idx < len(labels):
            lbl = labels[idx]
            if ":" in lbl:
                parts = lbl.split(":")
                return f"{parts[0]}:{parts[1]}" if show_chain else parts[1]
            for i, c in enumerate(lbl):
                if c.isdigit():
                    return lbl if show_chain else lbl[i:]
            return lbl
        return ""

    ax_x.xaxis.set_major_formatter(FuncFormatter(format_label))

    for tick in ax_x.get_xticklabels():
        tick.set_rotation(0)
    ax_x.set_xlabel("Residue Number")

    if not (ss_vals is not None and len(ss_vals) == n):
        fig.tight_layout()
    plt.savefig(fname + "." + fmt, format=fmt, dpi=150, bbox_inches="tight")
    plt.close(fig)
