"""
CABSDock module for contact map analysis of trajectories.

Created on 4 June 2017 by Tymoteusz hert Oleniecki.
"""

from functools import reduce
import operator

import matplotlib.pyplot as plt

plt.switch_backend("agg")
import matplotlib.ticker
import numpy as np

from CABS.analysis.plots import mk_histo
from CABS.utils.utils import _chunk_lst


class ContactMapFactory:
    def __init__(self, chains1, chains2, temp):
        """Builder for ContactMap.

        Arguments:
        chains1 -- list or str; chars for 1st chain(s).
        chains2 -- list or str; chars for 2nd chain(s).
        temp -- CABS.atom.Atoms instance containing both given chains.

        """
        chs = {}
        for n, i in enumerate(temp.atoms):
            chs.setdefault(i.chid, []).append(n)
        self.dims = (
            sum(map(len, [chs.get(ch1, []) for ch1 in chains1])),
            sum(map(len, [chs.get(ch2, []) for ch2 in chains2])),
        )
        self.inds1 = reduce(operator.add, [chs.get(i, []) for i in chains1])
        self.inds2 = reduce(operator.add, [chs.get(i, []) for i in chains2])
        self.ats1 = [temp.atoms[i] for i in self.inds1]
        self.ats2 = [temp.atoms[i] for i in self.inds2]
        res1 = {}
        res2 = {}
        for i, at in enumerate(self.ats1):
            res1.setdefault(at.fmt(), []).append(i)
        for i, at in enumerate(self.ats2):
            res2.setdefault(at.fmt(), []).append(i)
        self.res1 = sorted(
            [(k, (v[0], v[-1])) for k, v in res1.items()], key=lambda x: x[1][0]
        )
        self.res2 = sorted(
            [(k, (v[0], v[-1])) for k, v in res2.items()], key=lambda x: x[1][0]
        )

    def mk_cmap(self, traj, thr, frames=None, replicas=None):
        """Creates map of contacts between two given chains.

        Arguments:
        traj -- np.array of proper shape, i.e. Nreplicas x Nframes x Natoms x 3.
        thr -- float; threshold for side chain distance contact.
        frames -- tuple of ints; indexes of frames to be taken. All frames are taken by default.
        replicas -- tuple of replicas' indexes to be taken. All replicas are taken by default.

        Returns list of ContactMap for each replica in trajectory.
        """

        if replicas is None:
            replicas = range(traj.shape[0])
        if frames is None:
            fstf = 0
            frames = slice(1, None)
        else:
            frames = list(frames)
            fstf = frames.pop(0)
        resl = []
        for rep in traj[replicas,]:
            cmtx = self.mk_cmtx(self.mk_dmtx(rep[fstf]), thr)
            nframes = 1
            for fra in rep[frames,]:
                ncmtx = self.mk_cmtx(self.mk_dmtx(fra), thr)
                cmtx += ncmtx
                nframes += 1
            resl.append(
                ContactMap(
                    cmtx, [i[0] for i in self.res1], [i[0] for i in self.res2], nframes
                )
            )
        return resl

    def mk_cmtx(self, mtx, thr):
        """Returns boolean np.array of contacts from given mtx of distances.

        Arguments:
        mtx -- np.array of distances between atoms.
        thr -- thresohld below which atoms are in contact.
        """
        mtx = np.clip(np.sign(-mtx + thr), 0, 1)
        if len(self.ats1) == len(self.res1) and len(self.ats2) == len(self.res2):
            return mtx
        nmtx = np.zeros((len(self.res1), len(self.res2)))
        for i, (r1, inds1) in enumerate(self.res1):
            for j, (r2, inds2) in enumerate(self.res2):
                mtx_slice = mtx[inds1[0] : inds1[1], inds2[0] : inds2[1]]
                if mtx_slice.size > 0:
                    nmtx[i, j] = np.max(mtx_slice)
                else:
                    nmtx[i, j] = 0.0
        return nmtx

    def mk_dmtx(self, vec):
        """Returns 2D np.array of distances between atoms from vector of coordinates.

        Arguments:
        vec -- slice of trajectory.

        Calculates distances between atoms from given chains.
        """
        m1 = vec[self.inds1,].reshape(-1, 1, 3)[:, self.dims[1] * (0,)]
        m2 = vec[self.inds2,].reshape(1, -1, 3)[self.dims[0] * (0,), :]
        mtx = m1 - m2
        return (mtx * mtx).sum(axis=2) ** 0.5


class ContactMap:
    def __init__(self, mtx, nms1, nms2, n):
        """Contact map init.

        Arguments:
        mtx -- 2D np.array of distances between (pseudo)atoms.
        atoms1, atoms2 -- CABS.atom.Atoms instance; template for cmap.
        n -- number of frames.
        """
        self.cmtx = mtx
        self.s1 = nms1
        self.s2 = nms2
        self.n = n

    def zero_diagonal(self):
        np.fill_diagonal(self.cmtx, 0)

    def save_fig(
        self,
        fname,
        fmt="svg",
        norm_n=False,
        break_long_x=0,  # We default to 0 for a single continuous map
        colors_lst=None,
    ):
        """Saves cmap as matrix plot resembling the webserver style."""
        if self.cmtx.shape[0] == 0 or self.cmtx.shape[1] == 0 or len(self.s1) == 0 or len(self.s2) == 0:
            fig = plt.figure()
            plt.savefig(fname + "." + fmt, format=fmt)
            plt.close(fig)
            return

        # Always use the premium blue gradient for contact maps
        colors_lst = ["#ffffff", "#e0eafc", "#bfe3fd", "#93c5fd", "#60a5fa", "#3b82f6", "#2563eb", "#1d4ed8", "#1e3a8a", "#0b2c7a"]

        fig, sfig = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)

        if norm_n and self.n > 0:
            plot_mtx = self.cmtx.T / float(self.n)
            vmax = 1.0
        else:
            plot_mtx = self.cmtx.T
            vmax = np.max(self.cmtx) if self.cmtx.size > 0 else 1.0
            if vmax < 1:
                vmax = 1.0

        colors = matplotlib.colors.LinearSegmentedColormap.from_list(
            "blue_gradient",
            ["#" + color if "#" not in color else color for color in colors_lst],
        )

        # Plot matrix with origin='lower' to get diagonal from bottom-left to top-right
        im = sfig.imshow(
            plot_mtx,
            cmap=colors,
            vmin=0.0,
            vmax=vmax,
            origin="lower",
            aspect="equal",
        )

        # Extract unique chains in s1 and s2
        def get_unique_chains(labels_list):
            if not labels_list:
                return []
            chains = []
            for lbl in labels_list:
                ch = lbl.split(":")[0] if ":" in lbl else (lbl[0] if lbl else "A")
                if ch not in chains:
                    chains.append(ch)
            return chains

        chains1 = get_unique_chains(self.s1)
        chains2 = get_unique_chains(self.s2)
        show_chain_x = len(chains1) > 1
        show_chain_y = len(chains2) > 1

        # Configure X and Y ticks to show clean residue numbers (with chain ID if multiple chains exist)
        def clean_ticks(labels_list, show_chain_id, n_ticks=6):
            if not labels_list:
                return [], []
            inds = np.linspace(0, len(labels_list) - 1, n_ticks).astype(int)
            short_labels = []
            for idx in inds:
                lbl = labels_list[idx]
                if ":" in lbl:
                    parts = lbl.split(":")
                    short_labels.append(f"{parts[0]}:{parts[1]}" if show_chain_id else parts[1])
                else:
                    for i, c in enumerate(lbl):
                        if c.isdigit():
                            short_labels.append(lbl if show_chain_id else lbl[i:])
                            break
                    else:
                        short_labels.append(lbl)
            return inds, short_labels

        x_inds, x_lbls = clean_ticks(self.s1, show_chain_x)
        y_inds, y_lbls = clean_ticks(self.s2, show_chain_y)

        sfig.set_xticks(x_inds)
        sfig.set_xticklabels(x_lbls)
        sfig.set_yticks(y_inds)
        sfig.set_yticklabels(y_lbls)

        sfig.tick_params(
            bottom=True,
            top=False,
            labelbottom=True,
            labeltop=False,
            labelsize=9,
        )

        # Dynamic axis labels supporting single or multiple chains
        ch1_str = ", ".join(chains1)
        ch2_str = ", ".join(chains2)
        xlabel_prefix = "Chains" if len(chains1) > 1 else "Chain"
        ylabel_prefix = "Chains" if len(chains2) > 1 else "Chain"
        sfig.set_xlabel(f"Residues in {xlabel_prefix} {ch1_str}", fontsize=10, labelpad=8)
        sfig.set_ylabel(f"Residues in {ylabel_prefix} {ch2_str}", fontsize=10, labelpad=8)

        # Set title
        sfig.set_title(f"Contact Map: {ch1_str} vs. {ch2_str}" if ch1_str != ch2_str else f"Contact Map: {ch1_str} vs. {ch1_str}", fontsize=11, pad=10)

        # Remove top/right spines
        sfig.spines["top"].set_visible(False)
        sfig.spines["right"].set_visible(False)

        # Add clean vertical colorbar on the right
        cbar = fig.colorbar(im, ax=sfig, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=9)
        cbar.ax.set_title("Frequency" if self.n > 1 else "Contact", fontsize=9, pad=8)

        plt.savefig(fname + "." + fmt, format=fmt, dpi=150)
        plt.close(fig)

    def save_histo(self, fname, all_inds_stc2=True, fmt="svg"):
        """
        Saves histogram of contact counts for each atom.

        Arguments:
        fname -- str; name of file to be created.
        titles -- dict int: str; keys are indexes of histos, values are title to be set.
        all_inds_stc2 -- bool; True by default
        """
        if self.cmtx.shape[0] == 0 or self.cmtx.shape[1] == 0 or len(self.s1) == 0 or len(self.s2) == 0:
            fig = plt.figure()
            plt.savefig(fname + "." + fmt, format=fmt)
            plt.close(fig)
            return

        max_bars = 15

        trg_vls_all = self.cmtx.sum(axis=1) / float(self.n)
        pep_vls = self.cmtx.sum(axis=0) / float(self.n)

        inds1 = np.nonzero(trg_vls_all)[0]
        if all_inds_stc2:
            inds2 = range(self.cmtx.shape[1])
        else:
            inds2 = np.nonzero(pep_vls)[0]

        trg_vls = list(np.array(trg_vls_all)[inds1,])
        trg_lbls = list(np.array(self.s1)[inds1,])

        pep_lbls = [self.s2[i] for i in inds2]

        max_y = max([0.05] + trg_vls_all)

        chunks = _chunk_lst(trg_vls, max_bars, extend_last=0.0)
        grid = plt.GridSpec(2 + len(chunks), 1)
        size = (10, 3 * len(chunks))
        fig = plt.figure(figsize=size)

        peptH = mk_histo(
            plt.subplot(grid[0, 0]), pep_vls, pep_lbls, ylim=(0, max([0.05] + pep_vls))
        )[0]
        sbplts = [plt.subplot(grid[i, 0]) for i in range(1, len(chunks) + 1)]
        targBH = mk_histo(
            sbplts,
            chunks,
            _chunk_lst(trg_lbls, max_bars, extend_last=""),
            ylim=(0, max_y),
        )
        targAH = mk_histo(
            plt.subplot(grid[-1, 0]), trg_vls_all, self.s1, ylim=(0, max_y)
        )[0]

        peptH.set_title("Histogram of peptide contacts")
        targBH[0].set_title("Histogram of receptor contacts - detailed analysis")
        targAH.set_title("Histogram of receptor contacts - summary analysis")

        for sfig in targBH + [peptH, targAH]:
            sfig.set_ylabel("Contact frequency")
            sfig.set_xlabel("Residue id")

        xloc = targAH.get_xticks()
        if len(xloc) > max_bars:
            inds = np.linspace(0, len(xloc) - 1, max_bars).astype(int)
            targAH.set_xticks(xloc[inds,])
            targAH.set_xticklabels(np.array(self.s1)[inds,])

        grid.tight_layout(fig)
        plt.savefig(fname + "." + fmt, format=fmt)
        plt.close()

    def save_txt(self, stream):
        """Saves contact list in CSV format.

        Argument:
        stream -- file-like object; stream to which text will be passed.
        """
        inds1, inds2 = np.nonzero(self.cmtx)
        stream.write(f"# n={self.n}\n")
        for m1, m2, (c1, c2) in zip(
            [self.s1[i] for i in inds1], [self.s2[i] for i in inds2], zip(inds1, inds2)
        ):
            stream.write(f"{m1}\t{m2}\t{self.cmtx[c1, c2]:.3f}\n")

    def save_all(
        self,
        fname,
        norm_n=False,
        break_long_x=50,
        colors=["#ffffff", "#f2d600", "#4b8f24", "#666666", "#e80915", "#000000"],
    ):
        """Creates txt and png of given name."""
        with open(fname + ".txt", "w") as f:
            self.save_txt(f)
        self.save_fig(
            fname, norm_n=norm_n, break_long_x=break_long_x, colors_lst=colors
        )

    def __add__(self, other):
        """Sum contact maps matrices and return new object.

        Raises ValueError for cmaps of different particles.
        """
        if self.s1 != other.s1 or self.s2 != other.s2:
            raise ValueError("Cannot sum different particles' contact maps.")
        return ContactMap(self.cmtx + other.cmtx, self.s1, self.s2, self.n + other.n)
