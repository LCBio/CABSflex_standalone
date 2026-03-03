from abc import ABCMeta, abstractmethod
from functools import reduce
import operator
from os import remove
from subprocess import check_output
from tempfile import NamedTemporaryFile

import numpy as np

from CABS.config_loader import get_blosum62_indices, get_blosum62_matrix
from CABS.structures.atom import Atoms
from CABS.utils.utils import aa_to_short

# Load BLOSUM62 data from configuration
BLOSUM62 = np.array(get_blosum62_matrix())
B62h = get_blosum62_indices()


def raise_aerror_on(*errors):
    """Method wrapper that takes list of errors raised by method to be changed into AlignError."""

    def wrapper(method):
        def wrapped_mth(*args, **kwargs):
            try:
                return method(*args, **kwargs)
            except errors:
                raise AlignError

        return wrapped_mth

    return wrapper


def fmt_csv(atm):
    return f"{atm.chid}:{atm.resnum}{atm.icode.strip()}:{aa_to_short(atm.resname)}"


def save_csv(fname, stcs, aligned_mers):
    """Creates csv alignment file.

    Argument:
    fname -- str; name of file to be created.
    stcs -- tuple of structure names.
    aligned_mers -- sequence of tuples containing aligned CABS.atoms.Atom instances.
    """
    with open(fname, "w") as f:
        f.write("\t".join(stcs) + "\n")
        for mrs in aligned_mers:
            f.write("\t".join(map(fmt_csv, mrs)) + "\n")


def save_fasta(fname, stcs_names, stcs, aligned_mers):
    """Saves fasta file with alignment.

    Arguments:
    fname -- str; name of the file to be created.
    stcs_names -- tuple of strings; names of structures to be saved.
    stcs -- CABS.Atoms instances with atoms attribute. Full aligned structures.
    aligned_mers -- sequence of tuples containing only aligned CABS.atoms.Atom instances.
    """
    it1, it2 = map(iter, [i.atoms for i in stcs])
    txt1 = ""
    txt2 = ""
    for mrs in aligned_mers:
        while True:
            mer = next(it1)
            txt1 += aa_to_short(mer.resname)
            if mer not in mrs:
                txt2 += "-"
                continue
            break
        while True:
            mer = next(it2)
            txt2 += aa_to_short(mer.resname)
            if mer not in mrs:
                txt1 += "-"
                continue
            break
    with open(fname, "w") as f:
        for name, seq in zip(stcs_names, (txt1, txt2)):
            f.write(f">{name}\n{seq}\n")


def load_csv(fname, *stcs):
    """Loads pairwise alignment in csv format.

    Arguments:
    fname -- filelike object; stream containing alignment in csv format.
    stcs -- CABS.atom.Atoms instances.
    """
    dcts = [
        {fmt_csv(atm): atm for atm in stc.select("name CA and not HETERO")}
        for stc in stcs
    ]

    res = []
    for line in fname.read().split("\n")[1:]:
        if line.strip() == "":
            continue
        ms = line.split("\t")
        res.append([dct[mer] for mer, dct in zip(ms, dcts)])
    return res


def align_to(ref_stc, ref_chs, trg_stc, trg_chs, align_mth="SW", kwargs={}):
    """Calculates alignment of template to given reference structure.

    Arguments:
    ref_stc -- CABS.PDBlib.PDB instance of reference structure.
    ref_chs -- str; chain id(s) of reference selection.
    trg_stc -- CABS.PDBlib.PDB instance of target to be aligned with reference.
    trg_chs -- str; chain id(s) of trajectory structure selection.
    align_mth -- str; name of aligning method to be used. See CABS.align documentation for more information.
    kwargs -- as above, but used when aligning target protein.

    One needs to specify chains to be taken into account during alignment calculation.

    Returns two structures: reference and template -- both cropped to aligned parts only, and alignment as list of tuples.
    """
    mth = AbstractAlignMethod.get_subclass_dict()[align_mth]
    # aligning target
    mtch_mtx = np.zeros((len(ref_chs), len(trg_chs)), dtype=int)
    algs = {}
    key = 1
    # rch -- reference chain
    # tch -- template chain
    # in mtch_mtx rows are ref chs and cols are template chs
    for n, rch in enumerate(ref_chs):
        for m, tch in enumerate(trg_chs):
            ref = ref_stc.select("name CA and not HETERO and chain %s" % rch)
            tmp = trg_stc.select("name CA and not HETERO and chain %s" % tch)
            try:
                algs[key] = mth.execute(ref, tmp, **kwargs)
            except AlignError:
                continue
            mtch_mtx[n, m] = key
            key += 1

    # joining cabs chains
    pickups = []
    for n, refch in enumerate(mtch_mtx):
        inds = np.nonzero(refch)[0]
        picked_mers = set([])
        for ind in inds:
            trgmrs, dummy = zip(*algs[refch[ind]])
            if len(picked_mers.intersection(trgmrs)):
                mtch_mtx[n, ind] = 0
                continue
            picked_mers |= set(trgmrs)
            pickups.append(refch[ind])
            mtch_mtx[n + 1 :, ind] = 0

    try:
        trg_aln = reduce(operator.add, [algs.get(k, ()) for k in pickups])
    except TypeError:  # empty list of alignments --> no seq identity
        raise ValueError(
            "No sequential similarity between input and reference according to used alignment method (%s)."
            % align_mth
        )
    ref_mrs, tmp_mrs = zip(*trg_aln)
    ref_sstc = Atoms(arg=list(ref_mrs))
    tmp_sstc = Atoms(arg=list(tmp_mrs))
    # sstc -- selected substructure (only aligned part)
    return ref_sstc, tmp_sstc, trg_aln


class AlignError(Exception):
    pass


class AlnMeta(ABCMeta):
    def __init__(self, *args, **kwargs):
        super(AlnMeta, self).__init__(*args, **kwargs)
        try:
            self.methodname
        except AttributeError:
            self.methodname = self.__name__.lower()


class AbstractAlignMethod(metaclass=AlnMeta):
    @abstractmethod
    def execute(self, mers1, mers2, **kwargs):
        """Returns TUPLE of tuples containig subsequent aligned mers."""

    @classmethod
    def get_subclass_dict(cls):
        return {i.methodname: i() for i in cls.__subclasses__()}


class TrivialAlign(AbstractAlignMethod):
    methodname = "trivial"

    def execute(self, atoms1, atoms2, **kwargs):
        if len(atoms1) != len(atoms2):
            raise AlignError("Structure of different size passed to trivial alignment.")
        return tuple(zip(atoms1.atoms, atoms2.atoms))


class LoadCSVAlign(AbstractAlignMethod):
    methodname = "CSV"

    def execute(self, atoms1, atoms2, fname, **kwargs):
        try:
            with open(fname) as f:
                nms1, nms2 = zip(
                    *[i.split("\t") for i in map(str.strip, f.readlines()[1:])]
                )
        except TypeError:
            raise AlignError("No alignment file was given.")
        ats1 = [i for i in atoms1 if fmt_csv(i) in nms1]
        ats2 = [i for i in atoms2 if fmt_csv(i) in nms2]
        if 0 in map(len, (ats1, ats2)):
            raise AlignError("Empty alignment.")
        return zip(ats1, ats2)


class BLASTpAlign(AbstractAlignMethod):
    methodname = "blastp"

    @raise_aerror_on(StopIteration, IndexError)
    def execute(self, atoms1, atoms2, short=False, **kwargs):
        get_seq = lambda stc: [aa_to_short(r.resname) for r in stc.atoms]

        with NamedTemporaryFile(mode="w", delete=False) as f1:
            seq1 = get_seq(atoms1)
            f1.write(">tmp1\n%s" % "".join(seq1))
            tn1 = f1.name

        with NamedTemporaryFile(mode="w", delete=False) as f2:
            seq2 = get_seq(atoms2)
            f2.write(">tmp2\n%s" % "".join(seq2))
            tn2 = f2.name

        try:
            with open(tn1) as f1:
                with open(tn2) as f2:
                    task = ["-task", "blastp-short"] if short else []
                    res = check_output(
                        ["blastp", "-subject", f1.name, "-query", f2.name] + task
                    )  # ?? + kwargs
        finally:
            remove(tn1)
            remove(tn2)
        bhit = [
            i
            for i in map(str.strip, res.decode("utf-8").split("\n\n\n"))
            if i.startswith("Score")
        ][0].split("\n\n")[1:]
        mk_ind = lambda x: int(x) - 1
        aln_res = []
        for prt in bhit:
            query, match, subject = map(str.split, prt.split("\n"))
            indqs, indqe, indss, indse = map(
                mk_ind, (query[1], query[-1], subject[1], subject[-1])
            )
            it1 = iter(atoms1.atoms[indss : indse + 1])
            it2 = iter(atoms2.atoms[indqs : indse + 1])
            for i, j in zip(query[2], subject[2]):
                m1 = it1.next() if i.isalpha() else None
                m2 = it2.next() if i.isalpha() else None
                if None in (m1, m2):
                    continue
                aln_res.append((m1, m2))
        return tuple(aln_res)


class SmithWaterman(AbstractAlignMethod):
    methodname = "SW"

    def execute(self, atoms1, atoms2, ident_threshold=0.5, **kwargs):
        if len(atoms1) == 0 or len(atoms2) == 0:
            raise AlignError("Empty alignment.")
        # filling matrix
        mtx = np.zeros((len(atoms1) + 1, len(atoms2) + 1), dtype=int)
        for i, r1 in enumerate([aa_to_short(k.resname) for k in atoms1], 1):
            for j, r2 in enumerate([aa_to_short(kk.resname) for kk in atoms2], 1):
                dgn = mtx[i - 1, j - 1] + BLOSUM62[B62h[r1], B62h[r2]]
                dwn = mtx[i - 1, j] + BLOSUM62[B62h["null"], B62h[r2]]
                rght = mtx[i, j - 1] + BLOSUM62[B62h[r1], B62h["null"]]
                mtx[i, j] = max((dgn, dwn, rght, 0))
        # finding path
        if mtx.max() <= 0:
            raise AlignError("No sequential similarity.")
        i, j = np.unravel_index(np.argmax(mtx), mtx.shape)
        pickup = lambda ind1, ind2: (atoms1[ind1 - 1], atoms2[ind2 - 1])
        alg = [pickup(i, j)]
        while not (i == j == 0):
            di, dj = max(
                ((-1, -1), (-1, 0), (0, -1)),
                key=lambda x: float(mtx[i + x[0], j + x[1]]),
            )
            if 1 in (i, j):
                break  # ?? is it so
            i += di
            j += dj
            if di == dj:
                alg.append(pickup(i, j))
        result = tuple(reversed(alg))
        ident = (
            [i[0].resname == i[1].resname for i in result].count(True)
            * 1.0
            / min((len(atoms1), len(atoms2)))
        )
        if ident <= ident_threshold:
            raise AlignError("Identity below threshold.")
        return result
