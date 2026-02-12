"""Module for handling CABS simulation"""

from operator import attrgetter
import os
from random import randint
import re
from subprocess import PIPE, Popen
import subprocess
import sys
import tarfile
from threading import Thread
from typing import List, Tuple

import numpy as np
import numpy.typing as npt

try:
    from importlib.resources import as_file, files
except ImportError:
    import warnings

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", category=DeprecationWarning, module="pkg_resources"
        )
        from pkg_resources import resource_filename
from tempfile import mkdtemp
from time import strftime

from CABS import _JUNK
from CABS.core.trajectory import Trajectory
from CABS.io import logger
from CABS.structures.atom import Atoms
from CABS.structures.vector3d import Vector3d

_name = "CABS"


class CabsLattice:
    """
    This class represent a CABS-like lattice. It is initialized with:
    grid_spacing: distance between grid nodes, default 0.61 A
    r12: tuple with min and max allowed values for CA-CA pseudo-bond length
    r13: tuple with min and max allowed values for CA-CA-CA end distance
    """

    def __init__(
        self,
        grid_spacing: float = 0.61,
        r12: Tuple[float, float] = (3.28, 4.27),
        r13: Tuple[float, float] = (4.1, 7.35),
    ) -> None:
        self.grid: float = grid_spacing
        r12min = round((r12[0] / self.grid) ** 2)
        r12max = round((r12[1] / self.grid) ** 2)
        r13min = round((r13[0] / self.grid) ** 2)
        r13max = round((r13[1] / self.grid) ** 2)
        dim = int(r12max**0.5)

        self.vectors: List[Vector3d] = []
        for i in range(-dim, dim + 1):
            for j in range(-dim, dim + 1):
                for k in range(-dim, dim + 1):
                    _l = i * i + j * j + k * k
                    if r12min <= float(_l) <= r12max:
                        self.vectors.append(Vector3d(i, j, k))

        n = len(self.vectors)
        self.good: npt.NDArray[np.float64] = np.zeros((n, n))
        for i in range(n):
            vi = self.vectors[i]
            for j in range(n):
                vj = self.vectors[j]
                if r13min < (vi + vj).mod2() < r13max and vi.cross(vj).mod2():
                    self.good[i, j] = 1

    def cast(self, ch: "Atoms") -> List[Tuple[float, float, float]]:
        """
        Function that casts a single protein chain onto the lattice.
        Returns a list of tuples with (x, y, z) coordinates of CA atoms.
        """

        if len(ch.atoms) < 3:
            raise Exception("Protein chain too short!")

        prev = None
        coord = [
            Vector3d(
                round(ch.atoms[0].coord.x / self.grid),
                round(ch.atoms[0].coord.y / self.grid),
                round(ch.atoms[0].coord.z / self.grid),
            )
        ]

        for atom in ch.atoms[1:]:
            #  iterate over atoms
            min_dr = 1e12
            min_i = -1

            for i, v in enumerate(self.vectors):
                #  iterate over all possible vectors

                if len(coord) > 2 and self.good[prev, i] == 0:
                    continue

                new = coord[-1] + v
                dr = (self.grid * new - atom.coord).mod2()
                if dr < min_dr:
                    min_dr = dr
                    min_i = i

            if min_i < 0:
                raise Exception("Unsolvable geometric problem!")
            else:
                coord.append(coord[-1] + self.vectors[min_i])
                prev = min_i

        coord.insert(0, coord[0] + coord[1] - coord[2])
        coord.append(coord[-1] + coord[-2] - coord[-3])

        return coord


class CabsRun(Thread):
    """
    Class representing single cabs run.
    """

    LATTICE = (
        CabsLattice()
    )  # static object CabsLattice used to convert structures to CABS representation
    FORCE_FIELD = (
        4.0,
        1.0,
        2.0,
        0.125,
        -2.0,
        0.375,
    )  # parameters of the CABS force field
    FORTRAN_COMMAND = "gfortran -O2 -fno-trapping-math"
    CABS_DIR_FMT = "%y%m%d%H%M%S"

    def __init__(
        self,
        protein_complex,
        restraints,
        work_dir,
        replicas,
        replicas_dtemp,
        binding_interactions,
        mc_annealing,
        mc_cycles,
        mc_steps,
        temperature,
        ca_rest_weight,
        sc_rest_weight,
        excluding_distance,
        pairmod,
    ):
        """
        Initialize CabsRun object.
        :param protein_complex: ProteinComplex object with initial conformation of the complex (many replicas)
        :param restraints: Restraints object with complete list of CA-CA and SG-SG restraints
        """

        Thread.__init__(self)
        logger.debug(module_name=_name, msg="Loading structures...")
        fchains, seq, ids = CabsRun.load_structure(protein_complex)
        logger.debug(module_name=_name, msg="Loading restraints...")
        restr, maxres = CabsRun.load_restraints(
            restraints.update_id(ids), ca_rest_weight, sc_rest_weight
        )
        exclude = CabsRun.load_excluding(
            protein_complex.protein.exclude, excluding_distance, ids
        )

        ndim = max(protein_complex.chain_list.values()) + 2
        nmols = len(protein_complex.chain_list)
        nreps = replicas
        random_number = randint(999, 10000)
        inp = CabsRun.make_inp(
            random_number=random_number,
            nmols=nmols,
            force_field=CabsRun.FORCE_FIELD,
            replicas=replicas,
            replicas_dtemp=replicas_dtemp,
            binding_interactions=binding_interactions,
            mc_annealing=mc_annealing,
            mc_cycles=mc_cycles,
            mc_steps=mc_steps,
            temperature=temperature,
        )

        cabs_dir = mkdtemp(prefix="." + strftime(self.CABS_DIR_FMT), dir=work_dir)

        _JUNK.append(cabs_dir)

        if pairmod:
            epfile = ""
            for k, v in pairmod.items():
                mol, res = ids[protein_complex.new_ids[k]]
                weight, shift = v
                epfile += f"{mol} {res} {weight:.3f} {shift:.3f}\n"
            with open(os.path.join(cabs_dir, "EPAIRMOD"), "w") as f:
                f.write(epfile)

        with open(os.path.join(cabs_dir, "FCHAINS"), "w") as f:
            f.write(fchains)
        with open(os.path.join(cabs_dir, "SEQ"), "w") as f:
            f.write(seq)
        with open(os.path.join(cabs_dir, "INP"), "w") as f:
            f.write(inp + restr + exclude)

        logger.debug(module_name=_name, msg="Building exe...")

        # Get data file path using modern importlib.resources or fallback
        try:
            with as_file(files("CABS") / "data" / "data0.dat") as data_file:
                src_path = str(data_file)
        except (ImportError, AttributeError, NameError):
            src_path = resource_filename("CABS", "data/data0.dat")

        run_cmd = CabsRun.build_exe(
            params=(ndim, nreps, nmols, maxres),
            src=src_path,
            exe="cabs",
            build_command=self.FORTRAN_COMMAND,
            destination=cabs_dir,
        )

        # Extract data files using modern importlib.resources or fallback
        try:
            with as_file(files("CABS") / "data" / "data1.dat") as data_file:
                with tarfile.open(str(data_file)) as f:
                    f.extractall(cabs_dir)
        except (ImportError, AttributeError, NameError):
            with tarfile.open(resource_filename("CABS", "data/data1.dat")) as f:
                f.extractall(cabs_dir)

        self.cfg = {"cwd": cabs_dir, "exe": run_cmd}

    @staticmethod
    def load_structure(protein_complex):
        fchains = None
        seq = ""
        cabs_ids = {}
        for model in protein_complex.models():
            chains = model.chains()
            if not fchains:
                fchains = [""] * len(chains)
                ch = 1
                res = 1
                chid = model[0].chid
                for atom in model:
                    if atom.chid != chid:
                        res = 1
                        ch += 1
                        chid = atom.chid
                    cabs_ids[atom.resid_id()] = (ch, res)
                    res += 1
                    seq += "%5i%1s %1s%3s %1s%3i%6.2f\n" % (
                        atom.resnum,
                        atom.icode,
                        atom.alt,
                        atom.resname,
                        atom.chid,
                        int(atom.occ),
                        atom.bfac,
                    )

            for i, chain in enumerate(chains):
                vectors = CabsRun.LATTICE.cast(chain)
                fchains[i] += (
                    str(len(vectors))
                    + "\n"
                    + "\n".join([f"{int(v.x)} {int(v.y)} {int(v.z)}" for v in vectors])
                    + "\n"
                )

        return "".join(fchains), seq, cabs_ids

    @staticmethod
    def load_restraints(restraints, ca_weight=(1.0, 1.0), sg_weight=(1.0, 1.0)):
        max_r = 0

        rest = [r for r in restraints.data if not r.sg]
        restr = f"{len(rest)} {ca_weight[0]:.2f} {ca_weight[1]:.2f}\n"
        if len(rest):
            rest.sort(key=attrgetter("id2"))
            rest.sort(key=attrgetter("id1"))
            all_ids = [r.id1 for r in rest] + [r.id2 for r in rest]
            rest_count = {i: all_ids.count(i) for i in all_ids}
            max_r = max(1, *rest_count.values())
            rest = [
                "%2i %3i %2i %3i %8.4f %8.4f %6.2f %6.2f\n"
                % (
                    r.id1[0],
                    r.id1[1],
                    r.id2[0],
                    r.id2[1],
                    r.distance,
                    r.width,
                    r.weight_min,
                    r.weight_max,
                )
                for r in rest
            ]
            restr += "".join(rest)

        rest = [r for r in restraints.data if r.sg]
        restr += f"{len(rest)} {sg_weight[0]:.2f} {sg_weight[1]:.2f}\n"
        if len(rest):
            rest.sort(key=attrgetter("id2"))
            rest.sort(key=attrgetter("id1"))
            all_ids = [r.id1 for r in rest] + [r.id2 for r in rest]
            rest_count = {i: all_ids.count(i) for i in all_ids}
            max_r = max(max_r, *rest_count.values())
            rest = [
                f"{r.id1[0]:2d} {r.id1[1]:3d} {r.id2[0]:2d} {r.id2[1]:3d} {r.distance:8.4f} {r.width:8.4f} {r.weight_min:6.2f} {r.weight_max:6.2f}\n"
                for r in rest
            ]
            restr += "".join(rest)

        return restr, max_r

    @staticmethod
    def load_excluding(excl, dist, cabs_ids):
        if excl:
            return f"{len(excl)} {dist}\n" + "\n".join(
                f"{cabs_ids[i][0]} {cabs_ids[i][1]} {cabs_ids[j][0]} {cabs_ids[j][1]}"
                for i, j in excl
            )
        else:
            return "0 0.0\n"

    @staticmethod
    def build_exe(
        params, src, exe="cabs", build_command=FORTRAN_COMMAND, destination="."
    ):
        with open(src) as f:
            lines = f.read()

        names = ["NDIM", "NREPS", "NMOLS", "MAXRES"]
        for name, value in zip(names, params):
            lines = re.sub(rf"{name}=\d+", rf"{name}=%i" % value, lines)

        out_source = os.path.join(destination, "cabs_source.dat")  # TEMPORARY
        with open(out_source, "w") as f:
            f.write(lines)

        run_cmd = os.path.join(destination, exe)
        cmd = build_command.split() + ["-o", run_cmd, "-x", "f77", "-"]
        try:
            build_proc = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            out, err = build_proc.communicate(input=lines.encode("utf-8"))
            
            if build_proc.returncode != 0:
                stdout_msg = out.decode('utf-8') if out else ""
                stderr_msg = err.decode('utf-8') if err else ""
                msg = f"CABS binary compilation failed with return code {build_proc.returncode}.\nSTDOUT: {stdout_msg}\nSTDERR: {stderr_msg}"
                logger.critical(module_name=_name, msg=msg)
                raise Exception("CABS binary compilation failed")
            
            # On macOS (Darwin), binaries must be ad-hoc signed to execute after being compiled or modified locally. 
            if sys.platform == "darwin":
                import subprocess
                subprocess.run(["codesign", "-s", "-", "-f", run_cmd], check=True)

            if err:
                decoded_err = err.decode("utf-8")
                warning_message = decoded_err.split("\n")[-2] if "\n" in decoded_err else decoded_err
                logger.warning(module_name=_name, msg=warning_message)
        except OSError:
            logger.exit_program(_name, "Missing FORTRAN compiler (%s)." % cmd[0])
        except Exception as e:
            logger.critical(module_name=_name, msg=str(e))
            raise
        return run_cmd

    @staticmethod
    def make_inp(
        random_number,
        nmols,
        force_field,
        temperature,
        mc_annealing,
        mc_cycles,
        mc_steps,
        replicas,
        replicas_dtemp,
        binding_interactions,
    ):
        return f"""{random_number}
{mc_annealing} {mc_cycles} {mc_steps} {replicas} {nmols}
{temperature[0]:.2f} {temperature[1]:.2f} {force_field[0]:.2f} {binding_interactions:.2f} {replicas_dtemp:.2f}
{force_field[1]:.3f} {force_field[2]:.3f} {force_field[3]:.3f} {force_field[4]:.3f} {force_field[5]:.3f}
"""

    def run(self):
        monitor = logger.CabsObserver(
            interval=0.2, progress_file=os.path.join(self.cfg["cwd"], "PROGRESS")
        )
        try:
            cabs_proc = Popen(
                self.cfg["exe"], cwd=self.cfg["cwd"], stderr=PIPE, stdin=PIPE
            )
            stdout, stderr = cabs_proc.communicate()
            if stderr:
                logger.warning(module_name=_name, msg=stderr)
        except OSError:
            logger.exit_program(_name, "CABS binary missing or cannot be executed.")
        monitor.exit()

    def get_trajectory(self):
        traf = os.path.join(self.cfg["cwd"], "TRAF")
        seq = os.path.join(self.cfg["cwd"], "SEQ")
        return Trajectory.read_trajectory(traf, seq)


if __name__ == "__main__":
    pass
