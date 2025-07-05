"""Module to handle pdb files."""

import re
import os
import gzip
import json
import requests as req
from typing import Tuple, List, Dict, Union, Optional, Any, TextIO
from typing_extensions import Literal
from copy import deepcopy

from tempfile import mkstemp
from time import sleep
from requests.exceptions import HTTPError, ConnectionError
from subprocess import Popen, PIPE
from collections import OrderedDict

from CABS.io import logger
from CABS.structures.atom import Atom, Atoms
from CABS.constants import AA_NAMES, AA_SUB_NAMES
from CABS.analysis.plots import drop_csv_file

_name = 'PDB' # module name for logger
# PDB_CACHE = os.path.join(os.path.expanduser('~'), '.cabsPDBcache')
# try:
#     os.makedirs(PDB_CACHE)
# except FileExistsError:
#     pass


class Pdb:
    """
    Pdb parser.
    """

    DSSP_COMMAND: str = 'mkdssp'

    class InvalidPdbInput(Exception):
        pass

    def __init__(
            self,
            source: str,
            selection: str = '',
            pdb_cache: str = '~',
            remove_alternative_locations: bool = True,
            fix_non_standard_aa: bool = True,
            remove_water: bool = True,
            remove_hetero: bool = True,
            verify: bool = False,
            no_exit: bool = False,  # does not exit on error, raises InvalidPdbInput instead
            create_from_aa: bool = False
    ) -> None:
        if not create_from_aa:
            logger.debug(_name, 'Creating Pdb object from {}'.format(source))
        self.atoms = Atoms()

        words = source.split(':')
        try:
            name, rec, pep = words
            chains = rec + pep
        except ValueError:
            try:
                name, chains = words
            except ValueError:
                name = words[0]
                chains = ''

        try:
            self.body = self.read(name)
            self.name = os.path.basename(name).split('.')[0]
        except IOError:
            try:
                self.body = self.read(self.fetch(name, pdb_cache))
                self.name = name
            except ConnectionError as e:
                if no_exit:
                    raise Pdb.InvalidPdbInput(str(e))
                else:
                    logger.exit_program(
                        module_name=_name,
                        msg='Cannot connect to the PDB database',
                        exc=e
                    )
            except HTTPError as e:
                if no_exit:
                    raise Pdb.InvalidPdbInput(str(e))
                else:
                    logger.exit_program(
                        module_name=_name,
                        msg='Invalid PDB code: {}'.format(name),
                        exc=e
                    )
            except IOError as e:
                if no_exit:
                    raise Pdb.InvalidPdbInput(str(e))  # removed message
                else:
                    logger.exit_program(
                        module_name=_name,
                        msg='File {} not found'.format(name),
                        exc=e
                    )

        try:
            if not create_from_aa:
                logger.debug(_name, 'Processing {}'.format(name))
            current_model = 0
            # Ensure body is a string (handle both Python 2 and 3 compatibility)
            if isinstance(self.body, bytes):
                self.body = self.body.decode('utf-8')
            new_body = 'HEADER' + 74 * ' ' + '\n'
            for line in self.body.split('\n'):
                match = re.match(r'(ATOM|HETATM)', line)
                if match:
                    self.atoms.append(Atom(line, current_model))
                    new_body += line + '\n'
                else:
                    match = re.match(r'MODEL\s+(\d+)', line)
                    if match:
                        current_model = int(match.group(1))
                        new_body += line + '\n'
                    else:
                        match = re.match(r'(TER|ENDMDL)', line)
                        if match:
                            new_body += line + '\n'
            self.body = new_body

            if not len(self.atoms):
                raise Exception('{} contains no atoms'.format(source))

            if chains:
                if not create_from_aa:
                    logger.debug(_name, 'Selected chains {}'.format(chains))

                actual_chains = ''.join(self.atoms.list_chains().keys())
                extra_chains = set(chains) - set(actual_chains)

                if extra_chains:
                    msg = f'The chain ID(s): {" ".extra_chains} are not present in {name}'
                    raise Exception(msg)

                chains_selection = 'chain {}'.format(','.join(chains))
                self.atoms = self.atoms.select(chains_selection)

            if remove_alternative_locations:
                if not create_from_aa:
                    logger.debug(_name, 'Removing alternative locations from {}'.format(name))
                self.atoms.remove_alternative_locations()

            if remove_water:
                if not create_from_aa:
                    logger.debug(_name, 'Removing water molecules from {}'.format(name))
                self.atoms = self.atoms.drop('resname HOH')

            if not len(self.atoms):
                raise Exception('{} contains no atoms'.format(source))

            if fix_non_standard_aa:
                if not create_from_aa:
                    logger.debug(_name, 'Scanning {} for non-standard amino acids'.format(name))
                aa_names = [AA_NAMES[k] for k in AA_NAMES]
                for model in self.atoms.models():
                    for residue in model.residues():
                        resname = residue[0].resname
                        if resname not in aa_names:
                            if resname not in AA_SUB_NAMES:
                                logger.warning(
                                    _name, 'Unknown residue {} at {} in {}'.format(
                                        resname, residue[0].resid_id(), name
                                    )
                                )
                            else:
                                sub_name = AA_SUB_NAMES[resname]
                                for atom in residue:
                                    atom.resname = sub_name
                                    atom.hetatm = False
                                logger.warning(
                                    _name, 'Replacing {} -> {} for {} in {}'.format(
                                        resname, sub_name, residue[0].resid_id(), name
                                    )
                                )

            if remove_hetero:
                if not create_from_aa:
                    logger.debug(_name, 'Removing heteroatoms from {}'.format(name))
                self.atoms = self.atoms.drop('hetero')

            self.all_atoms = deepcopy(self.atoms)

            if selection:
                if not create_from_aa:
                    logger.debug(_name, 'Selecting [{}] from {}'.format(selection, name))
                self.atoms = self.atoms.select(selection)

            if ' ' in set([i.chid for i in self.atoms]):
                raise ValueError('Atoms with empty chain ID in selected part of PDB file detected.')

            if not len(self.atoms):
                raise Exception('{} contains no atoms'.format(source))

            if chains and verify:
                actual_chains = ''.join(self.atoms.list_chains().keys())
                logger.debug(
                    module_name=_name,
                    msg='Matching declared [{}] with actual [{}] chain IDs in {}.'.format(chains, actual_chains, name)
                )
                if set(chains) != set(actual_chains):
                    msg = 'Mismatch in chain IDs in {}: {} differs from {}'.format(name, chains, actual_chains)
                    logger.warning(_name, msg)
                    raise Exception(msg)

        except Exception as e:
            if no_exit:
                raise Pdb.InvalidPdbInput(str(e))
            else:
                logger.exit_program(
                    module_name=_name,
                    msg=str(e),
                    exc=e
                )

    @staticmethod
    def fetch(pdb_code: str, pdb_cache: str, force_download: bool = False) -> str:
        
        # Strip .pdb extension if present
        if pdb_code.lower().endswith('.pdb'):
            pdb_code = pdb_code[:-4]

        if not re.match(r'[1-9][0-9A-Za-z]{3}', pdb_code):
            raise IOError

        PDB_CACHE = os.path.join(os.path.expanduser(pdb_cache), '.cabsPDBcache')
        try:
            os.makedirs(PDB_CACHE)
        except FileExistsError:
            pass
        except OSError:
            PDB_CACHE = os.path.join(os.path.expanduser('~'), '.cabsPDBcache')
            try:
                os.makedirs(PDB_CACHE)
            except FileExistsError:
                pass

        pdb_low = pdb_code.lower()
        path = os.path.join(PDB_CACHE, pdb_low[1:3])
        try:
            os.makedirs(path)
        except OSError:
            pass

        filename = os.path.join(path, '%s.pdb.gz' % pdb_low)

        if not os.path.isfile(filename) or force_download:
            logger.debug(_name, 'Downloading {}'.format(pdb_low))
            url = f'https://files.rcsb.org/download/{pdb_low}.pdb.gz'
            r = req.get(url)
            r.raise_for_status()
            with open(filename, 'wb') as f:
                f.write(r.content)

        return filename

    @staticmethod
    def read(filename: str) -> str:
        try:
            with gzip.open(filename, 'rb') as f:
                content = f.read()
        except IOError:
            with open(filename, 'rb') as f:
                content = f.read()
        return content.decode('utf-8') if isinstance(content, bytes) else content

    def run_dssp_command(self, command: List[str]) -> Tuple[Optional[str], Optional[str], int]:
        """Run a subprocess command and return the output, error, and return code."""
        try:
            proc = Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            stdout, stderr = proc.communicate(input=self.body.encode('utf-8'))
            out = stdout.decode('utf-8')
            err = stderr.decode('utf-8')
            if stderr:
                if err.find('output-format') == -1:
                    logger.warning(
                        module_name=_name,
                        msg='DSSP ERROR: %s' % err.replace('\n', ' ')
                    )
                return None, None, -1
            return out, err, proc.returncode
        except OSError:
            return None, None, -1

    def dssp(self, work_dir: str = '', dssp_from_aa: bool = False) -> None:
        """Runs dssp on the read pdb file and returns a dictionary with secondary structure"""

        commands_to_try = [
            [self.DSSP_COMMAND, '--output-format', 'dssp', '/dev/stdin'],
            [self.DSSP_COMMAND, '/dev/stdin'],
            ['mkdssp', '--output-format', 'dssp', '/dev/stdin'],
            ['dssp', '/dev/stdin'],
        ]

        out, err, return_code = None, None, -1
        for command in commands_to_try:
            out, err, return_code = self.run_dssp_command(command)
            if return_code == 0:
                if not dssp_from_aa:
                    logger.debug(_name, 'DSSP successful')
                break

        if return_code != 0:
            logger.warning(
                module_name=_name,
                msg='DSSP was not ran at all.'
            )
            return None

        if work_dir and logger.output_dssp():
            output_dssp = os.path.join(work_dir, 'output_data', 'DSSP_output.txt')
            odir = os.path.dirname(output_dssp)
            if not os.path.isdir(odir):
                os.makedirs(odir)
            logger.to_file(
                filename=output_dssp,
                content=out,
                msg='Saving DSSP output to %s' % output_dssp
            )

        sec = OrderedDict()
        p = '^([0-9 ]{5}) ([0-9 ]{4}.)([A-Z ]) ([A-Z])  ([HBEGITSP ])(.*)$'

        for line in out.split('\n'):
            m = re.match(p, line)
            if m:
                key = str(m.group(2).strip() + ':' + m.group(3))
                if m.group(5) in 'HGIP':
                    val = 'H'
                elif m.group(5) in 'BE':
                    val = 'E'
                elif m.group(5) in 'T':
                    val = 'T'
                else:
                    val = 'C'
                sec[key] = val

        return sec

    def mk_ss_header(self, dssp_from_aa=False):
        dssp_data = self.dssp(dssp_from_aa=dssp_from_aa)

        def identify_boundaries(ss_type):
            def is_start(residue_triplet):
                """Returns 1 if a residue starts an SS sequence of the specified type."""
                prev, current, _ = residue_triplet
                if current[1] != ss_type:
                    return 0
                if prev[1] != current[1]:
                    return 1
                return 0

            def is_end(residue_triplet):
                """Returns 1 if a residue ends an SS sequence of the specified type."""
                _, current, next_ = residue_triplet
                if current[1] != ss_type:
                    return 0
                if current[1] != next_[1]:
                    return 1
                return 0

            return is_start, is_end

        def sliding_window(sequence):
            """Creates a triplet sliding window over a sequence."""
            sequence = [None] + list(sequence) + [None]
            return zip(sequence, sequence[1:], sequence[2:])

        extract_middle_residue = lambda triplet: triplet[1][0]

        residue_triplets = list(sliding_window(dssp_data.items()))

        # Identify helices (H)
        helix_start, helix_end = identify_boundaries('H')
        helix_ranges = list(zip(
            map(extract_middle_residue, filter(helix_start, residue_triplets)),
            map(extract_middle_residue, filter(helix_end, residue_triplets))
        ))

        # Identify sheets (E)
        sheet_start, sheet_end = identify_boundaries('E')
        sheet_ranges = list(zip(
            map(extract_middle_residue, filter(sheet_start, residue_triplets)),
            map(extract_middle_residue, filter(sheet_end, residue_triplets))
        ))

        output_lines = []

        # Generate HELIX records
        serial_number = 0
        helix_id = ''
        helix_class = 1
        helix_comment = ''
        ca_atoms = self.atoms.select('NAME CA')

        for start_residue, end_residue in helix_ranges:
            serial_number += 1
            start_num, start_chain = start_residue.split(":")
            start_atom = max(
                self.atoms.select('RESNUM %s' % start_num).select('CHAIN %s' % start_chain).select('NAME CA'))
            end_num, end_chain = end_residue.split(":")
            end_atom = max(self.atoms.select('RESNUM %s' % end_num).select('CHAIN %s' % end_chain).select('NAME CA'))
            helix_length = ca_atoms.atoms.index(end_atom) - ca_atoms.atoms.index(start_atom)
            helix_record = (
                "HELIX", serial_number, helix_id, start_atom.resname, start_chain,
                start_atom.resnum, start_atom.icode, end_atom.resname, end_chain,
                end_atom.resnum, end_atom.icode, helix_class, helix_comment, helix_length
            )
            line = "%-6s %3i %3s %3s %1s %4i%1s %3s %1s %4i%1s%2i%30s %5i\n" % helix_record
            output_lines.append(line)

        # Generate SHEET records
        serial_number = 0
        sheet_id = ''
        num_strands = 1
        strand_sense = 0

        for start_residue, end_residue in sheet_ranges:
            serial_number += 1
            start_num, start_chain = start_residue.split(":")
            start_atom = max(
                self.atoms.select('RESNUM %s' % start_num).select('CHAIN %s' % start_chain).select('NAME CA'))
            end_num, end_chain = end_residue.split(":")
            end_atom = max(self.atoms.select('RESNUM %s' % end_num).select('CHAIN %s' % end_chain).select('NAME CA'))
            sheet_record = (
                "SHEET", serial_number, sheet_id, num_strands, start_atom.resname,
                start_chain, start_atom.resnum, start_atom.icode, end_atom.resname,
                end_chain, end_atom.resnum, end_atom.icode, strand_sense, ''
            )
            line = "%-6s %3i %3s%2i %3s %1s%4i%1s %3s %1s%4i%1s%2i %29s\n" % sheet_record
            output_lines.append(line)

        return ''.join(output_lines)

    @staticmethod
    def xssp(filename, server='https://www3.cmbi.umcn.nl/xssp'):
        url_api = server + '/api/%s/pdb_file/dssp/'

        files = {'file_': open(filename, 'rb')}

        r = req.post(url=url_api % 'create', files=files)
        r.raise_for_status()
        job_id = json.loads(r.content)['id']
        while True:
            r = req.get(url_api % 'status' + job_id)
            r.raise_for_status()
            status = json.loads(r.content)['status']

            if status == 'SUCCESS':
                r = req.get(url_api % 'result' + job_id)
                r.raise_for_status()
                out = json.loads(r.content)['result']
                err = ''
                break
            elif status in ['FAILURE', 'REVOKED']:
                err = json.loads(r.content)['message']
                out = ''
                break
            else:
                sleep(1)

        return out, err

    def save_initial_pdb(self, work_dir=''):
        if work_dir:
            initial_pdb = os.path.join(work_dir, 'output_pdbs', 'start_all.pdb')
            odir = os.path.dirname(initial_pdb)
            if not os.path.isdir(odir):
                os.makedirs(odir)
            self.all_atoms.save_to_pdb(initial_pdb)

    def __str__(self) -> str:
        return self.body

    def __repr__(self) -> str:
        return f"<PDB from {self.name}, {len(self.atoms)} atoms>"
