"""Module to handle pdb files."""

import re
import os
import gzip
import json
import requests as req
from typing import Tuple

from tempfile import mkstemp
from time import sleep
from requests.exceptions import HTTPError, ConnectionError
from subprocess import Popen, PIPE
from collections import OrderedDict

from CABS import logger
from CABS.atom import Atom, Atoms
from CABS.utils import AA_NAMES, AA_SUB_NAMES
from CABS.plots import drop_csv_file

_name = 'PDB' # module name for logger
# PDB_CACHE = os.path.join(os.path.expanduser('~'), '.cabsPDBcache')
# try:
#     os.makedirs(PDB_CACHE)
# except FileExistsError:
#     pass


class Pdb(object):
    """
    Pdb parser.
    """

    DSSP_COMMAND = 'mkdssp'

    class InvalidPdbInput(Exception):
        pass

    def __init__(
            self,
            source,
            selection='',
            pdb_cache='~',
            remove_alternative_locations=True,
            fix_non_standard_aa=True,
            remove_water=True,
            remove_hetero=True,
            verify=False,
            no_exit=False,  # does not exit on error, raises InvalidPdbInput instead
            create_from_aa=False
    ):
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

            if chains:
                logger.debug(_name, 'Selected chains {}'.format(chains))
                if selection:
                    selection = '({}) and chain {}'.format(selection, ','.join(chains))
                else:
                    selection = 'chain {}'.format(','.join(chains))

            if remove_alternative_locations:
                if not create_from_aa:
                    logger.debug(_name, 'Removing alternative locations from {}'.format(name))
                self.atoms.remove_alternative_locations()

            if remove_water:
                if not create_from_aa:
                    logger.debug(_name, 'Removing water molecules from {}'.format(name))
                self.atoms = self.atoms.drop('resname HOH')

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
    def fetch(pdb_code, pdb_cache, force_download=False):

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
            url = f'http://files.rcsb.org/download/{pdb_low}.pdb.gz'
            r = req.get(url)
            r.raise_for_status()
            with open(filename, 'wb') as f:
                f.write(r.content)

        return filename

    @staticmethod
    def read(filename):
        try:
            with gzip.open(filename, 'rb') as f:
                content = f.read()
        except IOError:
            with open(filename, 'rb') as f:
                content = f.read()
        return content

    def run_dssp_command(self, command):
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

    def dssp(self, output='', dssp_from_aa=False):
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

        if output and logger.log_files():
            output_dssp = os.path.join(output, 'output_data', 'DSSP_output_%s.txt' % self.name)
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

        if output and logger.output_dssp():
            csv_output = os.path.join(output, 'output_data', 'DSSP_output_%s' % self.name)
            csv_dir = os.path.dirname(csv_output)
            if not os.path.isdir(csv_dir):
                os.makedirs(csv_dir)
            ss = self.atoms.get_dssp(sec)
            drop_csv_file(csv_output, [list(ss.keys()), list(ss.values())], fmts=["%s", "%s"])

        return sec

    def mk_ss_header(self, dssp_from_aa=False):
        dct = self.dssp(dssp_from_aa=dssp_from_aa)

        def mark(symbol):
            def getStart(xyz):
                """Returns 1 for residues that are at the beginning of seq of ss elements of the same type."""
                x, y, z = xyz
                if y[1] != symbol:
                    return 0
                if x[1] != y[1]:
                    return 1
                return 0

            def getEnd(xyz):
                """Returns 1 for residues that are at the end of seq of ss elements of the same type."""
                x, y, z = xyz
                if y[1] != symbol:
                    return 0
                if y[1] != z[1]:
                    return 1
                return 0

            return getStart, getEnd

        def szip(lst):
            lst = [None] + list(lst) + [None]
            return zip(lst, lst[1:], lst[2:])

        pickMidSS = lambda it: it[1][0]

        sseq = szip(dct.items())
        getSt, getEn = mark('H')
        rngH = zip(map(pickMidSS, filter(getSt, sseq)), map(pickMidSS, filter(getEn, sseq)))
        getSt, getEn = mark('E')
        rngS = zip(map(pickMidSS, filter(getSt, sseq)), map(pickMidSS, filter(getEn, sseq)))

        out = []

        serNum = 0
        helixID = ''
        hlxClass = 1
        comment = ''
        cas = self.atoms.select('NAME CA')
        recName = 'HELIX'
        # ~ for st, en in zip(rngH[::2], rngH[1::2]):
        for st, en in rngH:
            serNum += 1
            stNm, stChID = st.split(":")
            atSt = max(self.atoms.select('RESNUM %s' % stNm).select('CHAIN %s' % stChID).select('NAME CA'))
            enNm, enChID = en.split(":")
            atEn = max(self.atoms.select('RESNUM %s' % enNm).select('CHAIN %s' % enChID).select('NAME CA'))
            length = cas.atoms.index(atEn) - cas.atoms.index(atSt)
            inp = (recName, serNum, helixID, atSt.resname, stChID, atSt.resnum, atSt.icode, atEn.resname, enChID, atEn.resnum, atEn.icode, hlxClass, comment, length)
            line = "%-6s %3i %3s %3s %1s %4i%1s %3s %1s %4i%1s%2i%30s %5i\n" % inp
            out.append(line)

        recName = 'SHEET'
        serNum = 0
        sheetID = ''
        numStrs = 1
        sense = 0
        for st, en in rngS:
            serNum += 1
            stNm, stChID = st.split(":")
            atSt = max(self.atoms.select('RESNUM %s' % stNm).select('CHAIN %s' % stChID).select('NAME CA'))
            enNm, enChID = en.split(":")
            atEn = max(self.atoms.select('RESNUM %s' % enNm).select('CHAIN %s' % enChID).select('NAME CA'))
            inp = (recName, serNum, sheetID, numStrs, atSt.resname, stChID, atSt.resnum, atSt.icode, atEn.resname, enChID, atEn.resnum, atEn.icode, sense, '')
            line = "%-6s %3i %3s%2i %3s %1s%4i%1s %3s %1s%4i%1s%2i %29s\n" % inp
            out.append(line)

        return ''.join(out)

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

    def __str__(self):
        return self.body

    def __repr__(self):
        return "<PDB from %s, %i atoms>" % (self.name, len(self.atoms))
