import os
import glob
import re
import sys

from tempfile import mkstemp
from contextlib import closing
from CABS import logger
import sys


_name = 'MODELLER'
_PIR_TEMPLATE = '\n'.join(
    ['>P1;%s', 'sequence:::::::::', '%s', '*', '', '>P1;model_ca', 'structure:%s:FIRST:@:END:@::::', '*']
)

if sys.platform == 'win32':
	_PIR_TEMPLATE = _PIR_TEMPLATE.replace('\n', '\r\n');

try:
    from modeller import *
    from modeller.automodel import *
except ImportError as e:
    logger.exit_program(_name, 'MODELLER NOT FOUND')


def ca2all(
        filename, output=None, iterations=1, work_dir='.',
        out_mdl=os.path.join(os.getcwd(), 'output_data', 'modeller_output_0.txt')
):
    """
    Rebuilds ca to all-atom
    """

    old_stdout = sys.stdout
    if logger.log_files():
        sys.stdout = open(out_mdl, 'w')
    else:
        sys.stdout = open('/dev/null', 'w')

    pdb = mkstemp(prefix='.', suffix='.pdb', dir=work_dir, text=True)[1]
    prefix = os.path.basename(pdb).rsplit('.', 1)[0]

    aa_names = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
        'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
        'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
        'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }

    aa_names = {v: k for k, v in aa_names.items()}

    atoms = []
    pattern = re.compile('ATOM.{9}CA .([A-Z]{3}) ([A-Z ])(.{5}).{27}(.{12}).*')

    try:
        with closing(filename) as f, open(pdb, 'w') as tmp:
            for line in f:
                if line.startswith('ENDMDL'):
                    break
                else:
                    match = re.match(pattern, line)
                    if match:
                        atoms.append(match.groups())
                        tmp.write(line)

        if not len(atoms):
            raise Exception('File %s contains no CA atoms' % filename)
        chains = [atoms[0][1]]
        seq = ''
        for a in atoms:
            s, c = a[:2]
            if c not in chains:
                chains += c
                seq += '/'
            seq += aa_names[s]

        pir = prefix + '.pir'
        with open(pir, 'w') as f:
            f.write(_PIR_TEMPLATE % (prefix, seq, pdb))

        env = environ()
        env.io.atom_files_directory = ['.']

        class MyModel(automodel):
            def special_patches(self, aln):
                self.rename_segments(segment_ids=chains)

        mdl = MyModel(
            env,
            alnfile=pir,
            knowns='model_ca',
            sequence=prefix,
            assess_methods=assess.DOPE
        )

        mdl.md_level = refine.slow
        mdl.auto_align(matrix_file=prefix + '.mat')
        mdl.starting_model = 1
        mdl.ending_model = int(iterations)
        mdl.final_malign3d = True
        mdl.make()

        models = [m for m in mdl.outputs if m['failure'] is None]
        cmp_key = 'DOPE score'
        models.sort(lambda x, y: cmp(x[cmp_key], y[cmp_key]))
        final = models[0]['name'].rsplit('.', 1)[0] + '_fit.pdb'

        sys.stdout.close()
        sys.stdout = old_stdout

        if output:
            outfile = open(output, 'w')
        else:
            outfile = sys.stdout
        with open(final) as f:
            a = iter(atoms)
            current = ch = r = t = nl = None
            for line in f:
                if line.startswith('ATOM'):
                    res = line[21:27]
                    if not current or current != res:
                        current = res
                        ch, r, t = a.next()[1:]
                    nl = line[:21] + ch + r + line[27:54] + t
                    if len(line) > 66:
                        nl += line[66:]
                    outfile.write(nl)
                elif line.startswith('TER '):
                    outfile.write(line[:22] + nl[22:27] + '\n')
                else:
                    outfile.write(line)
    finally:
        junk = glob.glob(prefix + '*')
        try:
            map(os.remove, junk)
        except OSError as err:
            logger.warning(_name, err.message)
