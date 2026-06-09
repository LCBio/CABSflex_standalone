import os
from string import Template

PML_RESTRAINTS_TEMPLATE = Template('''\
load $start_pdb, start
hide all
show cartoon, start
$restraints_lines
''')

PML_MODELS_TEMPLATE = Template('''\
$load_lines
$align_lines
hide all
dss
show cartoon
# Optionally, group the models
# group models, model_*
''')

PML_COLOR_SS_TEMPLATE = Template('''\
# Color by CABS secondary structure values stored in start_secstr.pdb B-factors
# (1=Coil/gray, 2=Helix/purple, 3=Turn/green, 4=Sheet/orange)
# dss is intentionally omitted: it would overwrite the CABS B-factor SS codes
# with PyMOL's own geometry-based assignment.
load $target_pdb, secondary_structure_parameters
hide all
show cartoon
spectrum b, blue_white_red, secondary_structure_parameters, minimum=1, maximum=4
''')

PML_COLOR_RMSF_TEMPLATE = Template('''\
# Color by B-factor (RMSF/flexibility)
$load_lines
hide all
dss
show cartoon
spectrum b, rainbow_rev, minimum=0.0
''')

PML_ANIMATE_TEMPLATE = Template('''\
# Load all models into a single object with multiple states
$load_states_lines
hide all
dss
show cartoon
smooth
mset 1 -$n_states
mplay
''')

def generate_pymol_scripts(work_dir, start_pdb_path, models_pdbs, restraints_file=None):
    """
    Generates PyMOL visualization scripts (.pml) in the given work directory.
    
    Args:
        work_dir (str): Working directory to save the scripts.
        start_pdb_path (str): Path to the starting/reference PDB file.
        models_pdbs (list of str): List of paths to the model PDB files.
        restraints_file (str, optional): Path to the restraints.txt file.
    """
    scripts_dir = work_dir # Or a subdirectory, but work_dir is fine for easy access
    
    _generate_load_models_script(scripts_dir, start_pdb_path, models_pdbs)
    _generate_color_ss_script(scripts_dir, models_pdbs)
    _generate_color_rmsf_script(scripts_dir, models_pdbs)
    _generate_animate_script(scripts_dir, models_pdbs)
    
    if restraints_file and os.path.exists(restraints_file):
        _generate_restraints_script(scripts_dir, start_pdb_path, restraints_file)


def _generate_load_models_script(scripts_dir, start_pdb_path, models_pdbs):
    load_lines = []
    align_lines = []
    
    for i, model_path in enumerate(models_pdbs):
        rel_path = os.path.relpath(model_path, scripts_dir)
        load_lines.append(f"load {rel_path}, model_{i}")
        if i > 0:
            align_lines.append(f"align model_{i}, model_0")
            
    content = PML_MODELS_TEMPLATE.substitute(
        load_lines="\n".join(load_lines),
        align_lines="\n".join(align_lines)
    )
    with open(os.path.join(scripts_dir, "load_models.pml"), "w") as f:
        f.write(content)

def _generate_color_ss_script(scripts_dir, models_pdbs):
    target_pdb = ""
    if models_pdbs:
        first_model_dir = os.path.dirname(models_pdbs[0])
        secstr_pdb = os.path.join(first_model_dir, "start_secstr.pdb")
        if os.path.exists(secstr_pdb):
            target_pdb = secstr_pdb
        else:
            target_pdb = models_pdbs[0]

    rel_path = os.path.relpath(target_pdb, scripts_dir) if target_pdb else ""
    content = PML_COLOR_SS_TEMPLATE.substitute(
        target_pdb=rel_path
    )
    with open(os.path.join(scripts_dir, "color_by_ss.pml"), "w") as f:
        f.write(content)

def _generate_color_rmsf_script(scripts_dir, models_pdbs):
    # Depending on CABS output, b-factor output might be in start_rmsf.pdb 
    # instead of the models directly. Let's load the model equivalents or b-factor ones if known.
    # For now, we assume the B-factor column in the standard models contains 
    # the relevant data or PyMOL will use existing B-factors.
    # We should use start_rmsf.pdb or start_bfac.pdb if possible, but let's stick to models if we don't pass the specific bf-pdbs.
    load_lines = []
    
    # Check if we have start_rmsf.pdb in the same dir as the first model
    if models_pdbs:
        first_model_dir = os.path.dirname(models_pdbs[0])
        rmsf_pdb = os.path.join(first_model_dir, "start_rmsf.pdb")
        bfac_pdb = os.path.join(first_model_dir, "start_bfac.pdb")
        
        target_pdb = ""
        if os.path.exists(rmsf_pdb):
            target_pdb = rmsf_pdb
        elif os.path.exists(bfac_pdb):
            target_pdb = bfac_pdb
            
        if target_pdb:
            rel_path = os.path.relpath(target_pdb, scripts_dir)
            load_lines.append(f"load {rel_path}, flexible_model")
        else:
            # Fallback to models
            for i, model_path in enumerate(models_pdbs):
                rel_path = os.path.relpath(model_path, scripts_dir)
                load_lines.append(f"load {rel_path}, model_{i}")
    
    content = PML_COLOR_RMSF_TEMPLATE.substitute(
        load_lines="\n".join(load_lines)
    )
    with open(os.path.join(scripts_dir, "color_by_rmsf.pml"), "w") as f:
        f.write(content)

def _generate_animate_script(scripts_dir, models_pdbs):
    # Load each model into the same object without an explicit state index so
    # PyMOL appends them sequentially.  Specifying state=N is unreliable across
    # PyMOL versions and can leave gaps when models_pdbs is not contiguous.
    load_states_lines = []
    for model_path in models_pdbs:
        rel_path = os.path.relpath(model_path, scripts_dir)
        load_states_lines.append(f"load {rel_path}, animation_sequence")

    n_states = max(len(models_pdbs), 1)
    content = PML_ANIMATE_TEMPLATE.substitute(
        load_states_lines="\n".join(load_states_lines),
        n_states=n_states,
    )
    with open(os.path.join(scripts_dir, "animate_models.pml"), "w") as f:
        f.write(content)

def _generate_restraints_script(scripts_dir, start_pdb_path, restraints_file):
    rel_start_pdb = os.path.relpath(start_pdb_path, scripts_dir)
    restraints_lines = []
    
    try:
        with open(restraints_file, "r") as f:
            for line_idx, line in enumerate(f):
                parts = line.split()
                if len(parts) >= 2:
                    # Format: res1:chain1 res2:chain2 distance weight something
                    res_chain1 = parts[0].split(':')
                    res_chain2 = parts[1].split(':')
                    
                    if len(res_chain1) == 2 and len(res_chain2) == 2:
                        res1, chain1 = res_chain1
                        res2, chain2 = res_chain2
                        
                        cmd = (f"distance rest_{line_idx}, "
                               f"start and chain {chain1} and resi {res1} and name CA, "
                               f"start and chain {chain2} and resi {res2} and name CA")
                        restraints_lines.append(cmd)
    except Exception as e:
        # Silently fail or log if needed, generating what we can
        pass
        
    content = PML_RESTRAINTS_TEMPLATE.substitute(
        start_pdb=rel_start_pdb,
        restraints_lines="\n".join(restraints_lines)
    )
    with open(os.path.join(scripts_dir, "load_restraints.pml"), "w") as f:
        f.write(content)
