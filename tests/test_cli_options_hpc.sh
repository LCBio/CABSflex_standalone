#!/bin/bash
# ==============================================================================
# Comprehensive Parallel Test Suite for CABSflex - True User Simulation (v8)
# ==============================================================================

# --- Configuration ---
GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'
# Dynamically read variables from install-hpc.sh
if [ -f "install-hpc.sh" ]; then
    ENV_NAME=$(grep '^VENV_NAME=' install-hpc.sh | cut -d'=' -f2 | tr -d '"' | tr -d "'")
    GCC_MOD=$(grep '^GCC_MODULE=' install-hpc.sh | cut -d'=' -f2 | tr -d '"' | tr -d "'")
    PY_MOD=$(grep '^PYTHON_MODULE=' install-hpc.sh | cut -d'=' -f2 | tr -d '"' | tr -d "'")
else
    # Fallback defaults
    ENV_NAME="cabs_021"
    GCC_MOD="GCCcore/13.2.0"
    PY_MOD="Python/3.11.5"
fi

CONDA_BASE=$(conda info --base)

# --- HPC Environment Setup ---
# Load necessary modules (dynamically from install-hpc.sh)
if command -v module &> /dev/null; then
    module purge
    module load "$GCC_MOD" "$PY_MOD"
    # Add other modules if needed
fi

# Activate environment
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate $ENV_NAME

# Ensure local CABS is used. Package was installed with pip install -e .
CABS_CMD="CABSflex"

SCENARIO_ROOT="tests/test_cli_options_hpc"
RESULTS_FILE="tests/test_cli_options_hpc_results.log"
rm -rf "$SCENARIO_ROOT" "$RESULTS_FILE"
mkdir -p "$SCENARIO_ROOT"

# NO SPEED-UPS. Using tool defaults (e.g. -a 50, -y 50, -s 10).
# Total steps = 25,000 per simulation.
# Note: Since simulations take longer, we use fewer workers for wait logic if needed,
# but we'll try standard backgrounding first.
BASE_PDB="tests/inputs/Helix_short.pdb"
DS_PDB="tests/inputs/3BCI_ds.pdb"
REST_FILE="tests/inputs/restraints_short.txt"
MOD_CMD="mod10.8"

# Function to run an individual test
run_exhaustive_test() {
    local alias_tested="$1"
    local cmd="$2"
    local name=$(echo "$alias_tested" | tr -d '-' | tr ' ' '_' | tr '/' '_')
    local out_dir="${SCENARIO_ROOT}/${name}"
    
    mkdir -p "$out_dir"
    
    # Run command. No speed args added.
    eval "$cmd -w $out_dir > ${out_dir}/stdout.log 2> ${out_dir}/stderr.log"
    local status=$?
    
    local warnings=""
    local errors=""
    
    # Analyze CABS.log if it exists
    if [ -f "${out_dir}/CABS.log" ]; then
        warnings=$(grep "WARNING" "${out_dir}/CABS.log" | head -n 2 | tr '\n' ';' | sed 's/;$//')
        errors=$(grep "ERROR" "${out_dir}/CABS.log" | head -n 2 | tr '\n' ';' | sed 's/;$//')
        critical=$(grep "CRITICAL" "${out_dir}/CABS.log" | head -n 2 | tr '\n' ';' | sed 's/;$//')
        [ -n "$critical" ] && errors="${errors};${critical}"
    fi

    # Check stderr for system/python errors
    if [ $status -ne 0 ]; then
        # Capture last line of stderr which usually has the error message
        stderr_err=$(tail -n 1 "${out_dir}/stderr.log" | head -c 150 | tr '\n' ' ' | sed 's/|/ /g')
        errors="${errors};$stderr_err"
    fi

    local result="PASSED"
    [ $status -ne 0 ] && result="FAILED"
    
    # Escape pipe for markdown
    warnings=$(echo "$warnings" | sed 's/|/\\|/g' | tr '\n' ' ')
    errors=$(echo "$errors" | sed 's/|/\\|/g' | tr '\n' ' ')

    echo "$alias_tested|$warnings|$errors|$result" > "${out_dir}/result.tmp"
}

echo -e "${BLUE}🚀 Starting CABSflex Exhaustive CLI Tests (v8) - True Simulation${NC}"
echo -e "${YELLOW}Warning: This will take several minutes as it uses default simulation lengths.${NC}"

# 1. Basic & Misc (Individual)
run_exhaustive_test "-i" "$CABS_CMD -i $BASE_PDB" &
run_exhaustive_test "--input-protein" "$CABS_CMD --input-protein $BASE_PDB" &
run_exhaustive_test "-c" "$CABS_CMD -i $BASE_PDB -c tests/inputs/config.ini" &
run_exhaustive_test "--config" "$CABS_CMD -i $BASE_PDB --config tests/inputs/config.ini" &
run_exhaustive_test "-v" "$CABS_CMD -i $BASE_PDB -v 4" &
run_exhaustive_test "--verbose" "$CABS_CMD -i $BASE_PDB --verbose 4" &
run_exhaustive_test "--log" "$CABS_CMD -i $BASE_PDB --log" &
run_exhaustive_test "--version" "$CABS_CMD --version" &
run_exhaustive_test "-h" "$CABS_CMD -h" &
run_exhaustive_test "--help" "$CABS_CMD --help" &

# 2. Protein Options
run_exhaustive_test "-g" "$CABS_CMD -i $BASE_PDB -g flexible 3 3.8 11.5" &
run_exhaustive_test "--protein-restraints" "$CABS_CMD -i $BASE_PDB --protein-restraints flexible 3 3.8 11.5" &
run_exhaustive_test "-N" "$CABS_CMD -i $BASE_PDB -N" &
run_exhaustive_test "--no-protein-restraints" "$CABS_CMD -i $BASE_PDB --no-protein-restraints" &
run_exhaustive_test "--protein-restraints-retain" "$CABS_CMD -i $BASE_PDB --protein-restraints-retain 50" &
run_exhaustive_test "--protein-plddt" "$CABS_CMD -i $BASE_PDB --protein-plddt tests/inputs/flex.txt" &
run_exhaustive_test "--protein-category" "$CABS_CMD -i $BASE_PDB --protein-category tests/inputs/flex.txt" &
run_exhaustive_test "-f" "$CABS_CMD -i $BASE_PDB -f 0.5" &
run_exhaustive_test "--protein-flexibility" "$CABS_CMD -i $BASE_PDB --protein-flexibility 0.5" &
run_exhaustive_test "--weighted-fit" "$CABS_CMD -i $BASE_PDB --weighted-fit ss" &
run_exhaustive_test "--gauss-iterations" "$CABS_CMD -i $BASE_PDB --weighted-fit gauss --gauss-iterations 5" &
run_exhaustive_test "--receptor-ss" "$CABS_CMD -i $BASE_PDB --receptor-ss A:CCCHHHHHHHHHHHCCC" &
run_exhaustive_test "--peptide-structure-prediction" "$CABS_CMD -i HKILHRLLQD --peptide-structure-prediction -k 1" &

# 3. Restraints Options
run_exhaustive_test "--ca-rest-add" "$CABS_CMD -i $BASE_PDB --ca-rest-add 44:A 50:A 7.0 1.0" &
run_exhaustive_test "--sc-rest-add" "$CABS_CMD -i $BASE_PDB --sc-rest-add 44:A 50:A 5.0 1.0" &
run_exhaustive_test "--ca-rest-weight" "$CABS_CMD -i $BASE_PDB --ca-rest-weight 1.0 1.5" &
run_exhaustive_test "--sc-rest-weight" "$CABS_CMD -i $BASE_PDB --sc-rest-weight 0.5 1.0" &
run_exhaustive_test "--ca-rest-file" "$CABS_CMD -i $BASE_PDB --ca-rest-file $REST_FILE" &
run_exhaustive_test "--sc-rest-file" "$CABS_CMD -i $BASE_PDB --sc-rest-file $REST_FILE" &
run_exhaustive_test "-F" "$CABS_CMD -i $DS_PDB -F 26:A 29:A" &
run_exhaustive_test "--disulfide-bonds" "$CABS_CMD -i $DS_PDB --disulfide-bonds 26:A 29:A" &
run_exhaustive_test "--backbone-cyclization" "$CABS_CMD -i $BASE_PDB --backbone-cyclization A" &

# 4. Simulation Options (TRUE USER LOGIC)
run_exhaustive_test "-a" "$CABS_CMD -i $BASE_PDB -a 10" &
run_exhaustive_test "--mc-annealing" "$CABS_CMD -i $BASE_PDB --mc-annealing 10" &
run_exhaustive_test "-y" "$CABS_CMD -i $BASE_PDB -y 10" &
run_exhaustive_test "--mc-cycles" "$CABS_CMD -i $BASE_PDB --mc-cycles 10" &
run_exhaustive_test "-s" "$CABS_CMD -i $BASE_PDB -s 2" &
run_exhaustive_test "--mc-steps" "$CABS_CMD -i $BASE_PDB --mc-steps 2" &
run_exhaustive_test "-r" "$CABS_CMD -i $BASE_PDB -r 2" &
run_exhaustive_test "--replicas" "$CABS_CMD -i $BASE_PDB --replicas 2" &
run_exhaustive_test "-D" "$CABS_CMD -i $BASE_PDB -D 0.1" &
run_exhaustive_test "--replicas-dtemp" "$CABS_CMD -i $BASE_PDB --replicas-dtemp 0.1" &
run_exhaustive_test "-t" "$CABS_CMD -i $BASE_PDB -t 2.0 1.0" &
run_exhaustive_test "--temperature" "$CABS_CMD -i $BASE_PDB --temperature 2.0 1.0" &
run_exhaustive_test "-z" "$CABS_CMD -i $BASE_PDB -z 123" &
run_exhaustive_test "--random-seed" "$CABS_CMD -i $BASE_PDB --random-seed 123" &
run_exhaustive_test "-b" "$CABS_CMD -i $BASE_PDB -b 1.1" &
run_exhaustive_test "--binding-interactions" "$CABS_CMD -i $BASE_PDB --binding-interactions 1.1" &
run_exhaustive_test "--disable-centro" "$CABS_CMD -i $BASE_PDB --disable-centro" &

# 5. AA Reconstruction
run_exhaustive_test "-A" "$CABS_CMD -i $BASE_PDB -A --aa-method cg2all" &
run_exhaustive_test "--aa-rebuild" "$CABS_CMD -i $BASE_PDB --aa-rebuild --aa-method cg2all" &
run_exhaustive_test "--aa-method" "$CABS_CMD -i $BASE_PDB -A --aa-method cg2all" &
run_exhaustive_test "-m" "$CABS_CMD -i $BASE_PDB -A --aa-method modeller -m 1" &
run_exhaustive_test "--modeller-iterations" "$CABS_CMD -i $BASE_PDB -A --aa-method modeller --modeller-iterations 1" &

# 6. Analysis Options
# Using defaults for -a, -y etc.
run_exhaustive_test "-R" "$CABS_CMD -i $BASE_PDB -R tests/inputs/Helix.pdb" &
run_exhaustive_test "--reference-pdb" "$CABS_CMD -i $BASE_PDB --reference-pdb tests/inputs/Helix.pdb" &
run_exhaustive_test "-k" "$CABS_CMD -i $BASE_PDB -k 2" &
run_exhaustive_test "--clustering-medoids" "$CABS_CMD -i $BASE_PDB --clustering-medoids 2" &
run_exhaustive_test "--clustering-iterations" "$CABS_CMD -i $BASE_PDB --clustering-iterations 10" &
run_exhaustive_test "-n" "$CABS_CMD -i $BASE_PDB -n 10" &
run_exhaustive_test "--filtering-count" "$CABS_CMD -i $BASE_PDB --filtering-count 10" &
run_exhaustive_test "--filtering-mode" "$CABS_CMD -i $BASE_PDB --filtering-mode each" &
run_exhaustive_test "-M" "$CABS_CMD -i $BASE_PDB -M" &
run_exhaustive_test "--contact-maps" "$CABS_CMD -i $BASE_PDB --contact-maps" &
run_exhaustive_test "-T" "$CABS_CMD -i $BASE_PDB -M -T 6.0" &
run_exhaustive_test "--contact-threshold" "$CABS_CMD -i $BASE_PDB -M --contact-threshold 6.0" &
run_exhaustive_test "--contact-threshold-aa" "$CABS_CMD -i $BASE_PDB -M --contact-threshold-aa 5.0" &
run_exhaustive_test "--contact-map-colors" "$CABS_CMD -i $BASE_PDB -M --contact-map-colors '#000000' '#111111' '#222222' '#333333' '#444444' '#555555'" &
run_exhaustive_test "--align" "$CABS_CMD -i $BASE_PDB -R tests/inputs/Helix.pdb --align SW" &
run_exhaustive_test "--align-options" "$CABS_CMD -i $BASE_PDB -R tests/inputs/Helix.pdb --align-options gap_open=10" &

# 7. Output Options
run_exhaustive_test "-S" "$CABS_CMD -i $BASE_PDB -S" &
run_exhaustive_test "--save-cabs-files" "$CABS_CMD -i $BASE_PDB --save-cabs-files" &
run_exhaustive_test "-C" "$CABS_CMD -i $BASE_PDB -C" &
run_exhaustive_test "--save-config" "$CABS_CMD -i $BASE_PDB --save-config" &
run_exhaustive_test "-o" "$CABS_CMD -i $BASE_PDB -o RM" &
run_exhaustive_test "--pdb-output" "$CABS_CMD -i $BASE_PDB --pdb-output RM" &
run_exhaustive_test "--pdb-bfac-output" "$CABS_CMD -i $BASE_PDB --pdb-bfac-output BS" &
run_exhaustive_test "--csv-output" "$CABS_CMD -i $BASE_PDB --csv-output BS" &
run_exhaustive_test "--json-output" "$CABS_CMD -i $BASE_PDB --json-output" &
run_exhaustive_test "--dssp-output" "$CABS_CMD -i $BASE_PDB --dssp-output" &
run_exhaustive_test "--ss-output" "$CABS_CMD -i $BASE_PDB --ss-output" &
run_exhaustive_test "--restraints-output" "$CABS_CMD -i $BASE_PDB --restraints-output" &
run_exhaustive_test "--renumber-residues-to-original" "$CABS_CMD -i $BASE_PDB --renumber-residues-to-original" &

# 8. Misc
run_exhaustive_test "-w" "$CABS_CMD -i $BASE_PDB -w tests/scenarios_exhaustive/w_test" &
run_exhaustive_test "--work-dir" "$CABS_CMD -i $BASE_PDB --work-dir tests/scenarios_exhaustive/workdir_test" &
run_exhaustive_test "--dssp-command" "$CABS_CMD -i $BASE_PDB --dssp-command mkdssp" &
run_exhaustive_test "--fortran-command" "$CABS_CMD -i $BASE_PDB --fortran-command gfortran" &
run_exhaustive_test "--image-file-format" "$CABS_CMD -i $BASE_PDB -M --image-file-format png" &
run_exhaustive_test "--pdb-cache-dir" "$CABS_CMD -i $BASE_PDB --pdb-cache-dir ." &
run_exhaustive_test "--nsp3-model-path" "$CABS_CMD -i $BASE_PDB --nsp3-model-path dummy_weights" &

# 9. Previously Missing Options (User Requested)
run_exhaustive_test "-L" "$CABS_CMD -i $BASE_PDB -L tests/inputs/Test.cbs" &
run_exhaustive_test "--load-cabs-files" "$CABS_CMD -i $BASE_PDB --load-cabs-files tests/inputs/Test.cbs" &

# Dynamically find cg2all path for testing the CLI override option
TEST_CG2ALL_PATH=$(python3 -c "from CABS.config_loader import get_cg2all_env_prefix; print(get_cg2all_env_prefix() or '/dummy/path')" 2>/dev/null)
run_exhaustive_test "--cg2all-env-prefix" "$CABS_CMD -i $BASE_PDB -A --aa-method cg2all --cg2all-env-prefix $TEST_CG2ALL_PATH" &

# NOTE: Options like -p, -P, -e, -d, --pairmod are CABSdock specific and not valid for CABSflex.

echo -e "${YELLOW}Waiting for all exhaustive tests (v8) to complete. Parallel workers: unlimited.${NC}"
wait

cat ${SCENARIO_ROOT}/*/result.tmp > "$RESULTS_FILE"

echo -e "\n${BLUE}========================================${NC}"
echo -e "${BLUE}CABSflex CLI EXHAUSTIVE ALIAS REPORT (v8)${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "| Option Tested | Warnings | Errors | Passed/Failed |"
echo -e "|---------------|----------|--------|---------------|"

sort "$RESULTS_FILE" | while IFS='|' read -r opt warn err res; do
    if [ "$res" == "FAILED" ]; then
        echo -e "| $opt | $warn | ${RED}$err${NC} | ${RED}$res${NC} |"
    else
        echo -e "| $opt | $warn | $err | ${GREEN}$res${NC} |"
    fi
done
echo -e "${BLUE}========================================${NC}"