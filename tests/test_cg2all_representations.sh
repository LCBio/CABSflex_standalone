#!/bin/bash
# ==============================================================================
#           cg2all representation tests for CABSflex and CABSdock
# ==============================================================================

# --- Configuration ---
set +e
GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

# Dynamically read ENV_NAME from install.sh
if [ -f "install.sh" ]; then
    ENV_NAME=$(grep '^ENV_NAME=' install.sh | cut -d'=' -f2 | tr -d '"' | tr -d "'")
else
    ENV_NAME="cabs" # Fallback
fi

# Enable job control for 'jobs' command to work in script
set -m

# Micromamba activation support
export MAMBA_EXE="$HOME/.local/bin/micromamba"
if [ -f "$MAMBA_EXE" ]; then
    eval "$($MAMBA_EXE shell hook -s bash)"
    micromamba activate $ENV_NAME
else
    # Fallback to conda if micromamba not found
    CONDA_BASE=$(conda info --base)
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate $ENV_NAME
fi

CABSFLEX_CMD="CABSflex"
CABSDOCK_CMD="CABSdock"

SCENARIO_ROOT="tests/test_cg2all_representations"
RESULTS_FILE="tests/test_cg2all_representations_results.log"
DEBUG_LOG="${SCENARIO_ROOT}/debug_progress.log"

rm -rf "$SCENARIO_ROOT" "$RESULTS_FILE"
mkdir -p "$SCENARIO_ROOT"
touch "$DEBUG_LOG"

# Full simulations only: do not pass -a, -y, -s, -r, -k, or other speed-up options.
FLEX_PDB="tests/inputs/2BZ6.pdb"
DOCK_PDB="tests/inputs/1A2K.pdb"

# All output options enabled for the test
OUTPUT_OPTS="--save-cabs-files --save-config --json-output --dssp-output --ss-output --restraints-output --renumber-residues-to-original --generate-chimera-visualizations --generate-pymol-visualizations --generate-notebook --pdb-output A --pdb-bfac-output A --csv-output A --contact-maps -v 4"

# Concurrency limiting logic (from tests/test_cli_options.sh)
if [[ "$OSTYPE" == "darwin"* ]]; then
    NPROC=$(sysctl -n hw.ncpu)
else
    NPROC=$(nproc)
fi
# Use NPROC-1 workers, min 1
MAX_JOBS=$((NPROC - 1))
[ $MAX_JOBS -lt 1 ] && MAX_JOBS=1
echo "DEBUG: MAX_JOBS=$MAX_JOBS"

# Trap signals to ensure we don't just exit silently
trap 'echo "Received signal, killing children"; jobs -p | xargs -r kill; exit 1' SIGINT SIGTERM

limit_jobs() {
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 1
    done
    true
}

run_test() {
    local case_name="$1"
    local cmd="$2"
    local out_dir="${SCENARIO_ROOT}/${case_name}"
    
    mkdir -p "$out_dir"
    
    echo "Running: $case_name"
    # Execute the command
    eval "$cmd -w $out_dir > ${out_dir}/stdout.log 2> ${out_dir}/stderr.log"
    local status=$?
    
    local result="PASSED"
    [ $status -ne 0 ] && result="FAILED"
    
    # Simple result logging
    echo "$case_name|$result" > "${out_dir}/result.tmp"
    echo "DEBUG: FINISHED $case_name (Exit code: $status)" >> "$DEBUG_LOG"
}

echo -e "${BLUE}🚀 Starting cg2all representation tests (Full simulations).${NC}"
echo -e "${YELLOW}Testing --cg2all-representation: 'calpha' and 'calpha-sc'${NC}"

# FLEX Cases
for rep in "calpha" "calpha-sc"; do
    case_name="FLEX_rep_${rep}"
    cmd="$CABSFLEX_CMD -i $FLEX_PDB:LH -A --aa-method cg2all --cg2all-representation $rep $OUTPUT_OPTS"
    run_test "$case_name" "$cmd" & limit_jobs
done

# DOCK Cases
for rep in "calpha" "calpha-sc"; do
    case_name="DOCK_rep_${rep}"
    cmd="$CABSDOCK_CMD -i $DOCK_PDB:AB -p KYVATLGV:EECTTTTC -p TAGQEKFGGLRDGYYI:CCCHHHTCCCCHHHHC -A --aa-method cg2all --cg2all-representation $rep $OUTPUT_OPTS"
    run_test "$case_name" "$cmd" & limit_jobs
done

# Wait for all background jobs to finish
wait

# Compile results
cat ${SCENARIO_ROOT}/*/result.tmp > "$RESULTS_FILE"

echo -e "\n${BLUE}========================================${NC}"
echo -e "${BLUE}   cg2all Representation Test Report    ${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "| Case | Result |"
echo -e "|------|--------|"

sort "$RESULTS_FILE" | while IFS='|' read -r case res; do
    if [ "$res" == "FAILED" ]; then
        echo -e "| $case | ${RED}$res${NC} |"
    else
        echo -e "| $case | ${GREEN}$res${NC} |"
    fi
done
echo -e "${BLUE}========================================${NC}"
