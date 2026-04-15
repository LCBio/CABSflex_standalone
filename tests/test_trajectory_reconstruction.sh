#!/bin/bash
export PYTHONPATH=$(pwd):$PYTHONPATH
# ==============================================================================
#           Trajectory Reconstruction Tests for CABSflex
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

# Define local commands to ensure workspace code is used
CABSFLEX_CMD="micromamba run -n $ENV_NAME CABSflex"
CABSDOCK_CMD="micromamba run -n $ENV_NAME CABSdock"

SCENARIO_ROOT="tests/test_trajectory_reconstruction"
RESULTS_FILE="tests/test_trajectory_reconstruction_results.log"

rm -rf "$SCENARIO_ROOT" "$RESULTS_FILE"
mkdir -p "$SCENARIO_ROOT"

# Robust simulation options for testing reconstruction
TEST_OPTS="-y 10 -s 50" # robust cycles, default replicas
FLEX_PDB="tests/inputs/2BZ6.pdb"
BASE_OPTS="--aa-method cg2all --cg2all-representation calpha"

if [ ! -f "$FLEX_PDB" ]; then
    echo -e "${RED}Error: Input file $FLEX_PDB not found.${NC}"
    exit 1
fi

run_test() {
    local case_name="$1"
    local cmd="$2"
    local out_dir="${SCENARIO_ROOT}/${case_name}"
    
    mkdir -p "$out_dir"
    
    echo "Running: $case_name"
    eval "$cmd -w $out_dir > ${out_dir}/stdout.log 2> ${out_dir}/stderr.log"
    local status=$?
    
    local result="PASSED"
    [ $status -ne 0 ] && result="FAILED"
    
    # Specific checks for each case
    if [ "$result" == "PASSED" ]; then
        case "$case_name" in
            "BASIC")
                # Should have replica.pdb (now all-atom)
                if [ ! -f "${out_dir}/output_pdbs/replica.pdb" ]; then
                    result="FAILED (No output trajectory PDB)"
                else
                    # Verify it's all-atom (look for atoms other than CA)
                    if ! grep -q " N   " "${out_dir}/output_pdbs/replica.pdb"; then
                         result="FAILED (Trajectory is not all-atom)"
                    fi
                fi
                ;;
            "SAMPLING")
                if [ ! -f "${out_dir}/output_pdbs/replica.pdb" ]; then
                    result="FAILED (No output trajectory PDB)"
                else
                    # Verify model count (from FAST_OPTS: -y 10 with sample 5)
                    local model_count=$(grep -c "MODEL" "${out_dir}/output_pdbs/replica.pdb")
                    if [ "$model_count" -ne 2 ]; then
                        result="FAILED (Wrong model count: $model_count, expected 2)"
                    fi
                fi
                ;;
            "PARALLEL")
                if [ ! -f "${out_dir}/output_pdbs/replica_0.pdb" ] || [ ! -f "${out_dir}/output_pdbs/replica_1.pdb" ]; then
                    result="FAILED (No output trajectory PDBs)"
                fi
                # Check log for parallel mention
                if ! grep -q "Starting parallel reconstruction using" "${out_dir}/output_data/CABSflex.log" 2>/dev/null; then
                    result="FAILED (No parallel mention in log)"
                fi
                ;;
            "EARLY_VALIDATION_FAIL")
                # This one SHOULD fail
                result="FAILED (Should have failed early)"
                ;;
        esac
    else
        # If it failed, check if it was expected (EARLY_VALIDATION_FAIL)
        if [ "$case_name" == "EARLY_VALIDATION_FAIL" ]; then
            if grep -q "Failed to verify cg2all environment" "${out_dir}/stderr.log" || grep -q "Failed to verify cg2all environment" "${out_dir}/stdout.log"; then
                result="PASSED (Correctly failed early)"
            else
                result="FAILED (Failed but for wrong reason)"
            fi
        fi
    fi

    echo "$case_name|$result" > "${out_dir}/result.tmp"
}

echo -e "${BLUE}🚀 Starting Trajectory Reconstruction Tests${NC}"

# 1. BASIC
run_test "BASIC" "CABSflex -i $FLEX_PDB:LH $TEST_OPTS $BASE_OPTS --aa-rebuild-trajectory"

# 2. SAMPLING
run_test "SAMPLING" "CABSflex -i $FLEX_PDB:LH -a 1 -y 10 -s 10 --aa-rebuild-trajectory --aa-rebuild-trajectory-sample 5"

# 3. PARALLEL (Using CABSdock default of 10 replicas)
run_test "PARALLEL" "CABSdock -i $FLEX_PDB:LH -y 10 -s 10 $BASE_OPTS --aa-rebuild-trajectory --aa-rebuild-trajectory-parallel 4"

# 4. EARLY_VALIDATION_FAIL
run_test "EARLY_VALIDATION_FAIL" "CABSflex -i $FLEX_PDB:LH $TEST_OPTS $BASE_OPTS --aa-rebuild-trajectory --cg2all-env-prefix /tmp/nonexistent_env"

# Compile results
cat ${SCENARIO_ROOT}/*/result.tmp > "$RESULTS_FILE"

echo -e "\n${BLUE}========================================${NC}"
echo -e "${BLUE}   Trajectory Reconstruction Report      ${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "| Case | Result |"
echo -e "|------|--------|"

sort "$RESULTS_FILE" | while IFS='|' read -r case res; do
    if [[ "$res" == "FAILED"* ]]; then
        echo -e "| $case | ${RED}$res${NC} |"
    else
        echo -e "| $case | ${GREEN}$res${NC} |"
    fi
done
echo -e "${BLUE}========================================${NC}"
