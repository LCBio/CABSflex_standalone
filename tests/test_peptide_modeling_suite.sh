#!/bin/bash
# ==============================================================================
# Peptide Modeling Test Suite for CABS-flex (de novo and restrained)
# ==============================================================================

# --- Configuration ---
set +e
GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

# Dynamically read ENV_NAME from install.sh
if [ -f "install.sh" ]; then
    ENV_NAME=$(grep '^ENV_NAME=' install.sh | cut -d'=' -f2 | tr -d '"' | tr -d "'")
else
    ENV_NAME="cabs"
fi

# Enable job control for 'jobs' command to work in script
set -m

# Micromamba activation support
export MAMBA_EXE="$HOME/.local/bin/micromamba"
if [ -f "$MAMBA_EXE" ]; then
    eval "$($MAMBA_EXE shell hook -s bash)"
    micromamba activate $ENV_NAME
else
    # Fallback to conda
    CONDA_BASE=$(conda info --base)
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate $ENV_NAME
fi

# Use PYTHONPATH=. python3 -m CABS to use the local fixed code
# Use CABSflex command via micromamba
RUN_CMD="micromamba run -p cabs CABSflex"

SCENARIO_ROOT="tests/test_peptide_modeling"
RESULTS_FILE="${SCENARIO_ROOT}/test_results.log"
rm -rf "$SCENARIO_ROOT"
mkdir -p "$SCENARIO_ROOT"

# Test parameters (Default CABS-flex simulation depth)
MC_ANNEALING=50
MC_CYCLES=100
MC_STEPS=100
COMMON_OPTS="-A --aa-method cg2all -v 4 --log"

# Concurrency limiting
if [[ "$OSTYPE" == "darwin"* ]]; then
    NPROC=$(sysctl -n hw.ncpu)
else
    NPROC=$(nproc)
fi
MAX_JOBS=$((NPROC - 1))
[ $MAX_JOBS -lt 1 ] && MAX_JOBS=1

limit_jobs() {
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 1
    done
}

# Verification helpers
check_all_atom() {
    local file="$1"
    if [ ! -f "$file" ]; then
        echo "FAILED (Not found)"
        return 1
    fi
    # Check for atoms other than CA and SC (N, O)
    local n_count=$(grep -c " N   " "$file")
    local o_count=$(grep -c " O   " "$file")
    if [ "$n_count" -gt 0 ] && [ "$o_count" -gt 0 ]; then
        echo "PASSED ($n_count N, $o_count O)"
        return 0
    else
        echo "FAILED (No AA details)"
        return 1
    fi
}

check_cyclic() {
    local file="$1"
    local n_res="$2"
    if [ ! -f "$file" ]; then
        echo "FAILED (Not found)"
        return 1
    fi
    # Extract CA coords of first and last residue
    local first_ca=$(awk '$3=="CA" {print $7","$8","$9; exit}' "$file")
    local last_ca=$(awk -v nr="$n_res" '$3=="CA" && $6==nr {print $7","$8","$9; exit}' "$file")
    if [ -z "$first_ca" ] || [ -z "$last_ca" ]; then
        echo "FAILED (Coords not found)"
        return 1
    fi
    # Calculate distance using python
    local dist=$(python3 -c "import math; f=[$first_ca]; l=[$last_ca]; print(math.sqrt(sum((f[i]-l[i])**2 for i in range(3))))")
    if (( $(echo "$dist < 5.0" | bc -l) )); then
         echo "PASSED (Dist: ${dist:0:4}A)"
         return 0
    else
         return 1
    fi
}

check_distance() {
    local file="$1"
    local r1="$2"
    local r2="$3"
    local max_dist="$4"
    if [ ! -f "$file" ]; then
        echo "FAILED (Not found)"
        return 1
    fi
    local c1=$(awk -v r="$r1" '$3=="CA" && $6==r {print $7","$8","$9; exit}' "$file")
    local c2=$(awk -v r="$r2" '$3=="CA" && $6==r {print $7","$8","$9; exit}' "$file")
    if [ -z "$c1" ] || [ -z "$c2" ]; then
        echo "FAILED (Coords $r1/$r2 not found)"
        return 1
    fi
    local dist=$(python3 -c "import math; f=[$c1]; l=[$c2]; print(math.sqrt(sum((f[i]-l[i])**2 for i in range(3))))")
    if (( $(echo "$dist < $max_dist" | bc -l) )); then
         echo "PASSED (Dist: ${dist:0:4}A)"
         return 0
    else
         echo "FAILED (Dist: ${dist:0:4}A)"
         return 1
    fi
}

run_test() {
    local case_name="$1"
    local cmd="$2"
    local verify_expr="$3"
    local out_dir="${SCENARIO_ROOT}/${case_name}"
    mkdir -p "$out_dir"
    
    echo -e "${BLUE}Starting: $case_name${NC}"
    eval "$cmd -w $out_dir -a $MC_ANNEALING -y $MC_CYCLES -s $MC_STEPS > ${out_dir}/stdout.log 2> ${out_dir}/stderr.log"
    local status=$?
    
    local verification="N/A"
    if [ $status -eq 0 ] && [ -z "$verify_expr" ]; then
        verification="PASSED"
    elif [ $status -eq 0 ]; then
        # Replace OUT_DIR with actual path in verify_expr
        local actual_verify_expr=$(echo "$verify_expr" | sed "s|OUT_DIR|$out_dir|g")
        verification=$(eval "$actual_verify_expr")
        [ $? -ne 0 ] && status=1
    fi

    local result="PASSED"
    [ $status -ne 0 ] && result="FAILED"
    
    echo "$case_name|$result|$verification" >> "$RESULTS_FILE"
    echo -e "Finished $case_name: ${result}"
}

echo -e "${YELLOW}🚀 Starting Peptide Modeling Test Suite${NC}"

# Case 1: Linear Peptide de novo (Example from test_cabsdock_options.sh)
run_test "PEPTIDE_LINEAR" "$RUN_CMD --peptide EEWEVLCWTWETCER:CCCEEEECCCTTCCC $COMMON_OPTS" \
"check_all_atom OUT_DIR/output_pdbs/model_0.pdb" & limit_jobs

# Case 2: Peptide with Restraints (Cyclic-like)
run_test "PEPTIDE_RESTRAINED" "$RUN_CMD --peptide KYVATLGV:EECTTTTC --ca-rest-add 1:PEP 8:PEP 3.8 1.0 $COMMON_OPTS" \
"check_all_atom OUT_DIR/output_pdbs/model_0.pdb" & limit_jobs

# Case 3: Multiple Peptides de novo
run_test "PEPTIDE_MULTI" "$RUN_CMD --peptide KYVATLGV:EECTTTTC --peptide TAGQEKFGGLRDGYYI:CCCHHHTCCCCHHHHC $COMMON_OPTS" \
"check_all_atom OUT_DIR/output_pdbs/model_0.pdb" & limit_jobs

# Case 5: Cyclic Peptide (Backbone Cyclization)
run_test "PEPTIDE_CYCLIC" "$RUN_CMD --peptide EEWEVLCWTWETCER:CCCEEEECCCTTCCC --backbone-cyclization PEP1 $COMMON_OPTS" \
"check_cyclic OUT_DIR/output_pdbs/model_0.pdb 15" & limit_jobs

# Case 6: Disulfide Peptide (CEWEVLCWTWETCEC - Cys at 1 and 15)
run_test "PEPTIDE_DISULFIDE" "$RUN_CMD --peptide CEWEVLCWTWETCEC:CCCEEEECCCTTCCC --disulfide-bonds 1:PEP1 15:PEP1 $COMMON_OPTS" \
"check_distance OUT_DIR/output_pdbs/model_0.pdb 1 15 6.0" & limit_jobs

# Case 7: Cyclic + Disulfide (Combined backbone and sidechain cyclization)
run_test "PEPTIDE_HYBRID" "$RUN_CMD --peptide CEWEVLCWTWETCEC:CCCEEEECCCTTCCC --backbone-cyclization PEP1 --disulfide-bonds 1:PEP1 15:PEP1 $COMMON_OPTS" \
"check_distance OUT_DIR/output_pdbs/model_0.pdb 1 15 5.0" & limit_jobs

# Case 8: Peptide without Secondary Structure (Linear)
run_test "PEPTIDE_NO_SS" "$RUN_CMD --peptide EEWEVLCWTWETCER $COMMON_OPTS" \
"check_all_atom OUT_DIR/output_pdbs/model_0.pdb" & limit_jobs

# Case 9: Multi-Peptide without Secondary Structure
run_test "PEPTIDE_MULTI_NO_SS" "$RUN_CMD --peptide KYVATLGV --peptide TAGQEKFGGLRDGYYI $COMMON_OPTS" \
"check_all_atom OUT_DIR/output_pdbs/model_0.pdb" & limit_jobs

# Case 10: Impossible Disulfide (ALA at 1 and 15) - Correct failure is PASS
# We wrap it in a subshell to invert the status
(
    out_dir="${SCENARIO_ROOT}/PEPTIDE_DISULFIDE_FAIL"
    mkdir -p "$out_dir"
    echo -e "${BLUE}Starting: PEPTIDE_DISULFIDE_FAIL (Expect Failure)${NC}"
    eval "$RUN_CMD --peptide AEWEVLCWTWETCEA:CCCEEEECCCTTCCC --disulfide-bonds 1:PEP1 15:PEP1 $COMMON_OPTS -w $out_dir -a $MC_ANNEALING -y $MC_CYCLES -s $MC_STEPS > ${out_dir}/stdout.log 2> ${out_dir}/stderr.log"
    if [ $? -ne 0 ]; then
        echo "PEPTIDE_DISULFIDE_FAIL|PASSED|Failed as expected" >> "$RESULTS_FILE"
        echo -e "Finished PEPTIDE_DISULFIDE_FAIL: PASSED"
    else
        echo "PEPTIDE_DISULFIDE_FAIL|FAILED|Did not fail on impossible disulfide" >> "$RESULTS_FILE"
        echo -e "Finished PEPTIDE_DISULFIDE_FAIL: FAILED"
    fi
) & limit_jobs

wait

echo -e "\n${BLUE}========================================${NC}"
echo -e "${BLUE}Peptide Modeling Test Report${NC}"
echo -e "${BLUE}========================================${NC}"
printf "| %-20s | %-10s | %-20s |\n" "Scenario" "Result" "Verification"
echo "|----------------------|------------|----------------------|"
while IFS='|' read -r name res ver; do
    if [ "$res" == "FAILED" ]; then
        printf "| %-20s | ${RED}%-10s${NC} | %-20s |\n" "$name" "$res" "$ver"
    else
        printf "| %-20s | ${GREEN}%-10s${NC} | %-20s |\n" "$name" "$res" "$ver"
    fi
done < "$RESULTS_FILE"
echo -e "${BLUE}========================================${NC}"
