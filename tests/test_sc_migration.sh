#!/bin/bash
# ==============================================================================
#           Side-Chain Center (SC) Migration Tests for CABSflex/dock
# ==============================================================================

# --- Configuration ---
set +e
GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

ENV_NAME="cabs"
SCENARIO_ROOT="tests/test_sc_migration"
RESULTS_FILE="tests/test_sc_migration_results.log"
DEBUG_LOG="${SCENARIO_ROOT}/debug_progress.log"

rm -rf "$SCENARIO_ROOT" "$RESULTS_FILE"
mkdir -p "$SCENARIO_ROOT"
touch "$DEBUG_LOG"

# Inputs
FLEX_PDB="tests/inputs/2BZ6.pdb"
DOCK_PDB="tests/inputs/1A2K.pdb"
DOCK_REC="$DOCK_PDB:AB"

# Micromamba command prefix
RUN_CMD="micromamba run -n $ENV_NAME"

# Enable job control for 'jobs' command
set -m

# Concurrency limiting logic
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

# --- Verification Helpers ---

check_pdb_sc() {
    local file="$1"
    local expected_sc="$2" # true or false
    
    # If looking for replica_0.pdb but it doesn't exist, try replica.pdb (for CABSflex)
    if [[ "$file" == *"replica_0.pdb" ]] && [ ! -f "$file" ]; then
        local alt_file="${file%_0.pdb}.pdb"
        if [ -f "$alt_file" ]; then
            file="$alt_file"
        fi
    fi

    if [ ! -f "$file" ]; then
        echo "FAILED (File not found: $file)"
        return 1
    fi
    
    local sc_count=$(grep -c " SC " "$file")
    local ca_count=$(grep -c " CA " "$file")
    
    if [ "$expected_sc" == "true" ]; then
        if [ "$sc_count" -gt 0 ]; then
            # More precise check: SC + GLY == CA
            local gly_count=$(grep " CA " "$file" | grep -c "GLY")
            if [ $((sc_count + gly_count)) -eq "$ca_count" ]; then
                echo "PASSED (SC atoms present and interleaved correctly: $sc_count SC, $ca_count CA)"
                return 0
            else
                echo "FAILED (SC/CA ratio mismatch: $sc_count SC, $ca_count CA, $gly_count GLY)"
                return 1
            fi
        else
            echo "FAILED (No SC atoms found but expected)"
            return 1
        fi
    else
        if [ "$sc_count" -eq 0 ]; then
            echo "PASSED (No SC atoms found as expected)"
            return 0
        else
            echo "FAILED ($sc_count SC atoms found but expected NONE)"
            return 1
        fi
    fi
}

run_test() {
    local case_name="$1"
    local cmd="$2"
    local verify_func="$3"
    local out_dir="${SCENARIO_ROOT}/${case_name}"
    
    mkdir -p "$out_dir"
    echo -e "${BLUE}Running: $case_name${NC}"
    
    # Execute (using default simulation parameters)
    eval "$cmd -w $out_dir --save-config -v 4 -S > ${out_dir}/stdout.log 2> ${out_dir}/stderr.log"
    local status=$?
    
    local result="PASSED"
    local details=""
    if [ $status -ne 0 ]; then
        result="FAILED"
        details="Execution failed with exit code $status"
    else
        # Run verification
        details=$(eval "$verify_func $out_dir")
        if [[ "$details" == "FAILED"* ]]; then
            result="FAILED"
        fi
    fi
    
    echo "$case_name|$result|$details" > "${out_dir}/result.tmp"
    echo "DEBUG: FINISHED $case_name ($result)" >> "$DEBUG_LOG"
}

# --- Test Cases ---

echo -e "${BLUE}🚀 Starting Side-Chain Center (SC) Migration Tests.${NC}"

# 1. CABSflex Default (Should have SC in coarse-grained replica)
run_test "FLEX_Default" "$RUN_CMD CABSflex -i $FLEX_PDB:LH" \
"check_pdb_sc \$1/output_pdbs/replica_0.pdb true" & limit_jobs

# 2. CABSflex Disabled SC (Should NOT have SC)
run_test "FLEX_Disabled" "$RUN_CMD CABSflex -i $FLEX_PDB:LH --disable-side-chain-centers" \
"check_pdb_sc \$1/output_pdbs/replica_0.pdb false" & limit_jobs

# 3. CABSflex Optional Start SC (Should have SC in start.pdb)
run_test "FLEX_StartSC" "$RUN_CMD CABSflex -i $FLEX_PDB:LH --write-sc-start-pdbs" \
"check_pdb_sc \$1/output_pdbs/start.pdb true" & limit_jobs

# 4. CABSdock Default (Multiple Peptides, Should have SC)
run_test "DOCK_Default" "$RUN_CMD CABSdock -i $DOCK_REC -p KYVATLGV:EECTTTTC -p TAGQEKFGGLRDGYYI:CCCHHHTCCCCHHHHC" \
"check_pdb_sc \$1/output_pdbs/replica_0.pdb true" & limit_jobs

# 5. CABSdock Single Peptide A (Should have SC)
run_test "DOCK_PepA" "$RUN_CMD CABSdock -i $DOCK_REC -p KYVATLGV:EECTTTTC" \
"check_pdb_sc \$1/output_pdbs/replica_0.pdb true" & limit_jobs

# 6. CABSdock Single Peptide B (Should have SC)
run_test "DOCK_PepB" "$RUN_CMD CABSdock -i $DOCK_REC -p TAGQEKFGGLRDGYYI:CCCHHHTCCCCHHHHC" \
"check_pdb_sc \$1/output_pdbs/replica_0.pdb true" & limit_jobs

# 7. cg2all Reconstruction (calpha-sc) - Check replica for SC
run_test "CG2ALL_SC" "$RUN_CMD CABSflex -i $FLEX_PDB:LH -A --aa-method cg2all --cg2all-representation calpha-sc" \
"check_pdb_sc \$1/output_pdbs/replica_0.pdb true" & limit_jobs

# 8. cg2all Reconstruction (calpha) - SC should be AUTO-DISABLED for replica too
run_test "CG2ALL_CA" "$RUN_CMD CABSflex -i $FLEX_PDB:LH -A --aa-method cg2all --cg2all-representation calpha" \
"check_pdb_sc \$1/output_pdbs/replica_0.pdb false" & limit_jobs

# 9. Overriding calpha with SC centers enabled
run_test "CG2ALL_CA_SC_Override" "$RUN_CMD CABSflex -i $FLEX_PDB:LH -A --aa-method cg2all --cg2all-representation calpha --disable-side-chain-centers false" \
"check_pdb_sc \$1/output_pdbs/replica_0.pdb true" & limit_jobs

wait

# Special check for start.pdb in FLEX_Default (Should NOT have SC by default)
details=$(check_pdb_sc "${SCENARIO_ROOT}/FLEX_Default/output_pdbs/start.pdb" false)
if [[ "$details" == "FAILED"* ]]; then
    echo "FLEX_Default_Start|FAILED|$details" > "${SCENARIO_ROOT}/FLEX_Default_Start.tmp"
else
    echo "FLEX_Default_Start|PASSED|$details" > "${SCENARIO_ROOT}/FLEX_Default_Start.tmp"
fi

# Compile results
cat ${SCENARIO_ROOT}/*/result.tmp ${SCENARIO_ROOT}/*.tmp 2>/dev/null > "$RESULTS_FILE"

# Report
echo -e "\n${BLUE}========================================${NC}"
echo -e "${BLUE}   SC Migration Test Report             ${NC}"
echo -e "${BLUE}========================================${NC}"
printf "| %-25s | %-10s | %-40s |\n" "Case" "Result" "Details"
echo -e "|---------------------------|------------|------------------------------------------|"

sort "$RESULTS_FILE" | while IFS='|' read -r case res det; do
    if [ "$res" == "FAILED" ]; then
        printf "| %-25s | ${RED}%-10s${NC} | %-40s |\n" "$case" "$res" "$det"
    else
        printf "| %-25s | ${GREEN}%-10s${NC} | %-40s |\n" "$case" "$res" "$det"
    fi
done
echo -e "${BLUE}========================================${NC}"
