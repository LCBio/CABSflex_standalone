#!/bin/bash
# ==============================================================================
# Comprehensive Test Suite for CABSflex - Exhaustive CLI Coverage
# ==============================================================================
set -e

# --- Configuration ---
GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'
ENV_NAME="cabs"
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate $ENV_NAME

# Check environment
if ! command -v CABSflex &> /dev/null; then
    echo -e "${RED}❌ CABSflex command not found in environment '${ENV_NAME}'.${NC}"
    exit 1
fi

echo -e "${BLUE}🚀 Starting CABSflex Comprehensive Test Suite (Exhaustive CLI Check)${NC}"

SCENARIO_ROOT="tests/scenarios"
mkdir -p "$SCENARIO_ROOT"

# Function to run a test scenario
run_test() {
    local name=$1
    local cmd=$2
    local output_dir=$3
    local option_tested=$4
    
    echo -e "\n${YELLOW}🧪 Testing Scenario: ${name}${NC}"
    echo -e "${BLUE}Running: ${cmd}${NC}"
    
    rm -rf "$output_dir"
    
    eval "$cmd"
    local status=$?
    
    # Success check: Must exit 0 AND produce log OR output_pdbs
    if [ $status -eq 0 ] && ( ls "${output_dir}/CABS.log" &> /dev/null || ls "${output_dir}/output_pdbs" &> /dev/null ); then
        echo -e "${GREEN}✅ Scenario '${name}' PASSED${NC}"
        return 0
    else
        echo -e "${RED}❌ Scenario '${name}' FAILED.${NC}"
        echo -e "${RED}   -> Option/Feature Tested: ${option_tested}${NC}"
        echo -e "${RED}   -> Exit Code: ${status}${NC}"
        return 1
    fi
}

FAILED_TESTS=0
TEST_INPUT="-i tests/inputs/Helix.pdb" # Main test structure
RECPTOR_SS="A:CCCHHHHHHHHHHHCCC"

# --- SCENARIOS: Exhaustive Coverage of All Groups ---

# 2.1 BASIC OPTIONS
run_test "Basic Defaults" \
    "CABSflex $TEST_INPUT -w ${SCENARIO_ROOT}/basic_def -a 2 -y 5 -s 2" \
    "${SCENARIO_ROOT}/basic_def" "-i + Defaults" || ((FAILED_TESTS++))

run_test "Config File Override" \
    "CABSflex -c tests/inputs/config.ini -w ${SCENARIO_ROOT}/config_load -v 2 -i tests/inputs/Helix_short.pdb" \
    "${SCENARIO_ROOT}/config_load" "--config" || ((FAILED_TESTS++))

# 2.2 PROTEIN STRUCTURE OPTIONS
run_test "Protein Flexibility (File) -f" \
    "CABSflex $TEST_INPUT -f tests/inputs/flex.txt -w ${SCENARIO_ROOT}/flex_file" \
    "${SCENARIO_ROOT}/flex_file" "-f <file>" || ((FAILED_TESTS++))

run_test "Protein Flexibility (File) --protein-flexibility" \
    "CABSflex $TEST_INPUT --protein-flexibility tests/inputs/flex.txt -w ${SCENARIO_ROOT}/flex_file" \
    "${SCENARIO_ROOT}/flex_file" "--protein-flexibility <file>" || ((FAILED_TESTS++))

run_test "Protein Restraints (SS2) -g" \
    "CABSflex $TEST_INPUT -g ss2 3 3.8 8.0 -w ${SCENARIO_ROOT}/prot_rest_ss2" \
    "${SCENARIO_ROOT}/prot_rest_ss2" "-g" || ((FAILED_TESTS++))

run_test "Protein Restraints (SS2) --protein-restraints" \
    "CABSflex $TEST_INPUT --protein-restraints ss2 3 3.8 8.0 -w ${SCENARIO_ROOT}/prot_rest_ss2" \
    "${SCENARIO_ROOT}/prot_rest_ss2" "--protein-restraints" || ((FAILED_TESTS++))

run_test "No Protein Restraints -N" \
    "CABSflex $TEST_INPUT -N -w ${SCENARIO_ROOT}/prot_no_rest" \
    "${SCENARIO_ROOT}/prot_no_rest" "-N" || ((FAILED_TESTS++))

run_test "No Protein Restraints --no-protein-restraints" \
    "CABSflex $TEST_INPUT --no-protein-restraints -w ${SCENARIO_ROOT}/prot_no_rest" \
    "${SCENARIO_ROOT}/prot_no_rest" "--no-protein-restraints" || ((FAILED_TESTS++))
    
run_test "Protein Restraints Retain" \
    "CABSflex $TEST_INPUT --protein-restraints-retain 80 -w ${SCENARIO_ROOT}/prot_rest_retain" \
    "${SCENARIO_ROOT}/prot_rest_retain" "--protein-restraints-retain" || ((FAILED_TESTS++))

run_test "Receptor SS Assignment --receptor-ss" \
    "CABSflex $TEST_INPUT --receptor-ss tests/inputs/ss_assignment.txt -w ${SCENARIO_ROOT}/receptor_ss" \
    "${SCENARIO_ROOT}/receptor_ss" "--receptor-ss"  ${RECPTOR_SS}|| ((FAILED_TESTS++))
    
run_test "Weighted Fit (Gauss)" \
    "CABSflex $TEST_INPUT --weighted-fit gauss --gauss-iterations 5 -w ${SCENARIO_ROOT}/weight_gauss" \
    "${SCENARIO_ROOT}/weight_gauss" "--weighted-fit gauss" || ((FAILED_TESTS++))

run_test "Weighted Fit (flex)" \
    "CABSflex $TEST_INPUT --weighted-fit flex -w ${SCENARIO_ROOT}/weight_flex" \
    "${SCENARIO_ROOT}/weight_flex" "--weighted-fit flex" || ((FAILED_TESTS++))

run_test "Weighted Fit (ss)" \
    "CABSflex $TEST_INPUT --weighted-fit ss -w ${SCENARIO_ROOT}/weight_ss" \
    "${SCENARIO_ROOT}/weight_ss" "--weighted-fit ss" || ((FAILED_TESTS++))

run_test "Weighted Fit (off)" \
    "CABSflex $TEST_INPUT --weighted-fit off -w ${SCENARIO_ROOT}/weight_off" \
    "${SCENARIO_ROOT}/weight_off" "--weighted-fit off" || ((FAILED_TESTS++))

# 2.3 RESTRAINTS OPTIONS (Need to assume minimal files exist)
#run_test "Restraints (CA/SC/File/Disulfide)" \
#    "CABSflex $TEST_INPUT --ca-rest-add 1:A 17:A 25.0 1.0 --sc-rest-add 1:A 17:A 5.0 1.0 --ca-rest-file tests/inputs/restraints.txt --sc-rest-file tests/inputs/restraints.txt -F 1:A 17:A --backbone-cyclization A -w ${SCENARIO_ROOT}/all_rests" \
#    "${SCENARIO_ROOT}/all_rests" "All Restraint Flags" || ((FAILED_TESTS++))

# 2.4 SIMULATION OPTIONS
#run_test "Sim Flags (T, Cycles, Seed)" \
#    "CABSflex $TEST_INPUT -a 10 -y 100 -s 10 -r 2 -D 0.1 -t 3.0 1.0 -z 12345 -b 1.5 --disable-centro -w ${SCENARIO_ROOT}/sim_flags" \
#    "${SCENARIO_ROOT}/sim_flags" "All Sim Flags" || ((FAILED_TESTS++))

# 2.5 ALL-ATOM RECONSTRUCTION OPTIONS
run_test "AA Rebuild (Modeller)" \
    "CABSflex $TEST_INPUT --aa-rebuild --aa-method modeller --modeller-iterations 1 -w ${SCENARIO_ROOT}/aa_modeller" \
    "${SCENARIO_ROOT}/aa_modeller" "--aa-rebuild + -m" || ((FAILED_TESTS++))

run_test "AA Rebuild (cg2all)" \
    "CABSflex $TEST_INPUT --aa-rebuild --aa-method cg2all -w ${SCENARIO_ROOT}/aa_cg2all" \
    "${SCENARIO_ROOT}/aa_cg2all" "--aa-rebuild + -m" || ((FAILED_TESTS++))

run_test "AA Rebuild (default)" \
    "CABSflex $TEST_INPUT --aa-rebuild -w ${SCENARIO_ROOT}/aa_default" \
    "${SCENARIO_ROOT}/aa_default" "--aa-rebuild" || ((FAILED_TESTS++))

# 2.6 ANALYSIS OPTIONS
run_test "Analysis (Clustering & Maps)" \
    "CABSflex $TEST_INPUT -k 3 --filtering-count 50 --contact-maps -T 8.0 -R tests/inputs/Helix.pdb --align SW -w ${SCENARIO_ROOT}/analysis" \
    "${SCENARIO_ROOT}/analysis" "-k, -M, -R, --align" || ((FAILED_TESTS++))

# 2.7 OUTPUT OPTIONS
run_test "Output Flags" \
    "CABSflex $TEST_INPUT -S -C -o M --pdb-bfac-output BPS --csv-output S --json-output -w ${SCENARIO_ROOT}/output_flags" \
    "${SCENARIO_ROOT}/output_flags" "-S, -C, -o, --pdb-bfac-output, --csv-output, --json-output" || ((FAILED_TESTS++))

# 2.8 MISCELLANEOUS OPTIONS
run_test "Misc Options" \
    "CABSflex $TEST_INPUT -w ${SCENARIO_ROOT}/misc --image-file-format png --log --version" \
    "${SCENARIO_ROOT}/misc" "--work-dir + --image-format + --log" || ((FAILED_TESTS++))

# Scenario 9: Config File Override (Uses tests/scenarios/config_load/cabs_defaults.ini)
run_test "9. Config File Override" \
    "CABSflex $TEST_INPUT -c ${SCENARIO_ROOT}/config_load/cabs_defaults.ini -w ${SCENARIO_ROOT}/config_load -v 2" \
    "${SCENARIO_ROOT}/config_load" "--config" || ((FAILED_TESTS++))

# Scenario 10: Loading Previous Run (Requires a dummy .cbs file)
run_test "10. Load Files" \
    "CABSflex -L ${SCENARIO_ROOT}/helix_def/CABS.cbs -v 2 -R tests/inputs/Helix.pdb -w ${SCENARIO_ROOT}/load_files" \
    "${SCENARIO_ROOT}/load_files" "--load-cabs-files + -R" || ((FAILED_TESTS++))


# --- Final Summary ---
echo -e "\n${BLUE}========================================${NC}"
if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "${GREEN}⭐ ALL TESTS PASSED!${NC}"
    exit 0
else
    echo -e "${RED}💥 ${FAILED_TESTS} SCENARIO(S) FAILED. Review the output above.${NC}"
    exit 1
fi