#!/bin/bash

# Comprehensive Test Suite for CABSflex
# This script runs CABSflex with various inputs and options.

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

ENV_NAME="cabs"
CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate $ENV_NAME

# Check if CABSflex is available
if ! command -v CABSflex &> /dev/null; then
    echo -e "${RED}❌ CABSflex command not found in environment '${ENV_NAME}'.${NC}"
    exit 1
fi

echo -e "${BLUE}🚀 Starting CABSflex Comprehensive Test Suite${NC}"

# Test Directories
SCENARIO_ROOT="tests/scenarios"
mkdir -p "$SCENARIO_ROOT"

# Function to run a test scenario
run_test() {
    local name=$1
    local cmd=$2
    local output_dir=$3
    
    echo -e "\n${YELLOW}🧪 Testing Scenario: ${name}${NC}"
    echo -e "${BLUE}Running: ${cmd}${NC}"
    
    # Remove output dir if exists
    rm -rf "$output_dir"
    
    # Run command
    eval "$cmd"
    local status=$?
    
    # Check if any all_rmsds_*.txt file exists in output_data
    if [ $status -eq 0 ] && [ -d "$output_dir" ] && ( ls "${output_dir}/output_data/all_rmsds_"*.txt &> /dev/null || [ -f "${output_dir}/CABS.log" ] ); then
        echo -e "${GREEN}✅ Scenario '${name}' PASSED${NC}"
        return 0
    else
        echo -e "${RED}❌ Scenario '${name}' FAILED (Exit code: ${status})${NC}"
        return 1
    fi
}

FAILED_TESTS=0

# Scenario 1: Helix (Single Chain PDB)
run_test "Helix PDB" \
    "CABSflex -i tests/inputs/Helix.pdb -a 2 -y 5 -s 2 -w ${SCENARIO_ROOT}/helix -v 2" \
    "${SCENARIO_ROOT}/helix" || ((FAILED_TESTS++))

# Scenario 2: 8DTQ (Multi-chain mmCIF)
run_test "8DTQ mmCIF" \
    "CABSflex -i tests/inputs/8DTQ.cif -a 2 -y 5 -s 2 -w ${SCENARIO_ROOT}/8dtq -v 2" \
    "${SCENARIO_ROOT}/8dtq" || ((FAILED_TESTS++))

# Scenario 3: 5C2N_1 (Chain selection)
run_test "5C2N_1 Chain Selection" \
    "CABSflex -i tests/inputs/5C2N_1.pdb:BDFHJ -a 2 -y 5 -s 2 -w ${SCENARIO_ROOT}/5c2n_1 -v 2" \
    "${SCENARIO_ROOT}/5c2n_1" || ((FAILED_TESTS++))

# Scenario 4: Direct Fetching 5C2N:A
run_test "Fetch 5C2N:A" \
    "CABSflex -i 5C2N:A -a 2 -y 5 -s 2 -w ${SCENARIO_ROOT}/fetch_5c2n -v 2" \
    "${SCENARIO_ROOT}/fetch_5c2n" || ((FAILED_TESTS++))

# Scenario 5: Direct Fetching ModelArchive-like ID
run_test "Fetch 8zgr:LY" \
    "CABSflex -i pdb_00008zgr:LY -a 2 -y 5 -s 2 -w ${SCENARIO_ROOT}/fetch_8zgr -v 2" \
    "${SCENARIO_ROOT}/fetch_8zgr" || ((FAILED_TESTS++))

# Scenario 6: Helix Comprehensive (All Args)
run_test "Helix Comprehensive Args" \
    "CABSflex -i tests/inputs/Helix_short.pdb -a 3 -y 10 -s 5 -r 2 -D 0.5 -t 2.0 2.0 -k 5 -R tests/inputs/Helix_short.pdb -S -C -M -v 3 --dssp-output --ss-output --json-output --pdb-bfac-output BRS --csv-output BCPS --image-file-format png -w ${SCENARIO_ROOT}/helix_comprehensive" \
    "${SCENARIO_ROOT}/helix_comprehensive" || ((FAILED_TESTS++))

# Scenario 7: Flexibility and Restraints
run_test "Flexibility and Restraints" \
    "CABSflex -i tests/inputs/Helix_short.pdb -f tests/inputs/flex.txt -g plddt 3 3.8 11.5 --ca-rest-add 1:A 17:A 25.0 1.0 --ca-rest-file tests/inputs/restraints.txt --ca-rest-weight 1.0 1.0 --sc-rest-weight 1.0 1.0 -F 1:A 17:A --backbone-cyclization A --receptor-ss tests/inputs/ss_assignment.txt --protein-restraints-retain 80 --renumber-residues-to-original -a 2 -y 5 -s 2 -w ${SCENARIO_ROOT}/flex_restraints -v 2" \
    "${SCENARIO_ROOT}/flex_restraints" || ((FAILED_TESTS++))

# Scenario 8: Simulation and Analysis Focus
run_test "Simulation and Analysis" \
    "CABSflex -i tests/inputs/Helix_short.pdb -z 12345 -b 1.5 --disable-centro --clustering-iterations 10 --filtering-mode all -n 20 -T 8.0 --contact-threshold-aa 6.0 --contact-map-colors '#ff0000' '#00ff00' '#0000ff' '#ffff00' '#ff00ff' '#00ffff' -a 2 -y 5 -s 2 -w ${SCENARIO_ROOT}/simulation_analysis -v 2" \
    "${SCENARIO_ROOT}/simulation_analysis" || ((FAILED_TESTS++))

# Scenario 9: Config file
run_test "Config File" \
    "CABSflex -c ${SCENARIO_ROOT}/helix_comprehensive/config.ini -w ${SCENARIO_ROOT}/config_file -v 2" \
    "${SCENARIO_ROOT}/config_file" || ((FAILED_TESTS++))

# Scenario 10: Loading from previous run
# We use the output from Scenario 6
run_test "Load Files" \
    "CABSflex -L ${SCENARIO_ROOT}/helix_comprehensive/*.cbs -v 2 -R tests/inputs/Helix_short.pdb -w ${SCENARIO_ROOT}/load_files" \
    "${SCENARIO_ROOT}/load_files" || ((FAILED_TESTS++))

# Summary
echo -e "\n${BLUE}========================================${NC}"
if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "${GREEN}⭐ ALL TESTS PASSED${NC}"
    exit 0
else
    echo -e "${RED}💥 ${FAILED_TESTS} TEST(S) FAILED${NC}"
    exit 1
fi
