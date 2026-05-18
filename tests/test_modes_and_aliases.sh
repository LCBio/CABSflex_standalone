#!/bin/bash

# Test Script for Legacy and Modern Restraint Modes & Aliases in CABS-flex
# This script runs quick, short simulations to verify all options function correctly.

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Dynamically read ENV_NAME from install.sh
if [ -f "install.sh" ]; then
    ENV_NAME=$(grep '^ENV_NAME=' install.sh | cut -d'=' -f2 | tr -d '"' | tr -d "'")
else
    ENV_NAME="cabs"
fi

# Enable micromamba/conda environment
export MAMBA_EXE="$HOME/.local/bin/micromamba"
if [ -f "$MAMBA_EXE" ]; then
    eval "$($MAMBA_EXE shell hook -s bash)"
    micromamba activate $ENV_NAME
    RUN_CMD="CABSflex"
else
    # Fallback to conda
    CONDA_BASE=$(conda info --base 2>/dev/null)
    if [ -n "$CONDA_BASE" ] && [ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
        source "${CONDA_BASE}/etc/profile.d/conda.sh"
        conda activate $ENV_NAME
        RUN_CMD="CABSflex"
    else
        # Fallback to direct path
        RUN_CMD="CABSflex"
    fi
fi

echo -e "${BLUE}🚀 Starting CABSflex Legacy and Modern Restraint Modes Verification Test Suite (using CLI binary)${NC}"

# Test Directories
TEST_ROOT="tests/test_sim"
mkdir -p "$TEST_ROOT"

# Function to run a test mode
run_mode_test() {
    local mode=$1
    local expected_warning=$2
    local output_dir="${TEST_ROOT}/test_${mode}"
    
    echo -e "\n${YELLOW}🧪 Testing Restraint Mode/Alias: ${mode}${NC}"
    
    # Remove output dir if exists
    rm -rf "$output_dir"
    
    # Run command and capture output
    local cmd="${RUN_CMD} -i tests/inputs/Helix_short.pdb -g ${mode} -a 1 -y 2 -s 2 -w ${output_dir} -v 2"
    echo -e "${BLUE}Running: ${cmd}${NC}"
    
    local log_file="${TEST_ROOT}/${mode}_run.log"
    eval "$cmd" &> "$log_file"
    local status=$?
    
    # Verify execution status
    if [ $status -ne 0 ] || [ ! -d "$output_dir" ]; then
        echo -e "${RED}❌ Mode '${mode}' FAILED to run (Exit code: ${status})${NC}"
        cat "$log_file"
        return 1
    fi

    # Verify if alias warnings are triggered correctly
    if [ -n "$expected_warning" ]; then
        # Use case-insensitive grep with robust regex to support line wrapping
        if grep -iq "$expected_warning" "$log_file"; then
            echo -e "${GREEN}✅ Warning trigger verified: '${expected_warning}' was logged successfully.${NC}"
        else
            echo -e "${RED}❌ Warning NOT triggered. Expected message containing: '${expected_warning}'${NC}"
            echo -e "${YELLOW}--- Log content ---${NC}"
            cat "$log_file"
            echo -e "${YELLOW}--- End of Log ---${NC}"
            return 1
        fi
    fi
    
    echo -e "${GREEN}✅ Mode/Alias '${mode}' PASSED${NC}"
    return 0
}

FAILED_TESTS=0

# 1. Test Modern Mode 'flexible'
run_mode_test "flexible" "" || ((FAILED_TESTS++))

# 2. Test Modern Mode 'rigid'
run_mode_test "rigid" "" || ((FAILED_TESTS++))

# 3. Test Native Resurrected Mode 'ss1'
run_mode_test "ss1" "" || ((FAILED_TESTS++))

# 4. Test Legacy Alias 'ss2' (should map to 'flexible')
run_mode_test "ss2" "ss2.*legacy" || ((FAILED_TESTS++))

# 5. Test Legacy Alias 'all' (should map to 'rigid')
run_mode_test "all" "all.*legacy" || ((FAILED_TESTS++))

# Summary
echo -e "\n${BLUE}========================================${NC}"
if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "${GREEN}⭐ ALL MODES AND ALIASES VERIFIED SUCCESSFULLY VIA CLINARY COMMAND${NC}"
    # Cleanup test outputs
    rm -rf "$TEST_ROOT"
    exit 0
else
    echo -e "${RED}💥 ${FAILED_TESTS} MODE(S)/ALIAS(ES) FAILED VERIFICATION${NC}"
    exit 1
fi
