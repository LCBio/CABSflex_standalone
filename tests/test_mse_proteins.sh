#!/bin/bash
# Regression test for bug: CABSflex failing on proteins with MSE (selenomethionine)
# residues with:
#   [CRITICAL] No sequential similarity between input and reference
#              according to used alignment method (trivial).
#
# Affected PDB IDs: 2gzr, 2nvp, 2phz, 3r5s  (all contain MSE HETATM records)
# Fix: MSE added to extended_amino_acids.json (commit 61f7660)

set +e
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# --- Environment activation (mirrors test_cli_options.sh) ---
if [ -f "install.sh" ]; then
    ENV_NAME=$(grep '^ENV_NAME=' install.sh | cut -d'=' -f2 | tr -d '"' | tr -d "'")
else
    ENV_NAME="cabs"
fi

export MAMBA_EXE="$HOME/.local/bin/micromamba"
if [ -f "$MAMBA_EXE" ]; then
    eval "$($MAMBA_EXE shell hook -s bash)"
    micromamba activate $ENV_NAME
else
    CONDA_BASE=$(conda info --base)
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate $ENV_NAME
fi

CABS_CMD="CABSflex"

SCENARIO_ROOT="tests/test_mse_proteins"
FAILED=0
PASSED=0

mkdir -p "$SCENARIO_ROOT"

run_case() {
    local pdb_id=$1
    local work_dir="${SCENARIO_ROOT}/${pdb_id}"

    echo -e "\n${YELLOW}Testing ${pdb_id}...${NC}"
    rm -rf "$work_dir"

    $CABS_CMD -i "${pdb_id}:A" \
        -w "$work_dir" \
        --aa-method cg2all \
        --aa-rebuild M \
        --save-config \
        --json-output \
        --dssp-output \
        --ss-output \
        --restraints-output \
        --csv-output A \
        --pdb-bfac-output A \
        --generate-pymol-visualizations \
        --generate-chimera-visualizations \
        --generate-notebook \
        --contact-maps \
        --renumber-residues-to-original \
        -s 10
        2>&1 | tee "${work_dir}.log"

    local exit_code=${PIPESTATUS[0]}

    if [ $exit_code -ne 0 ]; then
        echo -e "${RED}FAILED${NC} ${pdb_id}: non-zero exit (${exit_code})"
        ((FAILED++))
        return
    fi

    if grep -q "No sequential similarity" "${work_dir}.log"; then
        echo -e "${RED}FAILED${NC} ${pdb_id}: trivial-alignment error still present"
        ((FAILED++))
        return
    fi

    echo -e "${GREEN}PASSED${NC} ${pdb_id}"
    ((PASSED++))
}

echo -e "${BLUE}MSE protein regression test${NC}"
echo "PDB IDs: 2gzr, 2nvp, 2phz, 3r5s"
echo "These structures all contain selenomethionine (MSE) HETATM records."
echo

for pdb_id in 2gzr 2nvp 2phz 3r5s; do
    run_case "$pdb_id"
done

echo
echo "Results: ${PASSED} passed, ${FAILED} failed"

if [ $FAILED -gt 0 ]; then
    echo -e "${RED}REGRESSION DETECTED${NC}"
    exit 1
fi

echo -e "${GREEN}All MSE protein tests passed.${NC}"
