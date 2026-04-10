#!/bin/bash
# Regression script for peptide RMSD fallback modes using 1D4T.

set +e
set -m

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

if [ -f "install.sh" ]; then
    ENV_NAME=$(grep '^ENV_NAME=' install.sh | cut -d'=' -f2 | tr -d '"' | tr -d "'")
else
    ENV_NAME="cabs"
fi

if [ -f "$HOME/.local/bin/micromamba" ]; then
    export MAMBA_EXE="$HOME/.local/bin/micromamba"
    eval "$($MAMBA_EXE shell hook -s bash)"
    micromamba activate "$ENV_NAME"
else
    CONDA_BASE=$(conda info --base)
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "$ENV_NAME"
fi

CABS_CMD="CABSdock"
TEST_ROOT="tests/test_peptide_rmsd_modes_1d4t"
RESULTS_FILE="${TEST_ROOT}/results.log"
INPUT_PDB="tests/inputs/1D4T.pdb"
REFERENCE="${INPUT_PDB}:A:B"

COMMON_FLAGS="-i ${INPUT_PDB}:A -R ${REFERENCE} -a 20 -y 50 -A -M -C -S -o A --renumber-residues-to-original --pdb-bfac-output A --csv-output A --json-output --dssp-output --ss-output --restraints-output --image-file-format svg --log"

rm -rf "$TEST_ROOT"
mkdir -p "$TEST_ROOT"

if [[ "$OSTYPE" == "darwin"* ]]; then
    NPROC=$(sysctl -n hw.ncpu)
else
    NPROC=$(nproc)
fi
MAX_JOBS=$((NPROC - 1))
[ $MAX_JOBS -lt 1 ] && MAX_JOBS=1

trap 'jobs -p | xargs -r kill; exit 1' SIGINT SIGTERM

limit_jobs() {
    while [ "$(jobs -r | wc -l)" -ge "$MAX_JOBS" ]; do
        sleep 1
    done
}

run_case() {
    local name="$1"
    local peptide="$2"
    local mode="$3"
    local expected="$4"

    local work_dir="${TEST_ROOT}/${name}"
    local mode_flag=""
    if [ -n "${mode}" ]; then
        mode_flag="--peptide-rmsd-mode ${mode}"
    fi
    local cmd="${CABS_CMD} ${COMMON_FLAGS} -p ${peptide} ${mode_flag} -w ${work_dir}"

    echo -e "${BLUE}Running ${name}${NC}"
    echo "${cmd}" > "${work_dir}.command"
    eval "${cmd} > ${work_dir}.stdout 2> ${work_dir}.stderr"
    local status=$?

    local summary="${work_dir}/output_data/peptide_alignment_summary.csv"
    local selected=""
    local detail="no-summary"

    if [ -f "${summary}" ]; then
        selected=$(awk -F',' 'NR==2 {print $3}' "${summary}" | tr -d '\r')
        detail=$(tail -n +2 "${summary}" | head -n 1 | tr -d '\r')
    fi

    local result="FAILED"
    if [ $status -eq 0 ] && [ "${selected}" = "${expected}" ]; then
        result="PASSED"
    fi

    printf "%-24s | %-10s | %-10s | %-6s | %s\n" \
        "${name}" "${expected}" "${selected:-missing}" "${result}" "${detail}" | tee -a "${RESULTS_FILE}"

    if [ "${result}" = "FAILED" ]; then
        echo -e "${RED}${name} failed${NC}"
        if [ -f "${work_dir}/CABS.log" ]; then
            tail -n 20 "${work_dir}/CABS.log"
        fi
        return 1
    fi
    return 0
}

echo "case                     | expected   | selected   | result | summary" > "${RESULTS_FILE}"
echo "-------------------------|------------|------------|--------|--------" >> "${RESULTS_FILE}"

# Exact native peptide from 1D4T chain B. Should use strict mode.
run_case "strict_native" "KSLTIYAQVQK" "strict" "strict" &
limit_jobs

# Default mode should also resolve to strict for the exact native peptide.
run_case "default_strict_native" "KSLTIYAQVQK" "" "strict" &
limit_jobs

# N-terminal overhang beyond the resolved crystal fragment. Should use overlap mode.
run_case "overlap_overhang" "QKSLTIYAQVQK" "overlap" "overlap" &
limit_jobs

# Default mode should fall back from strict to overlap for the overhang case.
run_case "default_overlap_overhang" "QKSLTIYAQVQK" "" "overlap" &
limit_jobs

# Screening-like variant with several substitutions. Should require mutational mode.
run_case "mutational_variant" "KALTVYVQIQR" "mutational" "mutational" &
limit_jobs

# Default mode should fall through to mutational for the screening-like variant.
run_case "default_mutational_variant" "KALTVYVQIQR" "" "mutational" &
limit_jobs

# Negative control: too dissimilar to pass the current thresholds.
run_case "rejected_control" "RALTVFVIVQR" "mutational" "skipped" &

# Default mode should reject the negative control after exhausting all fallbacks.
run_case "default_rejected_control" "RALTVFVIVQR" "" "skipped" &

FAILURES=0
for job in $(jobs -p); do
    wait "$job" || FAILURES=$((FAILURES + 1))
done

if [ $FAILURES -eq 0 ]; then
    echo -e "${GREEN}All 1D4T peptide RMSD mode tests passed.${NC}"
    exit 0
fi

echo -e "${RED}${FAILURES} test case(s) failed.${NC}"
exit 1
