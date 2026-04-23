#!/bin/bash
# ==============================================================================
#      Full Parallel/Batch Reconstruction Suite for CABSflex/dock
# ==============================================================================

# --- Configuration ---
set +e
GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

RUN_CABSFLEX=1
RUN_CABSDOCK=1
SHORT_RUN=0
ANNEALING_CYCLES=50
TRAJ_SAVE_FREQ=10
CLUSTER_FRAMES=10
SAVE_ALL_MODELS="-S"

print_usage() {
    cat <<EOF
Usage: $0 [--cabsflex-only] [--cabsdock-only] [--short] [--help]

Options:
  --cabsflex-only   Run only CABSflex scenarios.
  --cabsdock-only   Run only CABSdock scenarios.
  --short           Run a shorter, faster version of the suite.
  --help            Show this help message.
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --cabsflex-only)
            RUN_CABSFLEX=1
            RUN_CABSDOCK=0
            ;;
        --cabsdock-only)
            RUN_CABSFLEX=0
            RUN_CABSDOCK=1
            ;;
        --short)
            SHORT_RUN=1
            ;;
        --help|-h)
            print_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            print_usage >&2
            exit 1
            ;;
    esac
    shift
done

if [ "$SHORT_RUN" -eq 1 ]; then
    ANNEALING_CYCLES=10
    TRAJ_SAVE_FREQ=5
    CLUSTER_FRAMES=5
fi

ENV_NAME="cabs"
SCENARIO_ROOT="tests/test_parallel_native"
RESULTS_FILE="tests/test_parallel_native_results.log"
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

COMMON_OPTS="--save-config --json-output --dssp-output --ss-output --restraints-output --csv-output A --pdb-bfac-output A --generate-pymol-visualizations --generate-chimera-visualizations --generate-notebook --contact-maps --renumber-residues-to-original"

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

check_all_atom() {
    local file="$1"
    if [ ! -f "$file" ] && [ ! -d "$file" ]; then
        echo "FAILED (Not found: $file)"
        return 1
    fi
    # If it's a DCD, we just check existence for now as grep won't work easily on binary
    if [[ "$file" == *.dcd ]]; then
        echo "PASSED (DCD exists)"
        return 0
    fi
    # Check for atoms other than CA and SC (N, O)
    local n_count=$(grep -c " N   " "$file")
    local o_count=$(grep -c " O   " "$file")
    if [ "$n_count" -gt 0 ] && [ "$o_count" -gt 0 ]; then
        echo "PASSED ($n_count N, $o_count O)"
        return 0
    else
        echo "FAILED (No AA details: $n_count N, $o_count O)"
        return 1
    fi
}

check_logs() {
    local out_dir="$1"
    local pattern="$2"
    local log_file=$(ls ${out_dir}/output_data/*.log 2>/dev/null | head -n 1)
    if [ -z "$log_file" ] || [ ! -f "$log_file" ]; then
        echo "FAILED (Log missing)"
        return 1
    fi
    if grep -q "$pattern" "$log_file"; then
        echo "PASSED"
        return 0
    else
        echo "FAILED (Pattern not found)"
        return 1
    fi
}

run_test() {
    local case_name="$1"
    local cmd="$2"
    local verify_func="$3"
    local out_dir="${SCENARIO_ROOT}/${case_name}"
    mkdir -p "$out_dir"
    echo -e "${BLUE}Running: $case_name${NC}"
    
    eval "$cmd -w $out_dir -a $ANNEALING_CYCLES -y $TRAJ_SAVE_FREQ -s $CLUSTER_FRAMES $SAVE_ALL_MODELS > ${out_dir}/stdout.log 2> ${out_dir}/stderr.log"
    local status=$?
    
    local result="PASSED"
    local details=""
    if [ $status -ne 0 ]; then
        result="FAILED"
        details="Exit code $status"
    else
        # Use a temporary function to pass out_dir correctly
        verify_cmd=$(echo "$verify_func" | sed "s|OUT_DIR|$out_dir|g")
        details=$(eval "$verify_cmd")
        if [[ "$details" == "FAILED"* ]]; then result="FAILED"; fi
    fi
    echo "$case_name|$result|$details" > "${out_dir}/result.tmp"
}

# --- Test Cases ---

echo -e "${BLUE}🚀 Starting Native Parallel Reconstruction Suite.${NC}"

if [ "$RUN_CABSFLEX" -eq 1 ]; then
    # 1. FLEX Baseline (Medoids Only)
    run_test "FLEX_MEDOIDS" "$RUN_CMD CABSflex -i $FLEX_PDB:LH --aa-method cg2all $COMMON_OPTS" \
    "check_all_atom OUT_DIR/output_pdbs/model_0.pdb && check_logs OUT_DIR 'Saving final models (in AA representation)'" & limit_jobs

    # 2. FLEX Clusters (Parallel)
    run_test "FLEX_CLUSTERS" "$RUN_CMD CABSflex -i $FLEX_PDB:LH --aa-method cg2all -A C $COMMON_OPTS" \
    "check_all_atom OUT_DIR/output_pdbs/cluster_0.pdb && check_logs OUT_DIR 'parallel reconstruction for 10 clusters'" & limit_jobs

    # 2.5. FLEX Trajectories (Batch DCD)
    run_test "FLEX_TRAJECTORY" "$RUN_CMD CABSflex -i $FLEX_PDB:LH --aa-method cg2all -A T $COMMON_OPTS" \
    "check_all_atom OUT_DIR/output_pdbs/replica_all_atom.dcd && check_logs OUT_DIR 'batch reconstruction for 1 trajectories'" & limit_jobs

    # 3. FLEX All Reconstruction
    run_test "FLEX_ALL" "$RUN_CMD CABSflex -i $FLEX_PDB:LH --aa-method cg2all -A A $COMMON_OPTS" \
    "check_all_atom OUT_DIR/output_pdbs/model_0.pdb && check_all_atom OUT_DIR/output_pdbs/cluster_0.pdb && check_all_atom OUT_DIR/output_pdbs/replica_all_atom.dcd" & limit_jobs
fi

if [ "$RUN_CABSDOCK" -eq 1 ]; then
    # 4. DOCK Unified
    run_test "DOCK_UNIFIED" "$RUN_CMD CABSdock -i $DOCK_REC -p KYVATLGV:EECTTTTC --aa-method cg2all -A A $COMMON_OPTS" \
    "check_all_atom OUT_DIR/output_pdbs/model_0.pdb && check_all_atom OUT_DIR/output_pdbs/cluster_0.pdb && check_all_atom OUT_DIR/output_pdbs/replica_all_atom.dcd" & limit_jobs
fi

wait

# Compile and Report
echo -e "\n${BLUE}========================================${NC}"
echo -e "${BLUE}   Parallel Reconstruction Report       ${NC}"
echo -e "${BLUE}========================================${NC}"
printf "| %-25s | %-10s | %-40s |\n" "Case" "Result" "Details"
echo -e "|---------------------------|------------|------------------------------------------|"
if ls ${SCENARIO_ROOT}/*/result.tmp >/dev/null 2>&1; then
    sort ${SCENARIO_ROOT}/*/result.tmp | while IFS='|' read -r case res det; do
        [ "$res" == "FAILED" ] && COL=$RED || COL=$GREEN
        printf "| %-25s | ${COL}%-10s${NC} | %-40s |\n" "$case" "$res" "$det"
    done
else
    echo "| No cases were run.                                              |"
fi
echo -e "${BLUE}========================================${NC}"
