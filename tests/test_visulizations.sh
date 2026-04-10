#!/bin/bash
# ==============================================================================
#           Visualization CLI Integration Tests for CABSflex and CABSdock
# ==============================================================================

set +e
GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

set -m

if [ -f "install.sh" ]; then
    ENV_NAME=$(grep '^ENV_NAME=' install.sh | cut -d'=' -f2 | tr -d '"' | tr -d "'")
else
    ENV_NAME="cabs"
fi

CABSFLEX_CMD="micromamba run -n $ENV_NAME CABSflex"
CABSDOCK_CMD="micromamba run -n $ENV_NAME CABSdock"

SCENARIO_ROOT="tests/test_visulizations"
RESULTS_FILE="tests/test_visulizations_results.log"

rm -rf "$SCENARIO_ROOT" "$RESULTS_FILE"
mkdir -p "$SCENARIO_ROOT"

if [[ "$OSTYPE" == "darwin"* ]]; then
    NPROC=$(sysctl -n hw.ncpu)
else
    NPROC=$(nproc)
fi
MAX_JOBS=$((NPROC - 1))
[ $MAX_JOBS -lt 1 ] && MAX_JOBS=1

trap 'echo "Received signal, killing children"; jobs -p | xargs -r kill; exit 1' SIGINT SIGTERM

limit_jobs() {
    while [ $(jobs -r | wc -l) -ge $MAX_JOBS ]; do
        sleep 1
    done
    true
}

# Full simulations only: do not pass -a, -y, -s, -r, -k, or other speed-up options.
FLEX_PDB="tests/inputs/2BZ6.pdb"
DOCK_PDB="tests/inputs/1A2K.pdb"
FLEX_CASE="$CABSFLEX_CMD -i $FLEX_PDB:LH"
DOCK_CASE="$CABSDOCK_CMD -i $DOCK_PDB:AB -p KYVATLGV:EECTTTTC -p TAGQEKFGGLRDGYYI:CCCHHHTCCCCHHHHC"

escape_md() {
    echo "$1" | sed 's/|/\\|/g' | tr '\n' ' '
}

check_expected_files() {
    local out_dir="$1"
    shift
    local missing=""

    for path in "$@"; do
        if [ ! -e "${out_dir}/${path}" ]; then
            missing="${missing};missing ${path}"
        fi
    done

    echo "$missing"
}

run_visualization_test() {
    local name="$1"
    local cmd="$2"
    shift 2
    local out_dir="${SCENARIO_ROOT}/${name}"

    mkdir -p "$out_dir"

    echo -e "${BLUE}Running ${name}${NC}"
    eval "$cmd -w $out_dir > ${out_dir}/stdout.log 2> ${out_dir}/stderr.log"
    local status=$?

    local warnings=""
    local errors=""
    local critical=""
    local missing=""

    if [ -f "${out_dir}/CABS.log" ]; then
        warnings=$(grep "WARNING" "${out_dir}/CABS.log" | head -n 2 | tr '\n' ';' | sed 's/;$//')
        errors=$(grep "ERROR" "${out_dir}/CABS.log" | head -n 2 | tr '\n' ';' | sed 's/;$//')
        critical=$(grep "CRITICAL" "${out_dir}/CABS.log" | head -n 2 | tr '\n' ';' | sed 's/;$//')
        [ -n "$critical" ] && errors="${errors};${critical}"
    fi

    if [ $status -ne 0 ]; then
        local stderr_err
        stderr_err=$(tail -n 1 "${out_dir}/stderr.log" | head -c 200 | tr '\n' ' ' | sed 's/|/ /g')
        [ -z "$stderr_err" ] && stderr_err="Unknown error (exit code $status)"
        errors="${errors};${stderr_err}"
    fi

    if [ $status -eq 0 ]; then
        missing=$(check_expected_files "$out_dir" "$@")
        if [ -n "$missing" ]; then
            status=1
            errors="${errors}${missing}"
        fi
    fi

    local result="PASSED"
    [ $status -ne 0 ] && result="FAILED"

    warnings=$(escape_md "$warnings")
    errors=$(escape_md "$errors")

    echo "$name|$warnings|$errors|$result" >> "$RESULTS_FILE"
}

echo -e "${BLUE}Starting visualization integration tests.${NC}"
echo -e "${YELLOW}Warning: these use full default simulations and may take a while.${NC}"
echo -e "${YELLOW}Using Max Jobs: $MAX_JOBS${NC}"

run_visualization_test \
    "CABSflex_generate_pymol_visualizations" \
    "$FLEX_CASE --generate-pymol-visualizations" \
    "load_models.pml" \
    "color_by_ss.pml" \
    "color_by_rmsf.pml" \
    "animate_models.pml" \
    "load_restraints.pml" \
    "output_pdbs/start.pdb" \
    "output_pdbs/model_0.pdb" \
    "output_pdbs/start_rmsf.pdb" \
    "output_data/restraints.txt" \
    "output_data/ss.txt" \
    "contact_maps/all.txt" & limit_jobs

run_visualization_test \
    "CABSflex_generate_chimera_visualizations" \
    "$FLEX_CASE --generate-chimera-visualizations" \
    "color_rmsf.cxc" \
    "color_chain.cxc" \
    "color_ss.cxc" \
    "rmsf_worm.cxc" \
    "record_movie.cxc" \
    "load_restraints.cxc" \
    "output_pdbs/start.pdb" \
    "output_pdbs/model_0.pdb" \
    "output_pdbs/start_rmsf.pdb" \
    "output_data/restraints.txt" \
    "output_data/ss.txt" \
    "contact_maps/all.txt" & limit_jobs

run_visualization_test \
    "CABSflex_generate_notebook" \
    "$FLEX_CASE --generate-notebook --cg2all-representation calpha-sc" \
    "report.ipynb" \
    "report.html" \
    "output_pdbs/start.pdb" \
    "output_pdbs/model_0.pdb" \
    "output_pdbs/start_rmsf.pdb" \
    "output_data/restraints.txt" \
    "output_data/ss.txt" \
    "contact_maps/all.txt" & limit_jobs

run_visualization_test \
    "CABSdock_generate_pymol_visualizations" \
    "$DOCK_CASE --generate-pymol-visualizations" \
    "load_models.pml" \
    "color_by_ss.pml" \
    "color_by_rmsf.pml" \
    "animate_models.pml" \
    "load_restraints.pml" \
    "output_pdbs/start.pdb" \
    "output_pdbs/model_0.pdb" \
    "output_pdbs/start_rmsf.pdb" \
    "output_data/restraints.txt" \
    "output_data/ss.txt" & limit_jobs

run_visualization_test \
    "CABSdock_generate_chimera_visualizations" \
    "$DOCK_CASE --generate-chimera-visualizations" \
    "color_rmsf.cxc" \
    "color_chain.cxc" \
    "color_ss.cxc" \
    "rmsf_worm.cxc" \
    "record_movie.cxc" \
    "load_restraints.cxc" \
    "output_pdbs/start.pdb" \
    "output_pdbs/model_0.pdb" \
    "output_pdbs/start_rmsf.pdb" \
    "output_data/restraints.txt" \
    "output_data/ss.txt" & limit_jobs

run_visualization_test \
    "CABSdock_generate_notebook" \
    "$DOCK_CASE --generate-notebook --cg2all-representation calpha-sc" \
    "report.ipynb" \
    "report.html" \
    "output_pdbs/start.pdb" \
    "output_pdbs/model_0.pdb" \
    "output_pdbs/start_rmsf.pdb" \
    "output_data/restraints.txt" \
    "output_data/ss.txt" & limit_jobs

echo -e "${YELLOW}Waiting for all visualization tests to complete.${NC}"
wait

echo -e "\n${BLUE}========================================${NC}"
echo -e "${BLUE}CABS Visualization Integration Report${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "| Test | Warnings | Errors | Passed/Failed |"
echo -e "|------|----------|--------|---------------|"

sort "$RESULTS_FILE" | while IFS='|' read -r test_name warn err res; do
    if [ "$res" == "FAILED" ]; then
        echo -e "| $test_name | $warn | ${RED}$err${NC} | ${RED}$res${NC} |"
    else
        echo -e "| $test_name | $warn | $err | ${GREEN}$res${NC} |"
    fi
done

echo -e "${BLUE}========================================${NC}"

if grep -q '|FAILED$' "$RESULTS_FILE"; then
    exit 1
fi

exit 0
