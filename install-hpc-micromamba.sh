#!/bin/bash -l
# ==============================================================================
# CABS-flex HPC Micromamba Installer
# ==============================================================================
set -euo pipefail

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
# Override these with environment variables or edit this block for site defaults.
BASE_INSTALL_DIR="${BASE_INSTALL_DIR:-}"
ENV_NAME="${ENV_NAME:-cabs}"
CG2ALL_ENV_NAME="${CG2ALL_ENV_NAME:-${ENV_NAME}_cg2all}"

INSTALL_MODELLER="${INSTALL_MODELLER:-TRUE}"
MODELLER_KEY="${MODELLER_KEY:-}"
MODELLER_VERSION="${MODELLER_VERSION:-10.7}"
MODELLER_ARCH_INDEX="${MODELLER_ARCH_INDEX:-2}"

TEMP_ROOT="${TEMP_ROOT:-${SCRATCH:-/tmp}}"
MICROMAMBA_BIN_DIR="${MICROMAMBA_BIN_DIR:-${BASE_INSTALL_DIR}/bin}"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-${BASE_INSTALL_DIR}/micromamba}"

USE_MODULES="${USE_MODULES:-FALSE}"
GCC_MODULE="${GCC_MODULE:-}"
BZIP2_MODULE="${BZIP2_MODULE:-}"
EXTRA_MODULES="${EXTRA_MODULES:-}"

PYTHON_VERSION_MAIN="${PYTHON_VERSION_MAIN:-3.10}"
PYTHON_VERSION_CG2ALL="${PYTHON_VERSION_CG2ALL:-3.9}"
TORCH_VERSION="${TORCH_VERSION:-2.2.0+cpu}"
TORCHVISION_VERSION="${TORCHVISION_VERSION:-0.17.0+cpu}"
CG2ALL_COMMIT="${CG2ALL_COMMIT:-a789cb5}"

RECREATE_ENVS="${RECREATE_ENVS:-TRUE}"

# ------------------------------------------------------------------------------
# Derived paths
# ------------------------------------------------------------------------------
PROJECT_ROOT="$(pwd)"
MAIN_ENV_DIR="${BASE_INSTALL_DIR}/envs/${ENV_NAME}"
CG2ALL_ENV_DIR="${BASE_INSTALL_DIR}/envs/${CG2ALL_ENV_NAME}"
PIP_CACHE_DIR="${TEMP_ROOT}/cabsflex-pip-cache"
MAMBA_PKGS_DIR="${TEMP_ROOT}/cabsflex-mamba-pkgs"

# Output colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}      CABS-flex Standalone HPC Micromamba Installer             ${NC}"
echo -e "${BLUE}================================================================${NC}"

die() {
    echo -e "${RED}❌ $*${NC}" >&2
    exit 1
}

info() {
    echo -e "${YELLOW}📦 $*${NC}"
}

success() {
    echo -e "${GREEN}✅ $*${NC}"
}

require_tool() {
    local tool="$1"
    command -v "$tool" >/dev/null 2>&1 || die "Required tool not found: $tool"
}

sedi() {
    sed -i "$@"
}

load_modules() {
    [ "${USE_MODULES}" = "TRUE" ] || return 0

    if ! command -v module >/dev/null 2>&1; then
        return 0
    fi

    module purge
    if [ -n "${GCC_MODULE}" ]; then
        module load "${GCC_MODULE}"
    fi
    if [ -n "${BZIP2_MODULE}" ]; then
        module load "${BZIP2_MODULE}"
    fi
    if [ -n "${EXTRA_MODULES}" ]; then
        # shellcheck disable=SC2206
        local extra_modules=( ${EXTRA_MODULES} )
        if [ "${#extra_modules[@]}" -gt 0 ]; then
            module load "${extra_modules[@]}"
        fi
    fi
}

run_with_retry() {
    local attempt=1
    local max_attempts=5
    local delay=5

    while true; do
        if "$@"; then
            return 0
        fi
        if [ "${attempt}" -ge "${max_attempts}" ]; then
            die "Command failed after ${max_attempts} attempts: $*"
        fi
        attempt=$((attempt + 1))
        echo -e "${YELLOW}⚠️  Command failed. Retrying in ${delay}s (${attempt}/${max_attempts})...${NC}"
        sleep "${delay}"
    done
}

bootstrap_micromamba() {
    mkdir -p "${MICROMAMBA_BIN_DIR}" "${MAMBA_ROOT_PREFIX}" "${PIP_CACHE_DIR}" "${MAMBA_PKGS_DIR}"
    export PATH="${MICROMAMBA_BIN_DIR}:$PATH"
    export MAMBA_ROOT_PREFIX
    export MAMBA_PKGS_DIRS="${MAMBA_PKGS_DIR}"

    if command -v micromamba >/dev/null 2>&1; then
        MAMBA_EXE="$(command -v micromamba)"
        success "Using existing micromamba at ${MAMBA_EXE}"
        export MAMBA_EXE
        return 0
    fi

    local os_type arch_type platform url
    os_type="$(uname -s | tr '[:upper:]' '[:lower:]')"
    arch_type="$(uname -m)"

    case "${os_type}-${arch_type}" in
        linux-x86_64) platform="linux-64" ;;
        linux-aarch64) platform="linux-aarch64" ;;
        *)
            die "Unsupported platform for HPC micromamba bootstrap: ${os_type}-${arch_type}"
            ;;
    esac

    url="https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-${platform}"
    info "Downloading micromamba standalone binary"
    run_with_retry curl -Ls "${url}" -o "${MICROMAMBA_BIN_DIR}/micromamba"
    chmod +x "${MICROMAMBA_BIN_DIR}/micromamba"
    MAMBA_EXE="${MICROMAMBA_BIN_DIR}/micromamba"
    export MAMBA_EXE
    success "Micromamba bootstrapped to ${MAMBA_EXE}"
}

reset_env_dir() {
    local env_dir="$1"
    if [ -d "${env_dir}" ] && [ "${RECREATE_ENVS}" = "TRUE" ]; then
        echo -e "${YELLOW}⚠️  Removing existing environment: ${env_dir}${NC}"
        rm -rf "${env_dir}"
    fi
}

write_runtime_hooks() {
    local env_dir="$1"
    mkdir -p "${env_dir}/etc/conda/activate.d" "${env_dir}/etc/conda/deactivate.d"

    cat > "${env_dir}/etc/conda/activate.d/cabsflex_hpc.sh" <<EOF
#!/bin/bash
export PYTHONNOUSERSITE=1
export DGL_DISABLE_GRAPHBOLT=1
EOF
    if [ -f "${env_dir}/etc/cabsflex_modeller_env.sh" ]; then
        cat >> "${env_dir}/etc/conda/activate.d/cabsflex_hpc.sh" <<EOF
source "${env_dir}/etc/cabsflex_modeller_env.sh"
EOF
    fi
    if [ "${USE_MODULES}" = "TRUE" ]; then
        cat >> "${env_dir}/etc/conda/activate.d/cabsflex_hpc.sh" <<EOF
if command -v module >/dev/null 2>&1; then
    module purge
EOF
        if [ -n "${GCC_MODULE}" ]; then
            echo "    module load \"${GCC_MODULE}\"" >> "${env_dir}/etc/conda/activate.d/cabsflex_hpc.sh"
        fi
        if [ -n "${BZIP2_MODULE}" ]; then
            echo "    module load \"${BZIP2_MODULE}\"" >> "${env_dir}/etc/conda/activate.d/cabsflex_hpc.sh"
        fi
        if [ -n "${EXTRA_MODULES}" ]; then
            echo "    module load ${EXTRA_MODULES}" >> "${env_dir}/etc/conda/activate.d/cabsflex_hpc.sh"
        fi
        cat >> "${env_dir}/etc/conda/activate.d/cabsflex_hpc.sh" <<'EOF'
fi
EOF
    fi
    chmod +x "${env_dir}/etc/conda/activate.d/cabsflex_hpc.sh"

    cat > "${env_dir}/etc/conda/deactivate.d/cabsflex_hpc.sh" <<'EOF'
#!/bin/bash
unset DGL_DISABLE_GRAPHBOLT
EOF
    chmod +x "${env_dir}/etc/conda/deactivate.d/cabsflex_hpc.sh"
}

write_activation_helper() {
    local helper_path="$1"
    local env_dir="$2"

    cat > "${helper_path}" <<EOF
#!/bin/bash -l
export PATH="${MICROMAMBA_BIN_DIR}:\$PATH"
export MAMBA_EXE="${MAMBA_EXE}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX}"
export MAMBA_PKGS_DIRS="${MAMBA_PKGS_DIR}"
EOF
    if [ "${USE_MODULES}" = "TRUE" ]; then
        cat >> "${helper_path}" <<EOF
if command -v module >/dev/null 2>&1; then
    module purge
EOF
        if [ -n "${GCC_MODULE}" ]; then
            echo "    module load \"${GCC_MODULE}\"" >> "${helper_path}"
        fi
        if [ -n "${BZIP2_MODULE}" ]; then
            echo "    module load \"${BZIP2_MODULE}\"" >> "${helper_path}"
        fi
        if [ -n "${EXTRA_MODULES}" ]; then
            echo "    module load ${EXTRA_MODULES}" >> "${helper_path}"
        fi
        cat >> "${helper_path}" <<EOF
fi
EOF
    fi
    cat >> "${helper_path}" <<EOF
eval "\$("\${MAMBA_EXE}" shell hook -s bash)"
micromamba activate "${env_dir}"
export PYTHONNOUSERSITE=1
export DGL_DISABLE_GRAPHBOLT=1
EOF
    if [ -f "${env_dir}/etc/cabsflex_modeller_env.sh" ]; then
        cat >> "${helper_path}" <<EOF
source "${env_dir}/etc/cabsflex_modeller_env.sh"
EOF
    fi
    chmod +x "${helper_path}"
}

configure_cabs_paths() {
    mkdir -p "${PROJECT_ROOT}/CABS/data"
    printf '{"cg2all_env_prefix": "%s", "cabs_env_prefix": "%s"}\n' "${CG2ALL_ENV_DIR}" "${MAIN_ENV_DIR}" > "${PROJECT_ROOT}/CABS/data/cabs_paths.json"
    success "Configured cg2all environment path in CABS/data/cabs_paths.json"
}

install_modeller_source() {
    [ "${INSTALL_MODELLER}" = "TRUE" ] || return 0
    [ -n "${MODELLER_KEY}" ] || {
        echo -e "${YELLOW}ℹ️  MODELLER_KEY is empty. Skipping source fallback.${NC}"
        return 0
    }

    local install_dir="${MAIN_ENV_DIR}/modeller"
    local arch_name="x86_64-intel8"

    info "Installing Modeller ${MODELLER_VERSION} from source fallback"
    run_with_retry curl -L "https://salilab.org/modeller/${MODELLER_VERSION}/modeller-${MODELLER_VERSION}.tar.gz" -o "modeller.tar.gz" --fail
    tar -xzf "modeller.tar.gz"
    cd "modeller-${MODELLER_VERSION}"

    ./Install <<EOF
${MODELLER_ARCH_INDEX}
${install_dir}
${MODELLER_KEY}


EOF

    local site_pkgs
    site_pkgs="$("${MAIN_ENV_DIR}/bin/python" -c 'import site; print(site.getsitepackages()[0])')"

    printf '%s\n' "${install_dir}/modlib" > "${site_pkgs}/modeller.pth"
    printf '%s\n' "${install_dir}/lib/${arch_name}/python3.3" >> "${site_pkgs}/modeller.pth"
    printf 'export LD_LIBRARY_PATH="%s:$LD_LIBRARY_PATH"\n' "${install_dir}/lib/${arch_name}" \
        > "${MAIN_ENV_DIR}/etc/cabsflex_modeller_env.sh"
    success "Modeller source fallback installed"
}

configure_modeller_license() {
    [ -n "${MODELLER_KEY}" ] || return 0

    local mod_config=""

    info "Configuring Modeller license"
    mod_config="$("${MAMBA_EXE}" run -p "${MAIN_ENV_DIR}" python -c "import modlib.modeller.config as c; print(c.__file__)" 2>/dev/null || echo "")"

    if [ -z "${mod_config}" ]; then
        mod_config="$("${MAMBA_EXE}" run -p "${MAIN_ENV_DIR}" python -c "import modeller; import os; print(os.path.join(os.path.dirname(modeller.__file__), 'config.py'))" 2>/dev/null || echo "")"
    fi

    if [ -z "${mod_config}" ] || [ ! -f "${mod_config}" ]; then
        mod_config="$(find "${MAIN_ENV_DIR}" -name "config.py" | grep "/modeller/" | head -n 1 || echo "")"
    fi

    if [ -n "${mod_config}" ] && [ -f "${mod_config}" ]; then
        sedi "s/license = .*/license = r'${MODELLER_KEY}'/" "${mod_config}"
        success "Modeller license configured in ${mod_config}"
    else
        echo -e "${YELLOW}⚠️  Could not find Modeller config file in ${MAIN_ENV_DIR}.${NC}"
    fi
}

verify_import() {
    local env_dir="$1"
    local module_name="$2"
    "${env_dir}/bin/python" -c "import ${module_name}" >/dev/null 2>&1
}

verify_cli() {
    local bin_path="$1"
    local bin_name
    bin_name="$(basename "${bin_path}")"
    if [ -x "${bin_path}" ] && "${bin_path}" --help >/dev/null 2>&1; then
        success "${bin_name} confirmed and functional"
    else
        echo -e "${RED}❌ ${bin_name} not found or failed to execute${NC}"
        return 1
    fi
}

cleanup_temp_dir() {
    if [ -n "${INSTALL_TEMP_DIR:-}" ] && [ -d "${INSTALL_TEMP_DIR}" ]; then
        rm -rf "${INSTALL_TEMP_DIR}"
    fi
}

main() {
    [ -n "${BASE_INSTALL_DIR}" ] || die "BASE_INSTALL_DIR must be set to a persistent installation path (do not use scratch for the install root)"
    [ -f "${PROJECT_ROOT}/requirements-runtime.txt" ] || die "requirements-runtime.txt not found in ${PROJECT_ROOT}"
    [ -f "${PROJECT_ROOT}/pyproject.toml" ] || [ -f "${PROJECT_ROOT}/setup.py" ] || \
        die "Run this script from the root of the CABS-flex repository"

    require_tool git
    require_tool curl
    require_tool tar
    require_tool sed

    mkdir -p "${BASE_INSTALL_DIR}" "${BASE_INSTALL_DIR}/envs"
    load_modules
    bootstrap_micromamba

    # Isolate the build/install from any system or module-loaded Python state.
    # LD_LIBRARY_PATH is left alone since USE_MODULES / GCC_MODULE may depend on it.
    export PYTHONPATH=""
    export PYTHONHOME=""
    export PYTHONUSERBASE=""

    INSTALL_TEMP_DIR="$(mktemp -d -p "${TEMP_ROOT}" "cabs_mamba_install_XXXX")"
    trap cleanup_temp_dir EXIT
    cd "${INSTALL_TEMP_DIR}"

    reset_env_dir "${MAIN_ENV_DIR}"
    reset_env_dir "${CG2ALL_ENV_DIR}"

    configure_cabs_paths

    # gfortran/binutils provide a Fortran toolchain for deps that build from source
    # (mirrors install.sh). Every requirements-runtime.txt package with a C/C++/
    # Fortran extension (numpy, scipy, matplotlib+pillow, biopython, mdtraj, hdf5/
    # netcdf4/h5py) is installed here as a conda-forge binary rather than left for
    # pip to build, since old HPC clusters commonly have a glibc/system-library
    # baseline (e.g. HDF5 1.8.x, libpng 1.5) too old for manylinux wheels or for a
    # from-source build to succeed against.
    local main_specs=( "python=${PYTHON_VERSION_MAIN}" "pip" "openmm" "gfortran" "binutils" "numpy" "scipy" "matplotlib" "pillow" "biopython" "mdtraj" "hdf5" "libnetcdf" "netcdf4" "h5py" "pandas" "plotly" "jupyter" "nbconvert" )
    if [ "${INSTALL_MODELLER}" = "TRUE" ]; then
        main_specs+=( "modeller" )
    fi

    info "Creating main CABS-flex environment"
    "${MAMBA_EXE}" create -y -p "${MAIN_ENV_DIR}" \
        -c conda-forge -c bioconda -c salilab --override-channels \
        "${main_specs[@]}"
    write_runtime_hooks "${MAIN_ENV_DIR}"
    configure_modeller_license

    info "Installing runtime dependencies into main environment"
    run_with_retry "${MAMBA_EXE}" run -p "${MAIN_ENV_DIR}" pip install \
        --cache-dir "${PIP_CACHE_DIR}" -r "${PROJECT_ROOT}/requirements-runtime.txt"

    info "Cleaning stale build artifacts from ${PROJECT_ROOT}"
    rm -rf "${PROJECT_ROOT}/tests/test_cli_options" "${PROJECT_ROOT}/build" "${PROJECT_ROOT}/dist" "${PROJECT_ROOT}"/*.egg-info

    info "Installing local CABS-flex package"
    run_with_retry "${MAMBA_EXE}" run -p "${MAIN_ENV_DIR}" pip install \
        --cache-dir "${PIP_CACHE_DIR}" "${PROJECT_ROOT}"

    if [ "${INSTALL_MODELLER}" = "TRUE" ] && ! verify_import "${MAIN_ENV_DIR}" modeller; then
        echo -e "${YELLOW}⚠️  Modeller import failed after micromamba install. Trying source fallback.${NC}"
        install_modeller_source
        configure_modeller_license
        write_runtime_hooks "${MAIN_ENV_DIR}"
    fi

    info "Creating isolated cg2all reconstruction environment"
    "${MAMBA_EXE}" create -y -p "${CG2ALL_ENV_DIR}" \
        -c conda-forge --override-channels \
        "python=${PYTHON_VERSION_CG2ALL}" "pip" "c-compiler" "cxx-compiler" "make"
    write_runtime_hooks "${CG2ALL_ENV_DIR}"

    info "Installing PyTorch and TorchVision for cg2all"
    run_with_retry "${MAMBA_EXE}" run -p "${CG2ALL_ENV_DIR}" pip install \
        --cache-dir "${PIP_CACHE_DIR}" \
        "torch==${TORCH_VERSION}" \
        "torchvision==${TORCHVISION_VERSION}" \
        --index-url "https://download.pytorch.org/whl/cpu"

    info "Installing cg2all Python dependencies"
    run_with_retry "${MAMBA_EXE}" run -p "${CG2ALL_ENV_DIR}" pip install \
        --cache-dir "${PIP_CACHE_DIR}" "psutil>=5.8.0" tqdm
    run_with_retry "${MAMBA_EXE}" run -p "${CG2ALL_ENV_DIR}" pip install \
        --cache-dir "${PIP_CACHE_DIR}" --no-deps "dgl==1.1.3" \
        -f "https://data.dgl.ai/wheels/repo.html"
    run_with_retry "${MAMBA_EXE}" run -p "${CG2ALL_ENV_DIR}" pip install \
        --cache-dir "${PIP_CACHE_DIR}" --no-binary e3nn e3nn
    run_with_retry "${MAMBA_EXE}" run -p "${CG2ALL_ENV_DIR}" pip install \
        --cache-dir "${PIP_CACHE_DIR}" git+https://github.com/huhlim/mdtraj

    info "Installing patched SE3Transformer"
    local se3t_src="${INSTALL_TEMP_DIR}/se3t-src"
    run_with_retry git clone https://github.com/huhlim/SE3Transformer "${se3t_src}"
    cd "${se3t_src}"
    sed -i 's/python = "[^"]*"/python = ">=3.7"/' pyproject.toml
    sed -i 's/torch = "[^"]*"/torch = ">=2.1.0"/' pyproject.toml
    run_with_retry "${MAMBA_EXE}" run -p "${CG2ALL_ENV_DIR}" pip install \
        --cache-dir "${PIP_CACHE_DIR}" .

    info "Installing patched cg2all"
    local cg2all_src="${INSTALL_TEMP_DIR}/cg2all-src"
    run_with_retry git clone https://github.com/huhlim/cg2all.git "${cg2all_src}"
    cd "${cg2all_src}"
    git checkout "${CG2ALL_COMMIT}"
    sed -i 's/torch = "[^"]*"/torch = ">=2.1.0"/' pyproject.toml
    sed -i 's/numpy = "[^"]1"/numpy = ">=1.21"/' pyproject.toml
    run_with_retry "${MAMBA_EXE}" run -p "${CG2ALL_ENV_DIR}" pip install \
        --cache-dir "${PIP_CACHE_DIR}" .

    write_activation_helper "${BASE_INSTALL_DIR}/activate-cabs.sh" "${MAIN_ENV_DIR}"
    write_activation_helper "${BASE_INSTALL_DIR}/activate-cg2all.sh" "${CG2ALL_ENV_DIR}"

    info "Verifying main environment"
    "${MAMBA_EXE}" run -p "${MAIN_ENV_DIR}" python - <<'EOF'
import importlib

modules = ["Bio.PDB", "mdtraj"]
for name in modules:
    importlib.import_module(name)
print("Main environment imports verified.")
EOF

    if [ "${INSTALL_MODELLER}" = "TRUE" ]; then
        if verify_import "${MAIN_ENV_DIR}" modeller; then
            success "Modeller import verified"
        else
            echo -e "${YELLOW}⚠️  Modeller import still unavailable. Installation completed without verified Modeller support.${NC}"
        fi
    fi

    info "Verifying reconstruction environment"
    "${MAMBA_EXE}" run -p "${CG2ALL_ENV_DIR}" python - <<'EOF'
import cg2all
print("cg2all import verified.")
EOF

    [ -x "${CG2ALL_ENV_DIR}/bin/convert_cg2all" ] || die "convert_cg2all binary not found in ${CG2ALL_ENV_DIR}/bin"

    info "Verifying CLI entry points"
    verify_cli "${MAIN_ENV_DIR}/bin/CABSflex" || true
    verify_cli "${MAIN_ENV_DIR}/bin/CABSdock" || true

    success "Installation complete"
    echo "============================================================"
    echo "Main environment: ${MAIN_ENV_DIR}"
    echo "cg2all environment: ${CG2ALL_ENV_DIR}"
    echo "Micromamba root: ${MAMBA_ROOT_PREFIX}"
    echo
    echo "To start:"
    echo "  source ${BASE_INSTALL_DIR}/activate-cabs.sh"
    echo "  CABSflex --help"
    echo "  CABSdock --help"
    echo "============================================================"
}

main "$@"
