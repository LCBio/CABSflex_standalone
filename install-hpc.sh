#!/bin/bash -l
# ==============================================================================
# 🧬 CABS-flex Local Installer (Beta - Helios Optimized)
# ==============================================================================
set -e

# ==============================================================================
# --- 1. CONFIGURATION SECTION ---
# ==============================================================================
INSTALL_MODELLER="TRUE"
MODELLER_KEY=""   # <<< SET YOUR MODELLER LICENSE KEY HERE
MODELLER_VERSION="10.7"
MODELLER_ARCH_INDEX="2"      # 2 = x86_64-intel8

BASE_INSTALL_DIR=""   #<<< SET THE PATH TO INSTALL
VENV_NAME="cabs"
VENV_DIR="$BASE_INSTALL_DIR/$VENV_NAME"
CG2ALL_VENV_DIR="$BASE_INSTALL_DIR/${VENV_NAME}_reconstruct"
TEMP_ROOT="$SCRATCH" # Use scratch for high-I/O operations

# HPC Environment Modules
GCC_MODULE="GCCcore/13.2.0"
PYTHON_MODULE="Python/3.11.5"
BZIP2_MODULE="bzip2/1.0.8"
INTEL_MODULE="HDF5/1.14.3-serial"
NETCDF_MODULE="impi/2021.10.0 netCDF/4.9.3 "
HDF5_MODULE="intel-compilers/2023.2.1 HDF5/1.14.3-serial"
# Dependencies (Kept lean: MDTraj replaces DSSP binary)
CORE_DEPS=("numpy" "matplotlib" "requests" "biopython" "mdtraj" "biopandas" "tqdm" "scipy" "pandas" "plotly" "jupyter" "nbconvert" "ipymolstar")
TORCH_URL="https://download.pytorch.org/whl/cpu"
NSP3_REPO_URL="https://github.com/Eryk96/NetSurfP-3.0.git"

# Unlike the settings above, the two options below are read from the
# environment at run time (not just this file), so they don't need editing
# here — e.g. `E3NN_USE_WHEEL=TRUE bash install-hpc.sh`.
#
# e3nn is built from source (--no-binary) by default for compatibility with
# older cluster kernels/glibc. If that build fails, re-run with
# E3NN_USE_WHEEL=TRUE to install the prebuilt wheel instead. See the wiki
# Installation page, Section 6 (Troubleshooting & FAQ).
E3NN_USE_WHEEL="${E3NN_USE_WHEEL:-FALSE}"

# Set TRUE to reuse the existing main and cg2all reconstruction venvs
# instead of rebuilding them from scratch (each is only reused if actually
# found present and functional). Used to resume installation after manually
# placing a cg2all checkpoint file (see Section 6 of the wiki Installation
# page) without repeating the whole build.
RESUME_AFTER_CHECKPOINT="${RESUME_AFTER_CHECKPOINT:-FALSE}"
# ==============================================================================

# Output Colors
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

# Cross-platform MD5 helper (Linux: md5sum, macOS: md5)
md5_of() {
    if command -v md5sum &> /dev/null; then
        md5sum "$1" | awk '{print $1}'
    else
        md5 -q "$1"
    fi
}

# Download a cg2all checkpoint from Zenodo and verify its MD5, retrying a few
# times. If it still can't be verified, print manual-download instructions
# and return failure instead of installing a corrupt/blocked-response file.
fetch_checkpoint() {
    local name="$1" expected_md5="$2" model_home="$3"
    local dest="$model_home/$name"
    local url="https://zenodo.org/record/8393343/files/$name"
    local attempt=1 max=3 delay=5

    mkdir -p "$model_home"

    if [ -f "$dest" ] && [ "$(md5_of "$dest")" = "$expected_md5" ]; then
        echo -e "${GREEN}✅ $name already present and MD5-verified.${NC}"
        return 0
    fi

    while true; do
        echo -e "${YELLOW}📥 Downloading $name from Zenodo (attempt $attempt/$max)...${NC}"
        rm -f "$dest"
        if curl -fL -o "$dest" "$url" && [ "$(md5_of "$dest")" = "$expected_md5" ]; then
            echo -e "${GREEN}✅ $name downloaded and MD5-verified.${NC}"
            return 0
        fi
        rm -f "$dest"
        if [[ $attempt -ge $max ]]; then
            echo -e "${RED}❌ Could not obtain a valid $name after $max attempts.${NC}"
            echo -e "${YELLOW}ℹ️  This can happen when Zenodo blocks automated downloads from some networks.${NC}"
            echo -e "${YELLOW}   Please download it manually and place it at:${NC}"
            echo -e "${YELLOW}      $dest${NC}"
            echo -e "${YELLOW}   Source: $url${NC}"
            echo -e "${YELLOW}   Expected MD5: $expected_md5${NC}"
            echo -e "${YELLOW}   Then re-run with RESUME_AFTER_CHECKPOINT=TRUE bash install-hpc.sh${NC}"
            echo -e "${YELLOW}   See the wiki Installation page, Section 6 (Troubleshooting & FAQ).${NC}"
            return 1
        fi
        ((attempt++))
        sleep $delay
    done
}

echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}        CABS-flex Standalone Automated Installer (Beta)         ${NC}"
echo -e "${BLUE}================================================================${NC}"

# --- Environment Setup ---
module purge
# module load "$GCC_MODULE"
# module load "$BZIP2_MODULE"
# module load "$INTEL_MODULE"
# module load "$HDF5_MODULE"
# module load "$NETCDF_MODULE"
module load "$GCC_MODULE" "$PYTHON_MODULE"


CABS_FLEX_LOCAL_PATH=$(pwd)

if [ ! -f "$CABS_FLEX_LOCAL_PATH/requirements-runtime.txt" ]; then
    echo -e "${RED}❌ Error: requirements-runtime.txt not found in root directory.${NC}"
    exit 1
fi

# Workspace Setup
TEMP_DIR=$(mktemp -d -p "$TEMP_ROOT" "cabs_install_XXXX")
PIP_CACHE_DIR="$TEMP_ROOT/pip-cache"
mkdir -p "$PIP_CACHE_DIR"
cd "$TEMP_DIR"

# --- 1. Preparing CABS Virtual Environment Management ---
# RESUME_AFTER_CHECKPOINT also skips rebuilding this main env, but only if
# it is actually found present and functional (never skipped on the flag
# alone) — used when resuming after a manual checkpoint download: by the
# time that failure happens, the main env is already built, and it is set
# up before the reconstruction one, so redoing it on every retry would be
# pure waste. If the main env is not actually there, it is built normally
# regardless of the flag.
MAIN_ENV_READY=false
if [ "$RESUME_AFTER_CHECKPOINT" = "TRUE" ] && [ -x "$VENV_DIR/bin/CABSflex" ]; then
    echo -e "${BLUE}ℹ️  RESUME_AFTER_CHECKPOINT=TRUE and an existing, functional main environment was found at $VENV_DIR — skipping rebuild.${NC}"
    MAIN_ENV_READY=true
fi

if [ "$MAIN_ENV_READY" = true ]; then
    unset PYTHONPATH
    export PYTHONHOME=""
    export PYTHONUSERBASE=""
    export PYTHONNOUSERSITE=1
    source "$VENV_DIR/bin/activate"
    cd "$CABS_FLEX_LOCAL_PATH"
else
if [ -d "$VENV_DIR" ]; then
    echo -e "${YELLOW}⚠️  Existing venv found at $VENV_DIR. Recreating for clean install...${NC}"
    rm -rf "$VENV_DIR"
fi

unset PYTHONPATH
export PYTHONHOME=""
export PYTHONUSERBASE=""
export PYTHONNOUSERSITE=1
python3 -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"
module purge
module load "$GCC_MODULE"
module load "$BZIP2_MODULE"

# Load site HDF5/netCDF modules (if configured) and point pip's netcdf4/h5py
# builds at them, so they don't fall back to building against an old system
# HDF5 without netCDF-4 support (see requirements-runtime.txt: h5py, netcdf4).
if [ -n "$HDF5_MODULE" ]; then
    module load $HDF5_MODULE
fi
if [ -n "$NETCDF_MODULE" ]; then
    module load $NETCDF_MODULE
fi
export HDF5_DIR="${EBROOTHDF5:-${HDF5_DIR:-}}"
export NETCDF4_DIR="${EBROOTNETCDF:-${NETCDF4_DIR:-}}"

cat >> "$VENV_DIR/bin/activate" <<EOF

# ---- cg2all runtime safety ----
export PYTHONNOUSERSITE=1
export DGL_DISABLE_GRAPHBOLT=1

module purge
EOF

if [ -n "$GCC_MODULE" ]; then
    echo "module load $GCC_MODULE" >> "$VENV_DIR/bin/activate"
fi
if [ -n "$BZIP2_MODULE" ]; then
    echo "module load $BZIP2_MODULE" >> "$VENV_DIR/bin/activate"
fi

cat >> "$VENV_DIR/bin/activate" <<EOF
# --------------------------------

EOF


# --- 2. Install Core Dependencies ---
echo -e "${YELLOW}📦 Upgrading pip tools...${NC}"
pip install --cache-dir "$PIP_CACHE_DIR" --upgrade pip setuptools wheel

echo -e "${YELLOW}📦 Installing core runtime dependencies...${NC}"
cp "$CABS_FLEX_LOCAL_PATH/requirements-runtime.txt" .
pip install --cache-dir "$PIP_CACHE_DIR" -r requirements-runtime.txt


# --- 3. Optional MODELLER Installation (CONDITIONAL) ---
_install_modeller() {
    if [[ "$INSTALL_MODELLER" != "TRUE" ]]; then
        echo -e "${YELLOW}ℹ️  INSTALL_NSP3 is set to FALSE. Skipping ML Prediction setup.${NC}"
        return 0
    fi
    if [ -z "$MODELLER_KEY" ] ; then
        echo -e "${YELLOW}ℹ️  MODELLER_KEY is empty. Skipping Modeller installation.${NC}"
        return 0
    fi

    local install_dir="$VENV_DIR/modeller"
    local arch_name="x86_64-intel8" # Matches Index 2

    echo -e "${YELLOW}📥 Downloading Modeller $MODELLER_VERSION...${NC}"
    curl -L "https://salilab.org/modeller/${MODELLER_VERSION}/modeller-${MODELLER_VERSION}.tar.gz" -o "modeller.tar.gz" --fail
    tar -xzf "modeller.tar.gz"
    cd "modeller-${MODELLER_VERSION}"

    echo -e "${YELLOW}🤖 Running Automated Modeller Installer...${NC}"
    # Pipe the answers (Arch Index, Path, License Key) into the installer
    ./Install <<EOF
$MODELLER_ARCH_INDEX
$install_dir
$MODELLER_KEY


EOF

    echo -e "${YELLOW}🔗 Linking Modeller to Python Environment...${NC}"
    local site_pkgs=$(python3 -c 'import site; print(site.getsitepackages()[0])')

    # Create .pth file to point Python to Modeller libs
    echo "$install_dir/modlib" > "$site_pkgs/modeller.pth"
    echo "$install_dir/lib/$arch_name/python3.3" >> "$site_pkgs/modeller.pth"

    # Set LD_LIBRARY_PATH in the venv activate script
    echo "export LD_LIBRARY_PATH=\"$install_dir/lib/$arch_name:\$LD_LIBRARY_PATH\"" >> "$VENV_DIR/bin/activate"
    export LD_LIBRARY_PATH="$install_dir/lib/$arch_name:$LD_LIBRARY_PATH"

    echo -e "${GREEN}✅ Modeller installation script finished and linked.${NC}"
}
_install_modeller

# # --- 4. MDtraj ----
# echo -e "${YELLOW}📦 Installing mdtraj Package ...${NC}"
# pip install --cache-dir "$PIP_CACHE_DIR" git+https://github.com/mdtraj/mdtraj

# --- 4. CABS-flex Core ---
cd "$CABS_FLEX_LOCAL_PATH"
echo "{\"cg2all_env_prefix\": \"$CG2ALL_VENV_DIR\", \"cabs_env_prefix\": \"$VENV_DIR\"}" > "$CABS_FLEX_LOCAL_PATH/CABS/data/cabs_paths.json"
echo -e "${YELLOW}🧹 Cleaning up build artifacts from $CABS_FLEX_LOCAL_PATH...${NC}"
rm -rf "$CABS_FLEX_LOCAL_PATH/tests/test_cli_options" "$CABS_FLEX_LOCAL_PATH/build" "$CABS_FLEX_LOCAL_PATH/dist" "$CABS_FLEX_LOCAL_PATH"/*.egg-info
echo -e "${YELLOW}📦 Installing CABSflex from local source...${NC}"
pip install --cache-dir "$PIP_CACHE_DIR" .
fi
deactivate

# --- 5. Reconstruction (cg2all) in isolated environment (separate venv) ---
if [ "$RESUME_AFTER_CHECKPOINT" = "TRUE" ] && [ -x "$CG2ALL_VENV_DIR/bin/convert_cg2all" ]; then
    echo -e "${BLUE}ℹ️  RESUME_AFTER_CHECKPOINT=TRUE and an existing reconstruction environment was found at $CG2ALL_VENV_DIR — skipping rebuild.${NC}"
    source "$CG2ALL_VENV_DIR/bin/activate"
else
if [ -d "$CG2ALL_VENV_DIR" ]; then
    echo -e "${YELLOW}⚠️  Existing venv found at $CG2ALL_VENV_DIR. Recreating for clean install...${NC}"
    rm -rf "$CG2ALL_VENV_DIR"
fi

echo -e "${YELLOW}📦 Creating isolated cg2all environment...${NC}"
module purge
module load "$GCC_MODULE" "$PYTHON_MODULE"
unset PYTHONPATH
export PYTHONHOME=""
export PYTHONUSERBASE=""
export PYTHONNOUSERSITE=1
python3 -m venv "$CG2ALL_VENV_DIR"
source "$CG2ALL_VENV_DIR/bin/activate"

cat >> "$CG2ALL_VENV_DIR/bin/activate" <<EOF

# ---- cg2all runtime safety ----
export PYTHONNOUSERSITE=1
export DGL_DISABLE_GRAPHBOLT=1

module purge
EOF

if [ -n "$GCC_MODULE" ]; then
    echo "module load $GCC_MODULE" >> "$CG2ALL_VENV_DIR/bin/activate"
fi
if [ -n "$BZIP2_MODULE" ]; then
    echo "module load $BZIP2_MODULE" >> "$CG2ALL_VENV_DIR/bin/activate"
fi

cat >> "$CG2ALL_VENV_DIR/bin/activate" <<EOF
# --------------------------------

EOF

module purge
module load "$GCC_MODULE"
module load "$BZIP2_MODULE"
pip install --upgrade pip setuptools wheel

echo -e "${YELLOW}📦 Installing dependencies for cg2all package for reconstruction...${NC}"

echo -e "${YELLOW}📦 Installing torch and torchvision ...${NC}"
pip install --cache-dir "$PIP_CACHE_DIR" torch==2.2.0+cpu torchvision==0.17.0+cpu --index-url https://download.pytorch.org/whl/cpu

# pip install " --no-deps torchdata==0.6.1

echo -e "${YELLOW}📦 Installing dgl ...${NC}"
pip install --cache-dir "$PIP_CACHE_DIR" --no-deps dgl==1.1.3 -f https://data.dgl.ai/wheels/repo.html

install_e3nn() {
    local attempt=1 max=5 delay=5
    local pip_args=(--cache-dir "$PIP_CACHE_DIR" --no-binary e3nn e3nn)
    if [ "$E3NN_USE_WHEEL" = "TRUE" ]; then
        echo -e "${BLUE}ℹ️  E3NN_USE_WHEEL=TRUE: installing e3nn from the prebuilt wheel.${NC}"
        pip_args=(--cache-dir "$PIP_CACHE_DIR" e3nn)
    fi
    while true; do
        pip install "${pip_args[@]}" && return 0
        if [[ $attempt -lt $max ]]; then
            ((attempt++))
            echo -e "${YELLOW}⚠️  e3nn install failed. Attempt $attempt/$max. Retrying in $delay seconds...${NC}"
            sleep $delay
        else
            echo -e "${RED}❌ e3nn install failed after $max attempts.${NC}"
            if [ "$E3NN_USE_WHEEL" != "TRUE" ]; then
                echo -e "${YELLOW}ℹ️  By default this installer builds e3nn from source (--no-binary) for compatibility with older cluster kernels.${NC}"
                echo -e "${YELLOW}ℹ️  If this cluster doesn't need that, retry with the prebuilt wheel instead:${NC}"
                echo -e "${YELLOW}      E3NN_USE_WHEEL=TRUE bash install-hpc.sh${NC}"
            fi
            echo -e "${YELLOW}   See the wiki Installation page, Section 6 (Troubleshooting & FAQ).${NC}"
            exit 1
        fi
    done
}
echo -e "${YELLOW}📦 Installing e3nn ...${NC}"
install_e3nn

echo -e "${YELLOW}📦 Installing /huhlim/mdtraj  ...${NC}"
pip install --cache-dir "$PIP_CACHE_DIR" git+https://github.com/huhlim/mdtraj


echo -e "${YELLOW}📦 Installing /huhlim/SE3Transformer  ...${NC}"
SE3T_SRC="$TEMP_DIR/se3t-src"
git clone https://github.com/huhlim/SE3Transformer "$SE3T_SRC"
cd "$SE3T_SRC"
sed -i 's/python = "[^"]*"/python = ">=3.7"/' pyproject.toml
sed -i 's/torch = "[^"]*"/torch = ">=2.1.0"/' pyproject.toml
pip install --cache-dir "$PIP_CACHE_DIR" .

echo -e "${YELLOW}📦 Installing cg2all package for reconstruction...${NC}"

CG2ALL_SRC="$TEMP_DIR/cg2all-src"
git clone https://github.com/huhlim/cg2all.git "$CG2ALL_SRC"
cd "$CG2ALL_SRC"
git checkout a789cb5
sed -i 's/torch = "[^"]*"/torch = ">=2.1.0"/' pyproject.toml
sed -i 's/numpy = "[^"]1"/numpy = ">=1.21"/' pyproject.toml
pip install --cache-dir "$PIP_CACHE_DIR" --no-binary :all: .
fi

# --- 5b. Verify cg2all checkpoint files (see CABS/utils/utils.py
# CG2ALL_REPRESENTATIONS for the models CABS-flex actually uses) ---
echo -e "${YELLOW}🔍 Verifying cg2all checkpoint files...${NC}"
MODEL_HOME=$(python3 -c "import cg2all.lib.libconfig as c; print(c.MODEL_HOME)")
CKPT_OK=true
fetch_checkpoint "CalphaBasedModel.ckpt" "0b51f6fe4a12c878ec28b194e55099d3" "$MODEL_HOME" || CKPT_OK=false
fetch_checkpoint "CalphaSCModel.ckpt"    "d42f297f94b4ea33dafa4145b6495344" "$MODEL_HOME" || CKPT_OK=false
if [ "$CKPT_OK" = false ]; then
    echo -e "${RED}❌ Installation stopped: cg2all checkpoint verification failed (see instructions above).${NC}"
    exit 1
fi

deactivate

# --- 6. Final Verification ---
test_binary() {
    local bin_path="$1"
    local bin_name=$(basename "$bin_path")
    if [ -x "$bin_path" ]; then
        if "$bin_path" --help > /dev/null 2>&1; then
            echo -e "${GREEN}✅ $bin_name confirmed and functional.${NC}"
        else
            echo -e "${RED}❌ $bin_name exists but failed execution.${NC}"
        fi
    else
        echo -e "${RED}❌ $bin_name binary not found.${NC}"
    fi
}

echo -e "${YELLOW}🧪 Verifying Environment...${NC}"
source "$VENV_DIR/bin/activate"
python3 <<EOF
try:
    import Bio.PDB
    print(f"${GREEN}✅ BioPython dependencies ready.${NC}")
except ImportError:
    print("${RED}⚠️  BioPython module not found.${NC}")
try:
    import mdtraj
    print(f"${GREEN}✅ Mdtraj dependencies ready.${NC}")
except ImportError:
    print("${RED}⚠️  Mdtraj module not found.${NC}")

# Check for Modeller linkage (needs its lib path set)
try:
    import modeller
    print("${GREEN}✅ Modeller linked successfully.${NC}")
except ImportError:
    print("${YELLOW}⚠️  Modeller Python module not found/linked.${NC}")

EOF

deactivate

source "$CG2ALL_VENV_DIR/bin/activate"
python3 <<EOF
try:
    import cg2all
    print(f"${GREEN}✅ cg2all module is ready.${NC}")
except ImportError:
    print("${RED}⚠️  cg2all module not found.${NC}")

EOF
deactivate

echo -e "${BLUE}Checking Main Environment binaries:${NC}"
test_binary "$VENV_DIR/bin/CABSflex"
test_binary "$VENV_DIR/bin/CABSdock"

echo -e "${BLUE}Checking Reconstruction Environment binaries:${NC}"
test_binary "$CG2ALL_VENV_DIR/bin/convert_cg2all"

rm -rf "$TEMP_DIR"

echo -e "${GREEN}🎉 CABS-flex installation complete!${NC}"
echo "============================================================"
echo -e "${BLUE}To start:${NC} source $VENV_DIR/bin/activate"
echo "============================================================"
