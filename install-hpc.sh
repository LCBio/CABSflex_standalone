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
VENV_NAME="cabs_021"
VENV_DIR="$BASE_INSTALL_DIR/$VENV_NAME"
CG2ALL_VENV_DIR="$VENV_DIR/cg2all"
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
# ==============================================================================

# Output Colors
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

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
deactivate

# --- 5. Reconstruction (cg2all) in isolated environment (separate venv) ---
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
pip install --cache-dir "$PIP_CACHE_DIR" torch==2.1.2+cpu torchvision==0.16.2+cpu --index-url https://download.pytorch.org/whl/cpu

# pip install " --no-deps torchdata==0.6.1

echo -e "${YELLOW}📦 Installing dgl ...${NC}"
pip install --cache-dir "$PIP_CACHE_DIR" --no-deps dgl==1.1.3 -f https://data.dgl.ai/wheels/repo.html

echo -e "${YELLOW}📦 Installing e3nn ...${NC}"
pip install --cache-dir "$PIP_CACHE_DIR" --no-binary e3nn e3nn

echo -e "${YELLOW}📦 Installing /huhlim/mdtraj  ...${NC}"
pip install --cache-dir "$PIP_CACHE_DIR" git+https://github.com/huhlim/mdtraj


echo -e "${YELLOW}📦 Installing /huhlim/SE3Transformer  ...${NC}"
SE3T_SRC="$TEMP_DIR/se3t-src"
git clone https://github.com/huhlim/SE3Transformer "$SE3T_SRC"
cd "$SE3T_SRC"
sed -i 's/python = "[^"]*"/python = ">=3.7"/' pyproject.toml
sed -i 's/torch = "[^"]*"/torch = "=2.1.2"/' pyproject.toml
pip install --cache-dir "$PIP_CACHE_DIR" .

echo -e "${YELLOW}📦 Installing cg2all package for reconstruction...${NC}"

CG2ALL_SRC="$TEMP_DIR/cg2all-src"
git clone https://github.com/huhlim/cg2all.git "$CG2ALL_SRC"
cd "$CG2ALL_SRC"
git checkout a789cb5
sed -i 's/torch = "[^"]*"/torch = "=2.1.2"/' pyproject.toml
sed -i 's/numpy = "[^"]1"/numpy = ">=1.21"/' pyproject.toml
pip install --cache-dir "$PIP_CACHE_DIR" --no-binary :all: .

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
