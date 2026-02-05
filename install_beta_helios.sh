#!/bin/bash
# CABS-flex Local Installer (Beta - Helios Optimized)
# Replaces binary DSSP with MDTraj. Maintains interactive Modeller setup.
set -e

# ==============================================================================
# --- 1. CONFIGURATION SECTION (Cluster-Specific) ---
# ==============================================================================
# Installation Paths
BASE_INSTALL_DIR="${PLG_GROUPS_STORAGE}/plggmodel/NC/programs"
VENV_NAME="cabs_005"
VENV_DIR="$BASE_INSTALL_DIR/$VENV_NAME"
TEMP_ROOT="$SCRATCH"

# Module Names
GCC_MODULE="GCCcore/13.2.0"
PYTHON_MODULE="Python/3.11.5"

# Dependency Lists
# Standard scientific tools
CORE_DEPS=("numpy" "matplotlib" "requests" "biopython" "mdtraj" "biopandas" "tqdm")

# ML Package Hardware Configuration (HPC nodes often lack GPUs)
TORCH_URL="https://download.pytorch.org/whl/cpu"

# Modeller Version
MODELLER_VERSION="10.7"
# ==============================================================================

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}🧬 CABS-flex Local Installer (Beta - Helios Optimized)${NC}"
echo "================================================================"

# --- Prerequisite & Environment Checks ---
echo -e "${YELLOW}📋 Checking HPC environment...${NC}"
if [ -z "$PLG_GROUPS_STORAGE" ] || [ -z "$SCRATCH" ]; then
    echo -e "${RED}❌ Error: \$PLG_GROUPS_STORAGE and \$SCRATCH must be set on Helios.${NC}"
    exit 1
fi

# --- Module Loading ---
echo -e "${YELLOW}⚙️  Loading environment modules ($GCC_MODULE, $PYTHON_MODULE)...${NC}"
module purge
module load "$GCC_MODULE"
module load "$PYTHON_MODULE"

# --- Enforce Isolation ---
export PYTHONNOUSERSITE=1
CABS_FLEX_LOCAL_PATH=$(pwd)

# --- Setup Directories ---
mkdir -p "$BASE_INSTALL_DIR"
TEMP_DIR=$(mktemp -d -p "$TEMP_ROOT" "cabsflex-install-XXXXXXXX")
export TMPDIR=$TEMP_DIR
export PIP_BUILD=$TEMP_DIR
export PIP_TMPDIR=$TEMP_DIR
PIP_CACHE_DIR="$TEMP_ROOT/pip-cache"
mkdir -p "$PIP_CACHE_DIR"

echo -e "${YELLOW}📂 Working in temporary directory: $TEMP_DIR${NC}"

# --- Virtual Environment Management ---
if [ -d "$VENV_DIR" ]; then
    echo -e "${YELLOW}⚠️  Venv '${VENV_NAME}' already exists at $VENV_DIR${NC}"
    read -p "Do you want to REMOVE and recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf "$VENV_DIR"
    else
        echo -e "${BLUE}ℹ️  Aborted by user.${NC}"
        exit 0
    fi
fi

echo -e "${YELLOW}🔧 Creating new virtual environment at $VENV_DIR...${NC}"
python3 -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"

# --- 1. Core Pip Dependencies ---
echo -e "${YELLOW}📦 Upgrading pip tools...${NC}"
pip install --cache-dir "$PIP_CACHE_DIR" --upgrade pip setuptools wheel

echo -e "${YELLOW}📦 Installing CABS runtime dependencies (Pure Python)...${NC}"
pip install --cache-dir "$PIP_CACHE_DIR" "${CORE_DEPS[@]}"

# --- 2. ML & Reconstruction (cg2all) ---
echo -e "${YELLOW}📦 Installing ML Package (cg2all) and PyTorch (CPU)...${NC}"
pip install --cache-dir "$PIP_CACHE_DIR" torch torchvision --index-url "$TORCH_URL"
pip install --cache-dir "$PIP_CACHE_DIR" git+http://github.com/huhlim/cg2all

# --- 3. Modeller 10.7 (Interactive Setup) ---
_install_modeller() {
    echo -e "${YELLOW}⚙️  Starting Modeller $MODELLER_VERSION Interactive Setup...${NC}"
    echo -e "${RED}🚨 IMPORTANT: Modeller requires a license key!${NC}"
    echo -e "Register at: https://salilab.org/modeller/registration.html"
    echo ""

    local tarball="modeller-${MODELLER_VERSION}.tar.gz"
    local url="https://salilab.org/modeller/${MODELLER_VERSION}/$tarball"
    local src_dir="modeller-${MODELLER_VERSION}"
    local suggested_dir="$VENV_DIR/modeller"

    cd "$TEMP_DIR"
    echo -e "${YELLOW}📥 Downloading Modeller...${NC}"
    curl -L "$url" -o "$tarball" --fail

    echo -e "${YELLOW}📦 Unpacking...${NC}"
    tar -xzf "$tarball"
    cd "$src_dir"

    echo -e "${BLUE}================================================================${NC}"
    echo -e "Starting Modeller installer script."
    echo -e "Recommended installation path: $suggested_dir"
    echo -e "================================================================${NC}"

    # Run the standard installer provided by SaliLab
    if ! ./Install; then
        echo -e "${RED}❌ Modeller interactive installation failed.${NC}"
        return 1
    fi
    return 0
}

# Modeller is usually required for CABS v3 final reconstruction
_install_modeller || { echo -e "${RED}❌ Modeller installation failed.${NC}"; deactivate; exit 1; }

# --- 4. CABS-flex Local Install (Dev Mode) ---
cd "$CABS_FLEX_LOCAL_PATH"
echo -e "${YELLOW}📦 Installing CABS-flex from local source...${NC}"
pip install --cache-dir "$PIP_CACHE_DIR" -e .

# --- Final Tests ---
echo -e "${YELLOW}🧪 Verifying Installation...${NC}"
if CABSflex --version && CABSdock --version; then
    echo -e "${GREEN}✅ CABS-flex commands working!${NC}"
    if python3 -c "import mdtraj; print('MDTraj OK')" &> /dev/null; then
        echo -e "${GREEN}✅ Secondary Structure Engine (MDTraj) OK.${NC}"
    fi
else
    echo -e "${RED}❌ Warning: Commands not responding. Check your PATH.${NC}"
fi

# Cleanup
deactivate
rm -rf "$TEMP_DIR"

echo ""
echo -e "${GREEN}🎉 CABS-flex Local Beta installation complete!${NC}"
echo "============================================================"
echo -e "${BLUE}To use CABS-flex:${NC}"
echo "  source $VENV_DIR/bin/activate"
echo "  CABSflex --help"
echo ""
echo -e "${YELLOW}Note: Secondary structure is now 100% Python-based (MDTraj).${NC}"
echo -e "${YELLOW}No legacy 'mkdssp' binary is required for this installation.${NC}"
echo "============================================================"
