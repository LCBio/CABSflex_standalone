#!/bin/bash
# ==============================================================================
# 🧬 CABS-flex Local Installer (Beta - Helios Optimized)
# ==============================================================================
set -e

# ==============================================================================
# --- 1. CONFIGURATION SECTION ---
# ==============================================================================
MODELLER_KEY="MODELIRANJE"   # <<< SET YOUR ACADEMIC LICENSE KEY HERE
MODELLER_VERSION="10.7"
MODELLER_ARCH_INDEX="2"      # 2 = x86_64-intel8

INSTALL_NSP3="FALSE"         # <<< SET TO TRUE TO ENABLE FULL ML TIER INSTALLATION

BASE_INSTALL_DIR="${PLG_GROUPS_STORAGE}/plggmodel/NC/programs"
VENV_NAME="cabs_008"
VENV_DIR="$BASE_INSTALL_DIR/$VENV_NAME"
TEMP_ROOT="$SCRATCH" # Use scratch for high-I/O operations

# HPC Environment Modules
GCC_MODULE="GCCcore/13.2.0"
PYTHON_MODULE="Python/3.11.5"

# Dependencies (Kept lean: MDTraj replaces DSSP binary)
CORE_DEPS=("numpy" "matplotlib" "requests" "biopython" "mdtraj" "biopandas" "tqdm" "scipy")
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
module load "$GCC_MODULE" "$PYTHON_MODULE"
export PYTHONNOUSERSITE=1
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

# --- Virtual Environment Management ---
if [ -d "$VENV_DIR" ]; then
    echo -e "${YELLOW}⚠️  Existing venv found at $VENV_DIR. Recreating for clean install...${NC}"
    rm -rf "$VENV_DIR"
fi
python3 -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"

# --- 1. Install Core Dependencies ---
echo -e "${YELLOW}📦 Upgrading pip tools...${NC}"
pip install --cache-dir "$PIP_CACHE_DIR" --upgrade pip setuptools wheel

echo -e "${YELLOW}📦 Installing core runtime dependencies...${NC}"
cp "$CABS_FLEX_LOCAL_PATH/requirements-runtime.txt" .
pip install --cache-dir "$PIP_CACHE_DIR" -r requirements-runtime.txt


# --- 2. OPTIONAL: NetSurfP-3.0 Package & Model Download --
_install_nsp3() {
    if [[ "$INSTALL_NSP3" != "TRUE" ]]; then
        echo -e "${YELLOW}ℹ️  INSTALL_NSP3 is set to FALSE. Skipping ML Prediction setup.${NC}"
        return 0
    fi

    local nsp3_dir="$TEMP_DIR/nsp3_repo"
    local nsp3_model_dir="$CABS_FLEX_LOCAL_PATH/nsp3_model"

    echo -e "${YELLOW}⚙️  Setting up NetSurfP-3.0 (NSP3) Tier 1 Prediction...${NC}"

    # Clone the repository source code
    if [ ! -d "$nsp3_dir" ]; then
        echo -e "${YELLOW}📥 Cloning NetSurfP-3.0 source code...${NC}"
        git clone "$NSP3_REPO_URL" "$nsp3_dir"
    fi

    # Check for model URLs file
    local url_txt="$nsp3_dir/models/url.txt"
    if [ ! -f "$url_txt" ]; then
        echo -e "${RED}❌ Error: NSP3 URL list not found at $url_txt. Cannot download model weights.${NC}"
        return 1
    fi

    # Download models based on url.txt content
    mkdir -p "$nsp3_model_dir"
    echo -e "${YELLOW}⚙️  Downloading NSP3 model weights...${NC}"
    while IFS= read -r line; do
        # Assumes format: [NAME]: [URL]
        local url=$(echo "$line" | awk -F': ' '{print $2}')
        if [[ "$url" == *".pt"* ]]; then
            local filename=$(basename "$url")
            echo "   -> Downloading $filename..."
            curl -L -o "$nsp3_model_dir/$filename" "$url" --fail
        fi
    done < "$url_txt"
    echo -e "${GREEN}✅ NSP3 Model Weights downloaded to $nsp3_model_dir.${NC}"

    # Install the NSP3 package itself
    echo -e "${YELLOW}📦 Installing NetSurfP-3.0 package...${NC}"
    pip install --cache-dir "$PIP_CACHE_DIR" nsp3
    echo -e "${GREEN}✅ NSP3 package installed.${NC}"
}

_install_nsp3

# --- 2. ML Reconstruction (cg2all) ---
echo -e "${YELLOW}📦 Installing ML Package (cg2all) and PyTorch (CPU)...${NC}"
# pip install --cache-dir "$PIP_CACHE_DIR" torch torchvision --index-url "$TORCH_URL"
pip install --cache-dir "$PIP_CACHE_DIR" git+http://github.com/huhlim/cg2all

# --- 3. MODELLER Installation (CONDITIONAL) ---
_install_modeller() {
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

# --- 4. CABS-flex Core (Editable Mode) ---
cd "$CABS_FLEX_LOCAL_PATH"
echo -e "${YELLOW}📦 Installing CABSflex from local source...${NC}"
pip install --cache-dir "$PIP_CACHE_DIR" -e .

# --- 5. Final Verification ---
echo -e "${YELLOW}🧪 Verifying Environment...${NC}"
python3 <<EOF
import sys, torch, mdtraj, Bio.PDB
print(f"${GREEN}✅ Core Python dependencies ready.${NC}")
if hasattr(torch, 'compiler'):
    print(f"${GREEN}✅ ML Tier (Torch {torch.__version__}) is available.${NC}")
else:
    print("${YELLOW}⚠️  ML Tier (Torch) is missing. Structure reconstruction will use CG2ALL/Heuristics.${NC}")

# Check for Modeller linkage (needs its lib path set)
try:
    import modeller
    print("${GREEN}✅ Modeller linked successfully.${NC}")
except ImportError:
    print("${YELLOW}⚠️  Modeller Python module not found/linked.${NC}")

EOF

deactivate
rm -rf "$TEMP_DIR"

echo -e "${GREEN}🎉 CABS-flex installation complete!${NC}"
echo "============================================================"
echo -e "${BLUE}To start:${NC} source $VENV_DIR/bin/activate"
echo "============================================================"
