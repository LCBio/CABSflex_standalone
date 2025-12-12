#!/bin/bash
# CABS-flex Local Installer (Non-Conda, Runs from CABSflex Root, Auto-loads Modules, HPC-Optimized Paths)
# Usage: Run this script from the root directory of your cloned CABSflex repository.
# Example: cd /path/to/your/cabsflex/clone && ./install_beta_helios.sh

set -e

# --- Installation Paths Configuration ---
BASE_INSTALL_DIR="${PLG_GROUPS_STORAGE}/plggmodel/NC/programs"
VENV_NAME="cabs-venv"
VENV_DIR="$BASE_INSTALL_DIR/$VENV_NAME"
TEMP_ROOT="$SCRATCH"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}🧬 CABS-flex Local Installer (Non-Conda, HPC Paths)${NC}"
echo "================================================================"

# --- Environment Checks ---
echo -e "${YELLOW}📋 Checking environment variables...${NC}"
if [ -z "$PLG_GROUPS_STORAGE" ] || [ -z "$SCRATCH" ]; then
    echo -e "${RED}❌ Error: \$PLG_GROUPS_STORAGE and \$SCRATCH must be set.${NC}"
    exit 1
fi
echo -e "${GREEN}✅ Environment variables OK.${NC}"
echo ""

# --- Environment Module Loading ---
echo -e "${YELLOW}⚙️  Loading required environment modules (GCCcore, Python)...${NC}"
module purge
module load GCCcore/13.2.0 || { echo -e "${RED}❌ Failed to load GCCcore/13.2.0 module. Check 'module spider GCCcore' for available versions.${NC}"; exit 1; }
module load Python/3.11.5 || { echo -e "${RED}❌ Failed to load Python/3.11.5 module. Check 'module spider Python/3.11.5' for its specific dependencies.${NC}"; exit 1; }
echo -e "${GREEN}✅ Modules loaded.${NC}"
echo ""

# --- START: Enforce Environment Isolation ---
echo -e "${YELLOW}🛡️  Enforcing a clean, isolated Python environment...${NC}"
export PYTHONNOUSERSITE=1
echo -e "${GREEN}✅ User site-packages will be ignored to prevent contamination.${NC}"
echo ""


# --- Prerequisite Checks ---
CABS_FLEX_LOCAL_PATH=$(pwd)
echo -e "${YELLOW}📋 Checking core prerequisites...${NC}"
PYTHON_CMD=$(command -v python3)
if [ -z "$PYTHON_CMD" ] || ! "$PYTHON_CMD" -c "import sys; exit(not (sys.version_info.major == 3 and sys.version_info.minor >= 10))"; then
    echo -e "${RED}❌ Error: Python 3.10+ not found. Module load may have failed.${NC}"; exit 1;
fi
if [ ! -f "$CABS_FLEX_LOCAL_PATH/environment.yml" ]; then
    echo -e "${RED}❌ Error: 'environment.yml' not found. This script must be run from the CABSflex clone root.${NC}"; exit 1;
fi
echo -e "${GREEN}✅ Prerequisites OK (Using $PYTHON_CMD)${NC}"
echo ""

# Create base installation directory if it doesn't exist
echo -e "${YELLOW}📂 Ensuring base installation directory exists: $BASE_INSTALL_DIR${NC}"
mkdir -p "$BASE_INSTALL_DIR" || { echo -e "${RED}❌ Error: Failed to create base installation directory '$BASE_INSTALL_DIR'. Check permissions or path.${NC}"; exit 1; }

# --- Setup Directories ---
TEMP_DIR=$(mktemp -d -p "$TEMP_ROOT" "cabsflex-install-XXXXXXXX")
export TMPDIR=$TEMP_DIR
export PIP_BUILD=$TEMP_DIR
export PIP_TMPDIR=$TEMP_DIR
PIP_CACHE_DIR="$TEMP_ROOT/pip-cache"
mkdir -p $PIP_CACHE_DIR
cd "$TEMP_DIR"
echo -e "${YELLOW}📂 Working in temporary directory: $TEMP_DIR${NC}"

# Check if virtual environment already exists
if [ -d "$VENV_DIR" ]; then
    echo -e "${YELLOW}⚠️  Virtual environment '${VENV_NAME}' already exists at $VENV_DIR${NC}"
    read -p "Do you want to REMOVE and recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo -e "${YELLOW}🗑️  Removing existing virtual environment...${NC}"
        rm -rf "$VENV_DIR"
    else
        echo -e "${BLUE}ℹ️  Installation aborted by user.${NC}"
        exit 0
    fi
fi

# Create virtual environment
echo -e "${YELLOW}🔧 Creating new virtual environment at $VENV_DIR...${NC}"
"$PYTHON_CMD" -m venv "$VENV_DIR" || { echo -e "${RED}❌ Failed to create virtual environment.${NC}"; exit 1; }
source "$VENV_DIR/bin/activate"

# Install core pip packages
echo -e "${YELLOW}📦 Installing core pip dependencies...${NC}"
echo "This may take a few minutes..."
pip install --cache-dir "$PIP_CACHE_DIR" --upgrade pip setuptools wheel || { echo -e "${RED}❌ Error upgrading pip tools.${NC}"; deactivate; exit 1; }

echo -e "${YELLOW}📦 Installing core runtime dependencies from requirements-runtime.txt...${NC}"
# We copy the new requirements file from the clone root
cp "$CABS_FLEX_LOCAL_PATH/requirements-runtime.txt" .
pip install -r requirements-runtime.txt || { echo -e "${RED}❌ Error installing from requirements-runtime.txt.${NC}"; deactivate; exit 1; }

# --- START: Simplified and Corrected ML Package Installation ---
echo -e "${YELLOW}📦 Installing cg2all and its dependencies...${NC}"
pip install --cache-dir "$PIP_CACHE_DIR" torch torchvision --index-url https://download.pytorch.org/whl/cpu
pip install --cache-dir "$PIP_CACHE_DIR" git+http://github.com/huhlim/cg2all

# --- DSSP Installation Function (Modified to use pre-existing files) ---
_install_dssp() {
    local install_prefix="$1" # This is the VENV_DIR
    local bin_dest="${install_prefix}/bin/"

    echo -e "${YELLOW}⚙️  Searching for your pre-existing DSSP binary and library...${NC}"

    # 1. Handle the 'mkdssp' executable
    local mkdssp_path
    mkdssp_path=$(command -v mkdssp)

    if [ -n "$mkdssp_path" ]; then
        echo -e "${GREEN}✅ Found 'mkdssp' executable at: $mkdssp_path${NC}"
        echo -e "${YELLOW}   -> Copying to $bin_dest...${NC}"
        cp "$mkdssp_path" "$bin_dest" || { echo -e "${RED}❌ Error copying mkdssp executable.${NC}"; return 1; }
        chmod +x "${bin_dest}/mkdssp"
    else
        echo -e "${YELLOW}⚠️  Warning: 'mkdssp' executable not found in your PATH.${NC}"
        echo -e "${YELLOW}   The installation will continue, but CABSflex may fail if it needs this command.${NC}"
    fi

    echo -e "${GREEN}✅ DSSP check complete.${NC}"
    return 0
}

# Call DSSP build function
_install_dssp "$VENV_DIR" || { echo -e "${RED}❌ DSSP installation failed. See messages above.${NC}"; deactivate; exit 1; }

# --- Modeller 10.7 Installation Function (User-assisted) ---
_install_modeller() {
    echo -e "${YELLOW}⚙️  Attempting to install Modeller 10.7 (User-Assisted)...${NC}"
    echo -e "${RED}🚨 IMPORTANT: Modeller requires a license key and user interaction!${NC}"
    echo -e "${YELLOW}Please visit https://salilab.org/modeller/registration.html to obtain a license key BEFORE proceeding.${NC}"
    echo ""

    local modeller_version="10.7"
    local modeller_tarball="modeller-${modeller_version}.tar.gz"
    local modeller_url="https://salilab.org/modeller/${modeller_version}/$modeller_tarball"
    local modeller_src_dir="modeller-${modeller_version}"
    local suggested_install_dir="$VENV_DIR/$modeller_src_dir" # Suggest installing into venv dir

    # Check if Modeller is already extracted in temp dir for update scenario
    if [ -d "$modeller_src_dir" ]; then
        echo -e "${YELLOW}ℹ️  Modeller source directory already exists, skipping download and unpack.${NC}"
        cd "$modeller_src_dir"
    else
        # Download
        echo -e "${YELLOW}📥 Downloading Modeller from $modeller_url...${NC}"
        if ! curl -L "$modeller_url" -o "$modeller_tarball" --fail; then
            echo -e "${RED}❌ Error: Failed to download Modeller. Check URL or network.${NC}"
            return 1
        fi

        # Unpack
        echo -e "${YELLOW}📦 Unpacking Modeller...${NC}"
        if ! gunzip "$modeller_tarball"; then
            echo -e "${RED}❌ Error: Failed to gunzip Modeller tarball.${NC}"
            return 1
        fi
        local tar_file="${modeller_tarball%.gz}" # remove .gz
        if ! tar -xvf "$tar_file"; then
            echo -e "${RED}❌ Error: Failed to untar Modeller archive.${NC}"
            return 1
        fi

        if [ ! -d "$modeller_src_dir" ]; then
            echo -e "${RED}❌ Error: Modeller source directory '$modeller_src_dir' not found after unpacking.${NC}"
            return 1
        fi

        cd "$modeller_src_dir"
    fi

    echo -e "${BLUE}================================================================${NC}"
    echo -e "${BLUE}Starting Modeller's interactive installation script.${NC}"
    echo -e "${BLUE}When prompted for the installation directory, we recommend (but you can choose another):${NC}"
    echo -e "${BLUE}  $suggested_install_dir${NC}"
    echo -e "${BLUE}You WILL be asked for your Modeller license key.${NC}"
    echo -e "${BLUE}Press Enter to start the Modeller installer...${NC}"
    read -r # Wait for user to read

    # Run the interactive installer
    if ! ./Install; then
        echo -e "${RED}❌ Modeller interactive installation failed.${NC}"
        echo -e "${RED}Please check errors during the installation process, especially the license key and chosen directory.${NC}"
        cd .. # Go back to TEMP_DIR
        return 1
    fi

    cd .. # Go back to TEMP_DIR (from modeller_src_dir)
    echo -e "${GREEN}✅ Modeller installation script finished. Please verify its installation manually.${NC}"
    echo -e "${BLUE}Remember the installation path you chose (e.g., '$suggested_install_dir').${NC}"
    return 0
}

# Call Modeller installation function
_install_modeller || { echo -e "${RED}❌ Modeller installation failed or was aborted. See messages above.${NC}"; deactivate; exit 1; }


# Install CABSflex from local repository
echo -e "${YELLOW}📦 Installing CABSflex from local repository '$CABS_FLEX_LOCAL_PATH'...${NC}"
if ! pip install --cache-dir "$PIP_CACHE_DIR" "$CABS_FLEX_LOCAL_PATH"; then
    echo -e "${RED}❌ Error: Failed to install CABSflex from local path.${NC}"
    deactivate
    exit 1
fi

# --- Final Tests ---
echo -e "${YELLOW}🧪 Testing installation...${NC}"
if ! "$VENV_DIR/bin/CABSflex" --help > /dev/null 2>&1 || ! "$VENV_DIR/bin/CABSdock" --help > /dev/null 2>&1; then
    echo -e "${RED}❌ Warning: CABSflex or CABSdock command not working.${NC}"
else
    echo -e "${GREEN}✅ CABSflex and CABSdock commands working!${NC}"
fi
if ! "$VENV_DIR/bin/mkdssp" --version > /dev/null 2>&1; then
    echo -e "${GREEN}✅ mkdssp command working!${NC}"
else
    echo -e "${RED}❌ Warning: mkdssp command not working.${NC}"
fi

deactivate
cd - > /dev/null
rm -rf "$TEMP_DIR"

echo ""
echo -e "${GREEN}🎉 CABS-flex Local installation complete!${NC}"
echo "============================================================"
echo -e "${BLUE}To use CABS-flex:${NC}"
echo "  source $VENV_DIR/bin/activate"
echo "  # Remember to add your Modeller 'bin' directory to the PATH if needed."
echo "  CABSflex --help"
echo ""
echo -e "${BLUE}To update, simply delete the '$VENV_DIR' directory and run this script again.${NC}"
echo ""
