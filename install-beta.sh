#!/bin/bash
# CABS-flex Beta Installer
# Usage: curl -sSL https://raw.githubusercontent.com/LCBio/cabsflex/main/install-beta.sh | bash -s YOUR_BETA_TOKEN

set -e

BETA_TOKEN=$1
REPO_URL="https://github.com/LCBio/cabsflex"
ENV_NAME="cabs"
SOURCE_DIR=$(pwd)
IS_LOCAL=false

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}🧬 CABS-flex Beta Installer${NC}"
echo "================================="

if [ -z "$BETA_TOKEN" ]; then
    echo -e "${YELLOW}⚠️  No beta token provided. Assuming local installation from current directory.${NC}"
    IS_LOCAL=true
    
    # Check if we are in a valid repo
    if [ ! -f "$SOURCE_DIR/pyproject.toml" ] && [ ! -f "$SOURCE_DIR/setup.py" ]; then
        echo -e "${RED}❌ Error: No pyproject.toml or setup.py found in $SOURCE_DIR${NC}"
        echo "For local installation, please run this script from the root of the CABSflex repository."
        exit 1
    fi
    echo -e "${GREEN}✅ Local source found: $SOURCE_DIR${NC}"
else
    echo -e "${GREEN}✅ Beta token provided${NC}"
fi

echo -e "${YELLOW}📋 Checking prerequisites...${NC}"

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo -e "${RED}❌ Error: conda not found${NC}"
    echo ""
    echo "Please install Anaconda or Miniconda first:"
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Check if git is available
if ! command -v git &> /dev/null; then
    echo -e "${RED}❌ Error: git not found${NC}"
    echo "Please install git first."
    exit 1
fi

echo -e "${GREEN}✅ Prerequisites OK${NC}"

# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------
# Power user configuration
INSTALL_MODELLER="TRUE"
MODELLER_KEY=""   # <<< SET YOUR MODELLER LICENSE KEY HERE

# Ask for Modeller Key if we want to install it
if [ "$INSTALL_MODELLER" = "TRUE" ] && [ -z "$MODELLER_KEY" ]; then
    echo ""
    echo -e "${YELLOW}🔑 Modeller License Key (Optional)${NC}"
    echo "CABS-flex uses Modeller for reconstruction."
    echo "If you have a license key, enter just the key."
    echo "Press ENTER to skip Modeller configuration."
    read -p "Key: " INPUT_KEY
    if [ ! -z "$INPUT_KEY" ]; then
        MODELLER_KEY=$INPUT_KEY
    else
        echo -e "${YELLOW}⚠️  Skipping Modeller configuration (no key provided).${NC}"
        INSTALL_MODELLER="FALSE"
    fi
fi

# ---------------------------------------------------------
# Installation
# ---------------------------------------------------------

# Create temporary directory for building/downloading deps
TEMP_DIR=$(mktemp -d)
echo -e "${YELLOW}📂 Working in temp dir: $TEMP_DIR${NC}"

if [ "$IS_LOCAL" = true ]; then
    # Local Install Logic
    cp "$SOURCE_DIR/environment.yml" "$TEMP_DIR/environment.yml"
    cd "$TEMP_DIR"
else
    # Remote Install Logic
    cd "$TEMP_DIR"
    # Test token by trying to access the repo
    echo -e "${YELLOW}🔐 Validating beta token...${NC}"
    if ! git ls-remote https://$BETA_TOKEN@github.com/LCBio/cabsflex.git &> /dev/null; then
        echo -e "${RED}❌ Error: Invalid beta token or no repository access${NC}"
        echo "Please check your token or contact k.wroblewski7@uw.edu.pl"
        exit 1
    fi
    echo -e "${GREEN}✅ Token validated${NC}"

    # Download environment file using GitHub API
    echo -e "${YELLOW}📥 Downloading environment configuration...${NC}"
    if ! curl -H "Authorization: token $BETA_TOKEN" \
         -H "Accept: application/vnd.github.v3.raw" \
         -L "https://api.github.com/repos/LCBio/cabsflex/contents/environment.yml" \
         -o environment.yml --fail; then
        echo -e "${RED}❌ Error: Failed to download environment file${NC}"
        exit 1
    fi
fi

# Filter environment.yml to remove things we will install separately
# We remove cg2all, dgl, fair-esm (likely cg2all dep), and ensure clean install
sed -i '/dgl/d' environment.yml
sed -i '/cg2all/d' environment.yml
sed -i '/fair-esm/d' environment.yml

# Create environment
echo -e "${YELLOW}🔧 Creating/Updating conda environment '${ENV_NAME}'...${NC}"
if conda env list | grep -q "^${ENV_NAME} "; then
    echo -e "${BLUE}ℹ️  Updating existing environment...${NC}"
    conda env update -f environment.yml -n $ENV_NAME
else
    echo -e "${BLUE}ℹ️  Creating new environment...${NC}"
    conda env create -f environment.yml -n $ENV_NAME
fi

# Activate environment
eval "$(conda shell.bash hook)"
conda activate $ENV_NAME

# Install CABSflex
echo -e "${YELLOW}📦 Installing CABSflex...${NC}"

if [ "$IS_LOCAL" = true ]; then
    # Install from local source
    pip install --upgrade "$SOURCE_DIR"
else
    # Install from git
    pip install --upgrade git+https://$BETA_TOKEN@github.com/LCBio/cabsflex.git
fi

# Configure Modeller if key provided
if [ ! -z "$MODELLER_KEY" ]; then
    echo -e "${YELLOW}🔧 Configuring Modeller...${NC}"
    # Conda install of modeller (from salilab channel in env.yml) usually lacks config
    # We find the config.py and patch it
    MOD_CONFIG=$(python -c "import modlib.modeller.config as c; print(c.__file__)" 2>/dev/null || echo "")
    if [ ! -z "$MOD_CONFIG" ]; then
        # Replace license line
        sed -i "s/license = '.*'/license = '${MODELLER_KEY}'/" "$MOD_CONFIG"
        echo -e "${GREEN}✅ Modeller license configured${NC}"
    else
        echo -e "${RED}⚠️  Could not find Modeller config file. Please check installation.${NC}"
    fi
fi

# ---------------------------------------------------------
# Setup CG2ALL (Reconstruction) in Separate Environment
# ---------------------------------------------------------
echo -e "${YELLOW}🔧 Setting up Reconstruction Environment (cg2all)...${NC}"
CG2ALL_ENV_NAME="${ENV_NAME}_reconstruct"
# Create separate env for reconstruction to isolate dependencies
conda create -y -n $CG2ALL_ENV_NAME python=3.9 2>/dev/null || true

# Install dependencies in separate env
eval "$(conda shell.bash hook)"
conda activate $CG2ALL_ENV_NAME

echo -e "${BLUE}Installing PyTorch and DGL for reconstruction...${NC}"
# Use pip for specific versions as per Helios script
# Install Torch 2.1.2 (CPU)
pip install torch==2.1.2+cpu torchvision==0.16.2+cpu --index-url https://download.pytorch.org/whl/cpu --no-cache-dir
# Install DGL 1.1.3
pip install --no-deps dgl==1.1.3 -f https://data.dgl.ai/wheels/repo.html
# Install e3nn
pip install --no-cache-dir --no-binary e3nn e3nn
# Install custom mdtraj
pip install git+https://github.com/huhlim/mdtraj --no-cache-dir

# Install SE3Transformer (Source build with patches)
echo -e "${BLUE}Installing SE3Transformer...${NC}"
SE3T_SRC="$TEMP_DIR/se3t-src"
git clone https://github.com/huhlim/SE3Transformer "$SE3T_SRC"
pushd "$SE3T_SRC"
# Patch dependencies to match what we have
sed -i 's/python = "[^"]*"/python = ">=3.7"/' pyproject.toml
sed -i 's/torch = "[^"]*"/torch = ">=2.1.0"/' pyproject.toml
pip install . --no-cache-dir
popd

# Install cg2all (Source build with patches)
echo -e "${BLUE}Installing cg2all...${NC}"
CG2ALL_SRC="$TEMP_DIR/cg2all-src"
git clone https://github.com/huhlim/cg2all.git "$CG2ALL_SRC"
pushd "$CG2ALL_SRC"
# Checkout specific commit known to work
git checkout a789cb5
# Patch dependencies
sed -i 's/torch = "[^"]*"/torch = ">=2.1.0"/' pyproject.toml
sed -i 's/numpy = "[^"]1"/numpy = ">=1.21"/' pyproject.toml
pip install . --no-cache-dir
popd

# Get path to this environment
CG2ALL_PYTHON=$(which python)
CG2ALL_ENV_PATH=$(dirname $(dirname $CG2ALL_PYTHON))

# Switch back to main env
conda activate $ENV_NAME

# ---------------------------------------------------------
# Final Configuration
# ---------------------------------------------------------
echo -e "${YELLOW}📝 Writing configuration...${NC}"

# Find where CABS is installed
CABS_PATH=$(python -c "import CABS; print(CABS.__path__[0])")
DATA_DIR="$CABS_PATH/data"
mkdir -p "$DATA_DIR"

# Write paths json
echo "{\"cg2all_env_prefix\": \"$CG2ALL_ENV_PATH\"}" > "$DATA_DIR/cabs_paths.json"

# Test installation
echo -e "${YELLOW}🧪 Testing installation...${NC}"
if CABSflex --help > /dev/null 2>&1; then
    echo -e "${GREEN}✅ CABSflex command working!${NC}"
else
    echo -e "${RED}❌ Warning: CABSflex command not working${NC}"
fi

if CABSdock --help > /dev/null 2>&1; then
    echo -e "${GREEN}✅ CABSdock command working!${NC}"
else
    echo -e "${RED}❌ Warning: CABSdock command not working${NC}"
fi

# Cleanup
cd /
rm -rf $TEMP_DIR

echo ""
echo -e "${GREEN}🎉 CABS-flex Beta installation complete!${NC}"
echo "================================="
echo ""
echo -e "${BLUE}To use CABS-flex:${NC}"
echo "  conda activate $ENV_NAME"
echo "  CABSflex --help"
echo "  CABSdock --help"
echo ""
echo -e "${BLUE}For help and bug reports:${NC}"
echo "  📧 Email: k.wroblewski7@uw.edu.pl"
echo "  🐛 Issues: $REPO_URL/issues"
echo ""
echo -e "${YELLOW}Happy testing! 🧬${NC}"
