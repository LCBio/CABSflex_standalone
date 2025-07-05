#!/bin/bash
# CABS-flex Beta Installer
# Usage: curl -sSL https://raw.githubusercontent.com/LCBio/cabsflex/main/install-beta.sh | bash -s YOUR_BETA_TOKEN

set -e

BETA_TOKEN=$1
REPO_URL="https://github.com/LCBio/cabsflex"
ENV_NAME="cabs"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}🧬 CABS-flex Beta Installer${NC}"
echo "================================="

if [ -z "$BETA_TOKEN" ]; then
    echo -e "${RED}❌ Error: Beta token required${NC}"
    echo ""
    echo "Usage:"
    echo "  curl -sSL https://raw.githubusercontent.com/LCBio/cabsflex/main/install-beta.sh | bash -s YOUR_BETA_TOKEN"
    echo ""
    echo "Contact k.wroblewski7@uw.edu.pl to get your beta token."
    exit 1
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

# Create temporary directory
TEMP_DIR=$(mktemp -d)
cd $TEMP_DIR
echo -e "${YELLOW}📂 Working in: $TEMP_DIR${NC}"

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

# Check if environment already exists
if conda env list | grep -q "^${ENV_NAME} "; then
    echo -e "${YELLOW}⚠️  Environment '${ENV_NAME}' already exists${NC}"
    read -p "Do you want to remove and recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo -e "${YELLOW}�️  Removing existing environment...${NC}"
        conda env remove -n $ENV_NAME -y
    else
        echo -e "${BLUE}ℹ️  Updating existing environment...${NC}"
        conda env update -f environment.yml -n $ENV_NAME
        # Still need to update the package
        echo -e "${YELLOW}📦 Updating CABSflex...${NC}"
        eval "$(conda shell.bash hook)"
        conda activate $ENV_NAME
        pip install --upgrade git+https://$BETA_TOKEN@github.com/LCBio/cabsflex.git
        cd /
        rm -rf $TEMP_DIR
        echo -e "${GREEN}✅ CABSflex beta updated successfully!${NC}"
        echo ""
        echo -e "${BLUE}To use CABSflex:${NC}"
        echo "  conda activate $ENV_NAME"
        echo "  CABSflex --help"
        exit 0
    fi
fi

# Create environment
echo -e "${YELLOW}🔧 Creating conda environment '${ENV_NAME}'...${NC}"
echo "This may take a few minutes..."
if ! conda env create -f environment.yml -n $ENV_NAME; then
    echo -e "${RED}❌ Error: Failed to create conda environment${NC}"
    exit 1
fi

# Activate environment and install CABSflex
echo -e "${YELLOW}📦 Installing CABSflex from private repository...${NC}"
eval "$(conda shell.bash hook)"
conda activate $ENV_NAME

if ! pip install git+https://$BETA_TOKEN@github.com/LCBio/cabsflex.git; then
    echo -e "${RED}❌ Error: Failed to install CABSflex${NC}"
    exit 1
fi

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
