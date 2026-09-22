#!/bin/bash
# CABS-flex Beta Installer
# Usage: curl -sSL https://raw.githubusercontent.com/LCBio/CABSflex_standalone/main/install-beta.sh | bash -s YOUR_BETA_TOKEN

set -e

# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------
# Power user configuration
BETA_TOKEN=$1
REPO_URL="https://github.com/LCBio/CABSflex_standalone"
ENV_NAME="cabs"
INSTALL_SRC=$(pwd)
IS_LOCAL=false
INSTALL_MODELLER="TRUE"
MODELLER_KEY=""   # <<< SET YOUR MODELLER LICENSE KEY HERE

# e3nn is built from source (--no-binary) by default for compatibility with
# older Linux kernels/glibc. If that build fails on your system, re-run with
# E3NN_USE_WHEEL=TRUE to install the prebuilt wheel instead. See the wiki
# Installation page, Section 6 (Troubleshooting & FAQ).
E3NN_USE_WHEEL="${E3NN_USE_WHEEL:-FALSE}"

# Set TRUE to reuse an existing cg2all reconstruction environment instead of
# rebuilding it from scratch. Used to resume installation after manually
# placing a cg2all checkpoint file (see Section 6 of the wiki Installation
# page) without repeating the multi-minute torch/dgl/e3nn/cg2all build.
RESUME_AFTER_CHECKPOINT="${RESUME_AFTER_CHECKPOINT:-FALSE}"

# Root directory for Micromamba environments. Change this if you want Micromamba installed elsewhere.
export MAMBA_ROOT_PREFIX="$HOME/micromamba"


# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Cross-platform sed wrapper
sedi() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' "$@"
    else
        sed -i "$@"
    fi
}

echo -e "${BLUE}🧬 CABS-flex Beta Installer${NC}"
echo "================================="

if [ -z "$BETA_TOKEN" ]; then
    echo -e "${YELLOW}⚠️  No beta token provided. Assuming local installation from current directory.${NC}"
    IS_LOCAL=true

    # Check if we are in a valid repo
    if [ ! -f "$INSTALL_SRC/pyproject.toml" ] && [ ! -f "$INSTALL_SRC/setup.py" ]; then
        echo -e "${RED}❌ Error: No pyproject.toml or setup.py found in $INSTALL_SRC${NC}"
        echo "For local installation, please run this script from the root of the CABSflex repository."
        exit 1
    fi
    echo -e "${GREEN}✅ Local source found: $INSTALL_SRC${NC}"
else
    echo -e "${GREEN}✅ Beta token provided${NC}"
fi

echo -e "${YELLOW}📋 Checking prerequisites...${NC}"

# Check if git is available
if ! command -v git &> /dev/null; then
    echo -e "${RED}❌ Error: git not found${NC}"
    echo "Please install git first."
    exit 1
fi

# Check for other critical tools
for needed_tool in curl tar; do
    if ! command -v "$needed_tool" &> /dev/null; then
        echo -e "${RED}❌ Error: $needed_tool not found${NC}"
        echo "Please install $needed_tool first."
        exit 1
    fi
done

# ---------------------------------------------------------
# MICROMAMBA BOOTSTRAP (No Bzip2 / No Tar / No Conda needed)
# ---------------------------------------------------------
TEMP_DIR=$(mktemp -d)
BIN_DIR="$HOME/.local/bin"
mkdir -p "$BIN_DIR"
export PATH="$BIN_DIR:$PATH"

setup_micromamba() {
    if ! command -v micromamba &> /dev/null; then
        echo -e "${YELLOW}🚀 Downloading standalone micromamba binary...${NC}"

        # Detect OS and Arch
        OS_TYPE=$(uname -s | tr '[:upper:]' '[:lower:]')
        ARCH_TYPE=$(uname -m)

        # Map names to GitHub release binaries
        # GitHub names: micromamba-linux-64, micromamba-osx-arm64, etc.
        PLATFORM=""
        if [ "$OS_TYPE" = "linux" ]; then
            if [ "$ARCH_TYPE" = "x86_64" ]; then PLATFORM="linux-64";
            elif [ "$ARCH_TYPE" = "aarch64" ]; then PLATFORM="linux-aarch64"; fi
        elif [ "$OS_TYPE" = "darwin" ]; then
            if [ "$ARCH_TYPE" = "x86_64" ]; then PLATFORM="osx-64";
            elif [ "$ARCH_TYPE" = "arm64" ]; then PLATFORM="osx-arm64"; fi
        fi

        if [ -z "$PLATFORM" ]; then
            echo -e "${RED}❌ Unsupported architecture: $OS_TYPE $ARCH_TYPE${NC}"
            exit 1
        fi

        # Download the RAW binary (skips bzip2/tar dependency)
        MAMBA_URL="https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-${PLATFORM}"

        if ! curl -Ls "$MAMBA_URL" -o "$BIN_DIR/micromamba"; then
            echo -e "${RED}❌ Download failed.${NC}"
            exit 1
        fi

        chmod +x "$BIN_DIR/micromamba"
        export MAMBA_EXE="$BIN_DIR/micromamba"
        export MAMBA_ROOT_PREFIX="$MAMBA_ROOT_PREFIX"

        # Initialize shell
        eval "$($MAMBA_EXE shell hook -s bash)"
        echo -e "${GREEN}✅ Micromamba bootstrapped${NC}"
    else
        export MAMBA_EXE=$(which micromamba)
        echo -e "${GREEN}✅ Using existing micromamba${NC}"
    fi

    # Initialize shell (required for activate to work within script)
    eval "$($MAMBA_EXE shell hook -s bash)"
}

setup_micromamba

echo -e "${GREEN}✅ Micromamba ready${NC}"

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
# Environment Preparation
# ---------------------------------------------------------
if [ "$IS_LOCAL" = true ]; then
    # Local Install Logic
    cp "$INSTALL_SRC/environment.yml" "$TEMP_DIR/environment.yml"
    cd "$TEMP_DIR"
else
    # Remote Install Logic
    cd "$TEMP_DIR"
    # Test token by trying to access the repo
    echo -e "${YELLOW}🔐 Validating beta token...${NC}"
    if ! git ls-remote https://$BETA_TOKEN@github.com/LCBio/CABSflex_standalone.git &> /dev/null; then
        echo -e "${RED}❌ Error: Invalid beta token or no repository access${NC}"
        echo "Please check your token or contact the developers at $REPO_URL/issues"
        exit 1
    fi
    echo -e "${GREEN}✅ Token validated${NC}"

    # Download environment file using GitHub API
    echo -e "${YELLOW}📥 Downloading environment configuration...${NC}"
    if ! curl -H "Authorization: token $BETA_TOKEN" \
         -H "Accept: application/vnd.github.v3.raw" \
         -L "https://api.github.com/repos/LCBio/CABSflex_standalone/contents/environment.yml" \
         -o environment.yml --fail; then
        echo -e "${RED}❌ Error: Failed to download environment file${NC}"
        exit 1
    fi
fi

# Filter environment.yml to remove things we will install separately
# We remove cg2all, dgl, fair-esm (likely cg2all dep), and ensure clean install
sedi '/dgl/d' environment.yml
sedi '/cg2all/d' environment.yml
sedi '/fair-esm/d' environment.yml

# ---------------------------------------------------------
# Main Environment Installation (cabs)
# ---------------------------------------------------------
# RESUME_AFTER_CHECKPOINT also skips rebuilding this main env, but only if
# it is actually found present and functional (never skipped on the flag
# alone) — used when resuming after a manual checkpoint download: by the
# time that failure happens, the main env is already built, and it is set
# up before the reconstruction one, so redoing it on every retry would be
# pure waste. If the main env is not actually there, it is built normally
# regardless of the flag.
MAIN_ENV_PATH="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/$ENV_NAME"
MAIN_ENV_READY=false
if [ "$RESUME_AFTER_CHECKPOINT" = "TRUE" ] && [ -x "$MAIN_ENV_PATH/bin/CABSflex" ]; then
    echo -e "${BLUE}ℹ️  RESUME_AFTER_CHECKPOINT=TRUE and an existing, functional main environment was found at $MAIN_ENV_PATH — skipping rebuild.${NC}"
    MAIN_ENV_READY=true
fi

echo -e "${YELLOW}🔧 Creating/Updating environment '${ENV_NAME}'...${NC}"

# Isolate Micromamba solve from system Python packages
echo -e "${BLUE}ℹ️  Isolating installation from system Python...${NC}"
export PYTHONPATH=""
export PYTHONHOME=""
export PYTHONUSERBASE=""
unset LD_LIBRARY_PATH

if [ "$MAIN_ENV_READY" = false ]; then
# Use micromamba with --override-channels to bypass Anaconda ToS issues
if $MAMBA_EXE env list | grep -q "^${ENV_NAME} "; then
    echo -e "${BLUE}ℹ️  Updating existing environment...${NC}"
    "$MAMBA_EXE" install -v -y -n $ENV_NAME -c conda-forge -c bioconda -c salilab --override-channels \
        python=3.10 pip modeller dssp gfortran binutils openmm
else
    echo -e "${BLUE}ℹ️  Creating new environment...${NC}"
    "$MAMBA_EXE" create -v -y -n $ENV_NAME -c conda-forge -c bioconda -c salilab --override-channels \
        python=3.10 pip modeller dssp gfortran binutils openmm
fi
fi

# ---------------------------------------------------------
# Modeller Source Fallback Function
# ---------------------------------------------------------
install_modeller_source() {
    echo -e "${YELLOW}📥 Installing Modeller from source (fallback)...${NC}"
    local mod_version="10.7"
    local arch_index="2"
    [ "$(uname -m)" = "arm64" ] || [ "$(uname -m)" = "aarch64" ] && arch_index="10"

    local env_path=$($MAMBA_EXE info --envs | grep "^${ENV_NAME} " | awk '{print $NF}')
    local install_dir="$env_path/modeller"

    cd "$TEMP_DIR"
    curl -L "https://salilab.org/modeller/${mod_version}/modeller-${mod_version}.tar.gz" -o "modeller.tar.gz" --fail
    tar -xzf "modeller.tar.gz"
    cd "modeller-${mod_version}"

    echo -e "${YELLOW}🤖 Running Modeller Installer...${NC}"
    ./Install <<EOF
$arch_index
$install_dir
$MODELLER_KEY


EOF

    echo -e "${YELLOW}🔗 Linking Modeller to Python...${NC}"
    local site_pkgs=$("$env_path/bin/python" -c 'import site; print(site.getsitepackages()[0])')
    echo "$install_dir/modlib" > "$site_pkgs/modeller.pth"
    # Logic for LD_LIBRARY_PATH
    export LD_LIBRARY_PATH="$install_dir/lib/x86_64-intel8:$LD_LIBRARY_PATH"
    echo -e "${GREEN}✅ Modeller source installation complete${NC}"
}

# Verify Modeller installation and fallback if needed
if ! "$MAMBA_EXE" list -n $ENV_NAME | grep -q "modeller"; then
    echo -e "${YELLOW}⚠️  Modeller not found via Mamba. Attempting source installation...${NC}"
    if [ ! -z "$MODELLER_KEY" ]; then
        install_modeller_source
    else
        echo -e "${RED}❌ Error: Modeller key missing. Cannot install from source.${NC}"
         
    fi
fi

# Retry function for pip/git
run_with_retry() {
    local n=1
    local max=5
    local delay=5
    while true; do
        "$@" && break
        if [[ $n -lt $max ]]; then
            ((n++))
            echo -e "${YELLOW}⚠️  Command failed. Attempt $n/$max. Retrying in $delay seconds...${NC}"
            sleep $delay
        else
            echo -e "${RED}❌ Command failed after $max attempts.${NC}"
            exit 1
        fi
    done
}

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
            echo -e "${YELLOW}   Then re-run: RESUME_AFTER_CHECKPOINT=TRUE bash install.sh${NC}"
            echo -e "${YELLOW}   See the wiki Installation page, Section 6 (Troubleshooting & FAQ).${NC}"
            return 1
        fi
        ((attempt++))
        sleep $delay
    done
}

micromamba config append channels conda-forge
micromamba config remove channels defaults

# Activate environment
micromamba activate $ENV_NAME
echo -e "${BLUE}ℹ️  Using Python: $(which python)${NC}"

CG2ALL_ENV_NAME="${ENV_NAME}_reconstruct"
CG2ALL_ENV_PATH="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/$CG2ALL_ENV_NAME"

if [ "$IS_LOCAL" = false ]; then
    echo -e "${YELLOW}📥 Cloning repository for configuration injection...${NC}"
    INSTALL_SRC="$TEMP_DIR/cabs_src"
    [ -d "$INSTALL_SRC" ] || git clone https://$BETA_TOKEN@github.com/LCBio/CABSflex_standalone.git "$INSTALL_SRC"
fi

DATA_DIR="$INSTALL_SRC/CABS/data"
echo -e "${YELLOW}📝 Preparing CABS configuration in $DATA_DIR...${NC}"
mkdir -p "$DATA_DIR"
CABS_ENV_PATH="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/$ENV_NAME"
echo "{\"cg2all_env_prefix\": \"$CG2ALL_ENV_PATH\", \"cabs_env_prefix\": \"$CABS_ENV_PATH\"}" > "$DATA_DIR/cabs_paths.json"
echo -e "${GREEN}✅ Created cabs_paths.json configuration.${NC}"

echo -e "${YELLOW}🧹 Cleaning up build artifacts from $INSTALL_SRC...${NC}"
rm -rf "$INSTALL_SRC/tests/test_cli_options" "$INSTALL_SRC/build" "$INSTALL_SRC/dist" "$INSTALL_SRC"/*.egg-info

if [ "$MAIN_ENV_READY" = false ]; then
# ---------------------------------------------------------
# Dependency Installation via micromamba (No solver bounds)
# ---------------------------------------------------------

echo -e "${YELLOW}📦 Installing dependencies via micromamba...${NC}"
echo -e "${YELLOW}📦 Installing core libraries...${NC}"
micromamba install -n $ENV_NAME -y -c conda-forge --override-channels \
    python=3.10 numpy matplotlib requests tqdm scipy h5py netcdf4 \
    pandas plotly jupyter nbconvert

# echo -e "${YELLOW}📦 Installing dev tools...${NC}"
# micromamba install -n $ENV_NAME -y -c conda-forge --override-channels \
#     biopandas biopython mdtraj pytest pytest-cov black ruff mypy pre-commit \
#     pytest-mock pytest-benchmark pytest-html bandit

echo -e "${YELLOW}📦 Installing biopandas and biopython...${NC}"
micromamba install -n $ENV_NAME -y -c conda-forge --override-channels \
    biopandas biopython

echo -e "${YELLOW}📦 Installing CI/CD & Linting Tools...${NC}"
micromamba install -n $ENV_NAME -y -c conda-forge --override-channels \
    pytest pytest-cov black ruff mypy pre-commit \
    pytest-mock pytest-benchmark pytest-html bandit

echo -e "${YELLOW}📦 Installing mdtraj...${NC}"
run_with_retry micromamba run -n $ENV_NAME pip install --no-cache-dir mdtraj




# Install CABSflex from source (Local or the Temp Clone)
echo -e "${YELLOW}📦 Installing CABSflex from $INSTALL_SRC...${NC}"
micromamba run -n $ENV_NAME pip install --upgrade "$INSTALL_SRC"
fi

# Modeller Config Logic
if [ ! -z "$MODELLER_KEY" ]; then
    echo -e "${YELLOW}🔧 Configuring Modeller for environment '$ENV_NAME'...${NC}"
    MOD_CONFIG=$("$MAMBA_EXE" run -n "$ENV_NAME" python -c "import modlib.modeller.config as c; print(c.__file__)" 2>/dev/null || echo "")

    if [ -z "$MOD_CONFIG" ]; then
        MOD_CONFIG=$("$MAMBA_EXE" run -n "$ENV_NAME" python -c "import modeller; import os; print(os.path.join(os.path.dirname(modeller.__file__), 'config.py'))" 2>/dev/null || echo "")
    fi

    if [ -z "$MOD_CONFIG" ] || [ ! -f "$MOD_CONFIG" ]; then
        ENV_PATH=$("$MAMBA_EXE" info --envs | awk -v env="$ENV_NAME" '$1 == env {print $NF}')
        [ -z "$ENV_PATH" ] && ENV_PATH="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/$ENV_NAME"
        MOD_CONFIG=$(find "$ENV_PATH" -name "config.py" | grep "/modeller/" | head -n 1 || echo "")
    fi

    if [ ! -z "$MOD_CONFIG" ] && [ -f "$MOD_CONFIG" ]; then
        sedi "s/license = .*/license = r'${MODELLER_KEY}'/" "$MOD_CONFIG"
        echo -e "${GREEN}✅ Modeller license configured in $MOD_CONFIG${NC}"
    else
        echo -e "${RED}⚠️  Could not find Modeller config file in $ENV_NAME environment.${NC}"
    fi
fi

CG2ALL_READY=false
if [ "$RESUME_AFTER_CHECKPOINT" = "TRUE" ] && [ -x "$CG2ALL_ENV_PATH/bin/convert_cg2all" ]; then
    echo -e "${BLUE}ℹ️  RESUME_AFTER_CHECKPOINT=TRUE and an existing reconstruction environment was found at $CG2ALL_ENV_PATH — skipping rebuild.${NC}"
    CG2ALL_READY=true
fi

if [ "$CG2ALL_READY" = false ]; then
# ---------------------------------------------------------
# Setup CG2ALL (Reconstruction) in Separate Environment
# ---------------------------------------------------------
echo -e "${YELLOW}🔧 Setting up Reconstruction Environment (cg2all)...${NC}"

"$MAMBA_EXE" create -n $CG2ALL_ENV_NAME -y -c conda-forge --override-channels python=3.9 pip c-compiler cxx-compiler make

echo -e "${BLUE}Installing PyTorch, TorchVision, DGL, E3NN dependencies...${NC}"

# Dynamic check to force CPU on non-Mac, but allow Mac to resolve natively
if [[ "$OSTYPE" != "darwin"* ]]; then
    run_with_retry micromamba run -n $CG2ALL_ENV_NAME pip install \
        torch==2.2.0+cpu \
        torchvision==0.17.0+cpu \
        --index-url https://download.pytorch.org/whl/cpu \
        --no-cache-dir
else
    run_with_retry micromamba run -n $CG2ALL_ENV_NAME pip install \
        torch==2.2.0 \
        torchvision==0.17.0 \
        --index-url https://download.pytorch.org/whl/cpu \
        --no-cache-dir
fi

# Install psutil and tqdm
run_with_retry micromamba run -n $CG2ALL_ENV_NAME pip install \
    psutil>=5.8.0 \
    tqdm \
    --no-cache-dir

# Install dgl
run_with_retry micromamba run -n $CG2ALL_ENV_NAME pip install \
    --no-deps \
    dgl==1.1.3 \
    -f https://data.dgl.ai/wheels/repo.html

# Install e3nn
install_e3nn() {
    local n=1 max=5 delay=5
    local pip_args=(--no-cache-dir --no-binary e3nn e3nn)
    if [ "$E3NN_USE_WHEEL" = "TRUE" ]; then
        echo -e "${BLUE}ℹ️  E3NN_USE_WHEEL=TRUE: installing e3nn from the prebuilt wheel.${NC}"
        pip_args=(--no-cache-dir e3nn)
    fi
    while true; do
        micromamba run -n $CG2ALL_ENV_NAME pip install "${pip_args[@]}" && return 0
        if [[ $n -lt $max ]]; then
            ((n++))
            echo -e "${YELLOW}⚠️  e3nn install failed. Attempt $n/$max. Retrying in $delay seconds...${NC}"
            sleep $delay
        else
            echo -e "${RED}❌ e3nn install failed after $max attempts.${NC}"
            if [ "$E3NN_USE_WHEEL" != "TRUE" ]; then
                echo -e "${YELLOW}ℹ️  By default this installer builds e3nn from source (--no-binary) for compatibility with older Linux kernels.${NC}"
                echo -e "${YELLOW}ℹ️  If your system doesn't need that, retry with the prebuilt wheel instead:${NC}"
                echo -e "${YELLOW}      E3NN_USE_WHEEL=TRUE bash install.sh${NC}"
            fi
            echo -e "${YELLOW}   See the wiki Installation page, Section 6 (Troubleshooting & FAQ).${NC}"
            exit 1
        fi
    done
}
install_e3nn

# SE3Transformer
echo -e "${BLUE}Installing SE3Transformer...${NC}"
SE3T_SRC="$TEMP_DIR/se3t-src"
run_with_retry git clone https://github.com/huhlim/SE3Transformer "$SE3T_SRC"
pushd "$SE3T_SRC"
sedi 's/python = "[^"]*"/python = ">=3.7"/' pyproject.toml
sedi 's/torch = "[^"]*"/torch = ">=2.1.0"/' pyproject.toml
micromamba run -n $CG2ALL_ENV_NAME pip install . --no-cache-dir
popd

# Install cg2all
echo -e "${BLUE}Installing cg2all...${NC}"
CG2ALL_SRC="$TEMP_DIR/cg2all-src"
git clone https://github.com/huhlim/cg2all.git "$CG2ALL_SRC"
pushd "$CG2ALL_SRC"
# Checkout specific commit known to work
git checkout a789cb5
# Patch dependencies
sedi 's/torch = "[^"]*"/torch = ">=2.1.0"/' pyproject.toml
sedi 's/numpy = "[^"]1"/numpy = ">=1.21"/' pyproject.toml
micromamba run -n $CG2ALL_ENV_NAME pip install . --no-cache-dir
popd
fi

# ---------------------------------------------------------
# Verify cg2all checkpoint files (see CABS/utils/utils.py
# CG2ALL_REPRESENTATIONS for the models CABS-flex actually uses)
# ---------------------------------------------------------
echo -e "${YELLOW}🔍 Verifying cg2all checkpoint files...${NC}"
MODEL_HOME=$("$CG2ALL_ENV_PATH/bin/python" -c "import cg2all.lib.libconfig as c; print(c.MODEL_HOME)")
CKPT_OK=true
fetch_checkpoint "CalphaBasedModel.ckpt" "0b51f6fe4a12c878ec28b194e55099d3" "$MODEL_HOME" || CKPT_OK=false
fetch_checkpoint "CalphaSCModel.ckpt"    "d42f297f94b4ea33dafa4145b6495344" "$MODEL_HOME" || CKPT_OK=false
if [ "$CKPT_OK" = false ]; then
    echo -e "${RED}❌ Installation stopped: cg2all checkpoint verification failed (see instructions above).${NC}"
    exit 1
fi

micromamba activate $ENV_NAME

# ---------------------------------------------------------
# Detect Shell and Configure Profile
# ---------------------------------------------------------
DETECTED_SHELL=$(basename "$SHELL")
echo -e "${YELLOW}📝 Configuring shell: $DETECTED_SHELL...${NC}"

case "$DETECTED_SHELL" in
    zsh)  PROFILE_FILE="$HOME/.zshrc" ;;
    bash) PROFILE_FILE="$HOME/.bashrc"; [[ "$OSTYPE" == "darwin"* ]] && PROFILE_FILE="$HOME/.bash_profile" ;;
    csh|tcsh) PROFILE_FILE="$HOME/.cshrc" ;;
    *)    PROFILE_FILE="$HOME/.profile" ;;
esac

if [[ "$DETECTED_SHELL" == *"csh"* ]]; then
    if ! grep -q "setenv PATH.*$HOME/.local/bin" "$PROFILE_FILE"; then
        echo 'setenv PATH "$HOME/.local/bin:$PATH"' >> "$PROFILE_FILE"
        echo -e "${GREEN}✅ Added ~/.local/bin to $PROFILE_FILE${NC}"
    fi
else
    if ! grep -q "export PATH.*$HOME/.local/bin" "$PROFILE_FILE"; then
        echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$PROFILE_FILE"
        echo -e "${GREEN}✅ Added ~/.local/bin to $PROFILE_FILE${NC}"
    fi
fi

"$MAMBA_EXE" shell init -s "$DETECTED_SHELL" --root-prefix="$HOME/micromamba"
echo -e "${GREEN}✅ Shell configured in $PROFILE_FILE${NC}"

# ---------------------------------------------------------
# Shell-Aware Testing & Final Instructions
# ---------------------------------------------------------
echo -e "${YELLOW}🧪 Testing installation for $DETECTED_SHELL...${NC}"

MAIN_PREFIX=$($MAMBA_EXE info --envs | awk -v env="$ENV_NAME" '$1 == env {print $NF}')
MAIN_BIN_DIR="$MAIN_PREFIX/bin"
RECON_PREFIX=$($MAMBA_EXE info --envs | awk -v env="$CG2ALL_ENV_NAME" '$1 == env {print $NF}')
RECON_BIN_DIR="$RECON_PREFIX/bin"

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

echo -e "${BLUE}Checking Main Environment binaries:${NC}"
test_binary "$MAIN_BIN_DIR/CABSflex"
test_binary "$MAIN_BIN_DIR/CABSdock"

echo -e "${BLUE}Checking Reconstruction Environment binaries:${NC}"
test_binary "$RECON_BIN_DIR/convert_cg2all"

rm -rf $TEMP_DIR

echo ""
echo -e "${GREEN}🎉 CABS-flex Beta installation complete!${NC}"
echo "================================="
echo ""
echo -e "${BLUE}To use CABS-flex:${NC}"
echo "  source $PROFILE_FILE"
echo "  micromamba activate $ENV_NAME"
echo "  CABSflex --help"
echo "  CABSdock --help"
echo ""
echo -e "${BLUE}For help and bug reports:${NC}"
echo "  🐛 Issues: $REPO_URL/issues"
echo ""
