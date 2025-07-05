#!/bin/bash
# CABS-flex Beta Release Script
# This script packages a beta release for distribution to testers

set -e

# Configuration
VERSION=$(date +"%Y%m%d")
RELEASE_NAME="cabsflex-beta-${VERSION}"
ARCHIVE_NAME="${RELEASE_NAME}.tar.gz"

echo "Creating CABS-flex beta release: ${RELEASE_NAME}"

# Create temporary directory
rm -rf tmp_release
mkdir tmp_release
cd tmp_release

# Clone the repository (clean copy)
git clone .. ${RELEASE_NAME}
cd ${RELEASE_NAME}

# Remove development files not needed for beta testing
rm -rf .git
rm -rf .github/workflows  # Keep only essential workflows
rm -rf build/
rm -rf dist/
rm -rf *.egg-info/
rm -rf __pycache__/
find . -name "*.pyc" -delete
find . -name "__pycache__" -type d -delete

# Remove sensitive or unnecessary files
rm -f .pre-commit-config.yaml
rm -f codecov.yml
rm -f tox.ini
rm -f Makefile

# Create beta-specific README
cat > README_BETA.md << 'EOF'
# CABS-flex Beta Release

Thank you for participating in the CABS-flex beta testing!

## Quick Start

1. Extract this archive
2. cd cabsflex-beta-YYYYMMDD
3. conda env create -f environment.yml
4. conda activate cabsflex
5. pip install -e .
6. CABSflex --help

## Important Files

- README.md - Main documentation
- BETA_TESTING.md - Beta testing guide
- environment.yml - Conda environment

## Need Help?

- Check BETA_TESTING.md for detailed instructions
- Report issues on GitHub (private repo)
- Contact: k.wroblewski7@uw.edu.pl

Happy testing! 🧬
EOF

# Go back and create archive
cd ..
tar -czf ${ARCHIVE_NAME} ${RELEASE_NAME}

# Move to parent directory
mv ${ARCHIVE_NAME} ..

# Cleanup
cd ..
rm -rf tmp_release

echo "Beta release created: ${ARCHIVE_NAME}"
echo "You can now share this archive with beta testers."
echo ""
echo "Beta testers should:"
echo "1. Extract: tar -xzf ${ARCHIVE_NAME}"
echo "2. Follow instructions in README_BETA.md"
