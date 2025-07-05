# CABS-flex

Welcome to **CABS-flex** - a Python package for coarse-grained protein modeling and simulation!

CABS-flex provides tools for protein flexibility analysis, docking simulations, and structural predictions using the CABS force field.

> **⚠️ Private Beta**: This is currently a private repository for selected beta testers. The package will be made open source and available on conda-forge in the future.

## 🚀 Installation

### Prerequisites
- Conda or Miniconda installed
- Python 3.8-3.11

### Setup

```bash
# Create environment with all dependencies
conda env create -f environment.yml
conda activate cabs

# Install CABSflex in development mode
pip install -e .

# Test installation
CABSflex --help
CABSdock --help
```

## 🧪 Quick Start

```bash
# CABSflex - protein flexibility simulation
CABSflex -i 1CE1 -a 10 -y 20 -w output_dir

# CABSdock - peptide-protein docking
CABSdock -i protein.pdb -p PEPTIDESEQUENCE -a 10 -y 20 -w dock_output

# For more examples see:
CABSflex --help
CABSdock --help
```

## 🐛 Beta Testing

As a beta tester, please:
- Report any installation issues
- Test on your specific use cases
- Provide feedback on the API and usability
- Open GitHub issues for bugs or feature requests

See [`BETA_TESTING.md`](BETA_TESTING.md) for detailed testing guidelines.

## 🔮 Future Public Release

Once ready for public release:
```bash
# Future installation (will be available)
conda install -c conda-forge -c bioconda cabsflex
```

## 👨‍💻 For Developers

```bash
# Run tests
pytest tests/

# Code quality
ruff check CABS tests
black --check CABS tests
mypy CABS --ignore-missing-imports
```

**Contact**: k.wroblewski7@uw.edu.pl