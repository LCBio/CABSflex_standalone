# 🧬 CABS-flex
![CABS-flex logo](images/2430462593-CABS-flex-logo-1401.jpg)

[![GitHub](https://img.shields.io/badge/repo-GitHub-brightgreen.svg)](https://github.com/LCBio/cabsflex)
[![GitLab Mirror](https://img.shields.io/badge/mirror-GitLab-orange.svg)](https://gitlab.com/lcbio1/CABSflex_standalone.git)
[![Wiki](https://img.shields.io/badge/docs-Wiki-blue.svg)](https://github.com/LCBio/cabsflex/wiki)

**CABS-flex** is a powerful Python package designed for coarse-grained protein modeling and simulation. It leverages the efficient CABS force field to provide high-resolution insights into protein flexibility and structural dynamics.

---

## 🔬 The CABS Modeling Scheme

CABS-flex combines the efficient **CABS coarse-grained model** with structural clustering and all-atom reconstruction:
- **Efficiency**: A single amino acid is represented by 4 pseudo-atoms (CA, CB, Side-chain, and peptide bond center), significantly speeding up simulations.
- **Accuracy**: Default distance restraints and simulation settings are optimized to provide the best possible convergence with consensus MD fluctuations.
- **Pipeline**: Trajectories are clustered into representative models and subsequently reconstructed into all-atom representations using state-of-the-art tools like `cg2all` or `MODELLER`.

For a detailed explanation of the method, see the **[Modeling Scheme Wiki Page](https://github.com/LCBio/cabsflex/wiki/Modeling-Scheme)**.

---

## 🚀 Installation

CABS-flex uses a **multi-environment system** to ensure maximum stability and isolation for its dependencies, particularly the all-atom reconstruction tools.

### 1. Standard Installation (Conda)
Recommended for local machines (Linux, macOS, WSL 2).

```bash
# Clone the repository
git clone https://github.com/LCBio/cabsflex.git
cd cabsflex

# Run the installer
bash install.sh
```

### 2. HPC / Cluster Installation (Venv)
Optimized for High-Performance Computing clusters (e.g., Helios) using environment modules.

```bash
# Clone the repository
git clone https://github.com/LCBio/cabsflex.git
cd cabsflex

# Edit the configuration in install-hpc.sh (e.g., BASE_INSTALL_DIR, TEMP_ROOT)
# Run the HPC installer
bash install-hpc.sh
```
> [!TIP]
> Always run the installation from a worker node on HPC clusters to avoid login node restrictions.

### 3. Usage
Once installed, activate the environment as instructed by the installer:
- **Conda**: `conda activate cabs`
- **HPC**: `source path/to/your/venv/bin/activate`

---

## 🧪 Quick Start

### Protein Flexibility (CABS-flex)
Simulate protein structure flexibility for PDB `1CE1`:
```bash
CABSflex -i 1CE1 -a 10 -y 20 -w output_dir
```

### Protein-Peptide Docking (CABS-dock)
Flexible docking of a peptide sequence to a receptor structure:
```bash
CABSdock -i protein.pdb -p PEPTIDESEQUENCE -a 10 -y 20 -w dock_output
```

---

## 🔗 Resources & Documentation

- 📚 **[Full Documentation (Wiki)](https://github.com/LCBio/cabsflex/wiki)**
- 📖 **[Examples & Tutorials](https://github.com/LCBio/cabsflex/wiki/Examples)**
- ⚙️ **[Options Reference](https://github.com/LCBio/cabsflex/wiki/Options-Reference)**
- 🐛 **[Issue Tracker](https://github.com/LCBio/cabsflex/issues)**

---

## 👨‍💻 Development

```bash
# Run tests
pytest tests/

# Code quality
ruff check CABS tests
black --check CABS tests
mypy CABS --ignore-missing-imports
```

**Contact**: k.wroblewski7@uw.edu.pl  
**Laboratory**: Laboratory of Computational Biology, University of Warsaw