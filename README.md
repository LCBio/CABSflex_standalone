![CABS-flex logo](images/CABSflex3.png)

[![GitHub](https://img.shields.io/badge/repo-GitHub-brightgreen.svg)](https://github.com/LCBio/CABSflex_standalone)
[![GitLab Mirror](https://img.shields.io/badge/mirror-GitLab-orange.svg)](https://gitlab.com/lcbio1/CABSflex_standalone.git)
[![Wiki](https://img.shields.io/badge/docs-Wiki-blue.svg)](https://github.com/LCBio/CABSflex_standalone/wiki)

**CABS-flex standalone 3** is a unified command-line environment for coarse-grained protein modeling and simulation. It brings together protein flexibility simulations, de novo peptide modeling, and flexible peptide–protein docking into a single Python 3-based framework.

---

## 🔬 Core Workflows

CABS-flex 3 supports three main modeling pipelines:

*   **Protein Flexibility**: Fast simulations of protein conformational dynamics, providing fluctuation profiles that align with MD simulations and NMR ensembles.
*   **Peptide Modeling**: *De novo* structure prediction of linear and cyclic peptides (including backbone-cyclized and disulfide-bonded peptides).
*   **Peptide–Protein Docking**: Flexible docking of peptides to protein receptors, supporting cases with unknown binding sites or significant receptor flexibility.

For a detailed explanation of the methodology, see the **[Modeling Scheme](https://github.com/LCBio/CABSflex_standalone/wiki/Modeling-Scheme)**.

---

## 🚀 Installation

CABS-flex uses a **multi-environment system** to ensure maximum stability and isolation for its dependencies, particularly the all-atom reconstruction tools.

### 1. Standard Installation (Micromamba / Conda)
Recommended for local machines (Linux, macOS, WSL 2). The installer automatically bootstraps **Micromamba** if no environment manager is found.

```bash
# Clone from GitHub
git clone https://github.com/LCBio/CABSflex_standalone.git cabsflex
# OR from GitLab mirror:
# git clone https://gitlab.com/lcbio1/CABSflex_standalone.git cabsflex

# Enter the cloned directory
cd cabsflex

# Run the installer
bash install.sh
```

### 2. HPC / Cluster Installation
Optimized for High-Performance Computing clusters using environment modules or standard virtual environments.

```bash
# Clone from GitHub (specifying folder name 'cabsflex')
git clone https://github.com/LCBio/CABSflex_standalone.git cabsflex
# OR from GitLab mirror (specifying folder name 'cabsflex'):
# git clone https://gitlab.com/lcbio1/CABSflex_standalone.git cabsflex

# Enter the cloned directory
cd cabsflex

# Option A: HPC Micromamba installer (Recommended)
bash install-hpc-micromamba.sh

# Option B: Standard Venv installer (for restricted clusters)
bash install-hpc.sh
```

> [!TIP]
> See the **[Installation Guide](https://github.com/LCBio/CABSflex_standalone/wiki/Installation)** for detailed configuration tips for specific HPC environments.

---

## 🧪 Quick Start

Once installed, activate the environment (`micromamba activate cabs`) and run a short verification simulation:

### Protein Flexibility
```bash
CABSflex -i 1HPW -a 10 -y 20 -w output_dir
```

### Peptide–Protein Docking
```bash
CABSdock -i receptor.pdb -p PEPTIDESEQUENCE -a 10 -y 20 -w dock_output
```

---

## 🔗 Resources & Documentation

*   📚 **[Full Documentation (Wiki)](https://github.com/LCBio/CABSflex_standalone/wiki)**
*   📖 **[Examples & Case Studies](https://github.com/LCBio/CABSflex_standalone/wiki/Examples)**
*   ⚙️ **[Options Reference](https://github.com/LCBio/CABSflex_standalone/wiki/Options-Reference)**
*   🐛 **[Issue Tracker](https://github.com/LCBio/CABSflex_standalone/issues)**

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

**Laboratory**: [Laboratory of Computational Biology](https://lcbio.pl/), University of Warsaw