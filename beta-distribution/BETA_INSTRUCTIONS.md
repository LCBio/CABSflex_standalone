# 🧬 CABS-flex Beta Testing Instructions

Thank you for participating in the CABS-flex beta program!

## Quick Installation

1. **Download the beta archive**: `cabsflex-beta-YYYYMMDD.tar.gz`
2. **Extract it**:
   ```bash
   tar -xzf cabsflex-beta-YYYYMMDD.tar.gz
   cd cabsflex-beta-YYYYMMDD
   ```
3. **Install**:
   ```bash
   conda env create -f environment.yml
   conda activate cabs
   pip install -e .
   ```

## Test It Works

```bash
conda activate cabs
CABSflex --help
CABSdock --help
```

## Quick Test Run

```bash
# Test CABSflex with a simple protein (short simulation)
CABSflex -i 1CE1 -a 5 -y 10 -s 10 -w test_output

# Test CABSdock with a simple peptide docking (very short)
CABSdock -i 1CE1 -p ACDEFGHIKLMNPQRSTVWY -a 3 -y 5 -s 5 -w dock_test
```

## Report Issues

- 📧 **Email**: k.wroblewski7@uw.edu.pl jan.kulczycki1@gmail.com
- 🐛 **Subject**: [CABS-flex Beta] Your issue description
- 📝 **Include**: 
  - What you were trying to do
  - Error messages (if any)
  - Your operating system
  - Python version: `python --version`

## Documentation

- See `README.md` for detailed usage
- See `BETA_TESTING.md` for beta-specific notes

Thank you for testing! 🚀
