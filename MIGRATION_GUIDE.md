# Migration Guide: From covid-classifier-tool to DeepCovVar

This guide helps you migrate from the old `covid-classifier-tool` project to the new `DeepCovVar` package.

## What Changed

### Project Structure
- **Old**: `covid-classifier-tool/` with scattered files
- **New**: `DeepCovVar/` with proper Python package structure

### Package Name
- **Old**: `SARSuite` (internal name)
- **New**: `DeepCovVar` (consistent naming)

### Model Storage
- **Old**: Models stored in separate `final_models/` directory
- **New**: Models integrated into `DeepCovVar/models/` package directory
- **Optimization**: Now uses quantized models by default for space efficiency

## File Mapping

| Old Location | New Location | Notes |
|--------------|--------------|-------|
| `covid_classifier.py` | `DeepCovVar/covid_classifier.py` | Updated with relative imports |
| `features.py` | `DeepCovVar/features.py` | Updated with relative imports |
| `neminer_utils.py` | `DeepCovVar/neminer_utils.py` | Updated with relative imports |
| `feature_data.py` | `DeepCovVar/feature_data.py` | Updated with relative imports |
| `final_models/*.keras` | `DeepCovVar/models/*quantized.keras` | Quantized models only |
| `test*.fasta` | `DeepCovVar/test*.fasta` | Test data preserved |
| `p*.fasta` | `DeepCovVar/p*.fasta` | Phase-specific test data |

## Usage Changes

### Old Way (covid-classifier-tool)
```bash
cd covid-classifier-tool
python covid_classifier.py --input input.fasta --phase 1
```

### New Way (DeepCovVar)
```bash
# Option 1: Python module
python -m DeepCovVar -f input.fasta -o results -p 1

# Option 2: Shell script
./run_deepcovvar.sh input.fasta results 1

# Option 3: Python API
from DeepCovVar.covid_classifier import COVIDClassifier
classifier = COVIDClassifier()
prediction = classifier.predict_phase(sequence, phase=1)
```

## Installation

### Old Way
```bash
cd covid-classifier-tool
pip install -r requirements.txt
```

### New Way
```bash
cd DeepCovVar
pip install -r requirements.txt
pip install -e .
```

## Benefits of Migration

1. **Better Organization**: Proper Python package structure
2. **Space Optimization**: Quantized models reduce storage requirements
3. **Easier Installation**: Standard Python package installation
4. **Command-Line Interface**: Built-in CLI with proper argument parsing
5. **Better Testing**: Comprehensive test suite
6. **Documentation**: Improved README and examples
7. **Maintainability**: Cleaner code structure with relative imports

## Backward Compatibility

The core functionality remains the same:
- Same classification phases (1-5)
- Same model architecture
- Same feature extraction methods
- Same output format

## Troubleshooting

### Import Errors
If you encounter import errors, ensure you're running from the correct directory:
```bash
cd DeepCovVar
python -m DeepCovVar.test_deepcovvar
```

### Model Loading Issues
Verify quantized models are present:
```bash
ls DeepCovVar/models/*quantized.keras
```

### Dependencies
Install all required dependencies:
```bash
pip install -r requirements.txt
```

## Support

For issues or questions about the migration:
1. Check the test suite: `python -m DeepCovVar.test_deepcovvar`
2. Review the examples: `python DeepCovVar/example.py`
3. Check the README.md for usage instructions
4. Verify your Python environment has all required packages
