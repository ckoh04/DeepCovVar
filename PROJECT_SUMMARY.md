# DeepCovVar Project Summary

## Project Overview

DeepCovVar is a comprehensive deep learning-based tool for COVID-19 variant classification, reorganized from the original `covid-classifier-tool` project. The project has been restructured following modern Python package standards and optimized for space efficiency.

## Project Structure

```
DeepCovVar/
├── DeepCovVar/                    # Main package directory
│   ├── __init__.py               # Package initialization
│   ├── __main__.py               # Command-line interface
│   ├── covid_classifier.py       # Main classifier class
│   ├── features.py               # Feature extraction
│   ├── neminer_utils.py          # Utility functions
│   ├── feature_data.py           # Feature data handling
│   ├── example.py                # Usage examples
│   ├── test_deepcovvar.py        # Test suite
│   ├── models/                   # Quantized model files
│   │   ├── p1_40_final_model_quantized.keras
│   │   ├── p2_final_model_quantized.keras
│   │   ├── p3_final_model_quantized.keras
│   │   ├── p4_final_model_quantized.keras
│   │   └── p2_final_model_7_quantized.keras
│   ├── tests/                    # Test directory
│   └── docs/                     # Documentation directory
├── setup.py                      # Package installation
├── requirements.txt               # Dependencies
├── README.md                     # Project documentation
├── MIGRATION_GUIDE.md            # Migration instructions
├── pytest.ini                    # Test configuration
├── .gitignore                    # Git ignore rules
├── run_deepcovvar.sh             # Convenience shell script
└── test data files               # FASTA test sequences
```

## Key Features

### 1. Five-Phase Classification Pipeline
- **Phase 1**: Virus vs Non-virus classification
- **Phase 2**: (+)ssRNA (Class IV) vs Others
- **Phase 3**: Coronavirus vs Other ssRNA(+)
- **Phase 4**: SARS-CoV-2 vs Other Coronaviruses
- **Phase 5**: SARS-CoV-2 Variant Classification (PyTorch Transformer)

### 2. Model Architecture
- **Phases 1-4**: Keras-based neural networks with quantized weights
- **Phase 5**: PyTorch Transformer model (ESM-2) with quantized weights
- **Feature Extraction**: CKSAAP (Composition of k-spaced Amino Acid Pairs)
- **Optimization**: Quantized models for reduced memory footprint

### 3. Usage Modes
- **Command-Line Interface**: `python -m DeepCovVar -f input.fasta -o results -p 1`
- **Shell Script**: `./run_deepcovvar.sh input.fasta results 1`
- **Python API**: Direct import and usage in scripts

## Space Optimization

The project now uses quantized models exclusively:
- **Original models**: ~1.2GB each
- **Quantized models**: ~300-600MB each
- **Total space saved**: ~50-75% reduction
- **Models included**: All 5 phases with quantized versions

## Installation

```bash
# Clone and install
git clone <repository>
cd DeepCovVar
pip install -r requirements.txt
pip install -e .

# Test installation
python -m DeepCovVar.test_deepcovvar
```

## Usage Examples

### Basic Classification
```bash
python -m DeepCovVar -f input.fasta -o results -p 1
```

### With Verbose Logging
```bash
python -m DeepCovVar -f sequences.fa -o output -p 4 --verbose
```

### Python API
```python
from DeepCovVar.covid_classifier import COVIDClassifier

classifier = COVIDClassifier()
classifier.load_model(phase=1)
prediction = classifier.predict_phase(sequence, phase=1)
```

## Testing

The project includes a comprehensive test suite:
```bash
# Run all tests
python -m DeepCovVar.test_deepcovvar

# Run with pytest
pytest DeepCovVar/tests/
```

## Dependencies

### Core Requirements
- Python >= 3.8
- TensorFlow >= 2.13.0
- PyTorch >= 2.0.0
- Keras >= 2.13.1
- Biopython >= 1.81
- NumPy >= 1.24.0
- Pandas >= 2.1.0
- Transformers >= 4.20.0

### Optional
- GPU support (CUDA)
- Development tools (pytest, black, flake8)
- Documentation tools (Sphinx)

## Migration from covid-classifier-tool

The project has been completely reorganized from the original structure:
- **Old name**: SARSuite (internal)
- **New name**: DeepCovVar (consistent)
- **Old structure**: Scattered files in covid-classifier-tool
- **New structure**: Proper Python package with models integrated
- **Model optimization**: Quantized models for space efficiency

See `MIGRATION_GUIDE.md` for detailed migration instructions.

## Benefits of Reorganization

1. **Professional Structure**: Follows Python packaging best practices
2. **Space Efficiency**: Quantized models reduce storage requirements
3. **Easy Installation**: Standard pip installation process
4. **Better Testing**: Comprehensive test coverage
5. **Documentation**: Improved README and examples
6. **Maintainability**: Cleaner code organization
7. **CLI Interface**: Built-in command-line tool
8. **API Access**: Easy integration into other projects

## Future Enhancements

Potential areas for future development:
- Web interface
- REST API
- Docker containerization
- Additional model architectures
- Real-time classification
- Batch processing optimization
- Integration with bioinformatics pipelines

## License

This project is licensed under the GNU General Public License v3.0.

## Support

For issues, questions, or contributions:
1. Check the test suite first
2. Review the examples
3. Consult the migration guide
4. Verify dependencies are installed correctly
