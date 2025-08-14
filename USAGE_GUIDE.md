# DeepCovVar Usage Guide

## Problem Fixed

The original error `'COVIDClassifier' object has no attribute 'predict_phase'` occurred because:

1. **Method Name Mismatch**: The main script was calling `classifier.predict_phase()` but the actual method name is `classifier.predict()`
2. **Interface Inconsistency**: The main script expected a method that didn't exist in the COVIDClassifier class

## Solutions Implemented

### 1. Fixed Method Calls
- Changed `classifier.predict_phase()` to `classifier.predict()`
- Updated the main script to use the correct method signature

### 2. Added Complete Pipeline Support
- New `--all-phases` flag to run all phases from 1 to 5 automatically
- Leverages the existing `run_all_phases()` method in COVIDClassifier

### 3. Simplified Interface
- Made output directory optional (defaults to current directory)
- Created a simple wrapper script for minimal arguments

## New Usage Options

### Option 1: Simple Wrapper (Recommended for Beginners)
```bash
# Run complete pipeline with just input file
python run_deepcovvar_simple.py input.fasta

# Run with custom output directory
python run_deepcovvar_simple.py input.fasta results/

# Run with verbose output
python run_deepcovvar_simple.py --verbose input.fasta results/
```

### Option 2: Direct Module Usage
```bash
# Run complete pipeline (recommended)
python -m DeepCovVar -f input.fasta -o output_dir --all-phases

# Run specific phase
python -m DeepCovVar -f input.fasta -o output_dir -p 5

# Use default output directory
python -m DeepCovVar -f input.fasta --all-phases
```

### Option 3: Python API
```python
from DeepCovVar.covid_classifier import COVIDClassifier

classifier = COVIDClassifier()

# Run single phase
results = classifier.predict(phase=1, input_file="sequences.fasta")

# Run complete pipeline
all_results = classifier.run_all_phases(
    input_file="sequences.fasta",
    output_dir="results/"
)
```

## What Each Phase Does

1. **Phase 1**: Virus vs Non-virus classification
2. **Phase 2**: (+)ssRNA (Class IV) vs Others  
3. **Phase 3**: Coronavirus vs Other ssRNA(+)
4. **Phase 4**: SARS-CoV-2 vs Other Coronaviruses
5. **Phase 5**: SARS-CoV-2 Variant Classification (Omicron, Alpha, Delta, etc.)

## Output Files

When running all phases, you'll get:
- `input_phase_1_results.csv` - Phase 1 results
- `input_phase_2_results.csv` - Phase 2 results  
- `input_phase_3_results.csv` - Phase 3 results
- `input_phase_4_results.csv` - Phase 4 results
- `input_phase_5_results.csv` - Phase 5 results
- `input_pipeline_summary.txt` - Summary report

## Troubleshooting

### Common Issues
1. **Model files not found**: Ensure the `models/` directory contains the required model files
2. **Prodigal not found**: Install Prodigal for nucleotide sequence processing
3. **Memory issues**: The tool automatically uses CPU-only mode to avoid GPU memory problems

### Verbose Mode
Use `--verbose` flag to get detailed logging information for debugging.

## Performance Notes

- **Complete Pipeline**: Running all phases takes longer but provides comprehensive analysis
- **Single Phase**: Faster if you only need specific classification
- **Sequence Processing**: Nucleotide sequences are automatically converted to protein sequences using Prodigal
- **Model Loading**: Models are loaded on-demand for each phase

## Best Practices

1. **Use `--all-phases`** for comprehensive analysis
2. **Provide output directory** to organize results
3. **Check model files** exist before running
4. **Use verbose mode** for debugging
5. **Start with simple wrapper** if you're new to the tool
