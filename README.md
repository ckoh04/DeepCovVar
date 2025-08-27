# DeepCovVar: Deep Learning-Based COVID-19 Variant Classification Tool

A comprehensive bioinformatics package for predicting and classifying COVID-19 variants using state-of-the-art deep learning models.

## Features

- **Five-Phase Classification Pipeline**: Multi-stage classification from virus detection to variant identification
- **Automatic Sequence Processing**: Intelligent detection and conversion of nucleotide sequences to protein sequences
- **Prodigal Integration**: Seamless nucleotide-to-protein conversion using Prodigal gene prediction software
- **Advanced Deep Learning Models**: Keras and PyTorch-based models for optimal performance
- **Optimized Models**: Space-optimized models for efficient deployment
- **Comprehensive Sequence Analysis**: Feature extraction and preprocessing capabilities
- **Command-Line Interface**: Easy-to-use CLI for batch processing with pipeline automation
- **Python API**: Programmatic access for integration into workflows

## Installation

```bash
# Clone the repository
git clone https://github.com/ai-lab/DeepCovVar.git
cd DeepCovVar

# Install dependencies
pip install -r requirements.txt

# Install Prodigal for nucleotide sequence conversion
# Option 1: Using conda
conda install -c bioconda prodigal

# Option 2: From source
git clone https://github.com/hyattpd/Prodigal.git
cd Prodigal
make install

# Install the package
pip install .
```

## Quick Start

### Command Line Usage

```bash
# Run complete pipeline (all phases) - recommended
python -m deepcovvar -f input.fasta -o output_dir --all-phases

# Run specific phase
python -m deepcovvar -f input.fasta -o output_dir -p 5

# Use default output directory (current directory)
python -m deepcovvar -f input.fasta --all-phases

# Run with verbose output
python -m deepcovvar -f input.fasta -o output_dir --all-phases --verbose
```

### Python API Usage

```python
from deepcovvar.covid_classifier import COVIDClassifier

# Initialize classifier
classifier = COVIDClassifier()

# Run single phase
results = classifier.predict(phase=1, input_file="sequences.fasta")

# Run complete pipeline (all phases)
all_results = classifier.run_all_phases(
    input_file="sequences.fasta",
    output_dir="results/",
    base_filename="my_analysis"
)

# Manual sequence processing
processed_file, was_converted = classifier.process_input_sequences("input.fasta")
```

## Sequence Processing

DeepCovVar automatically detects and processes different sequence types:

- **Automatic Detection**: Intelligently identifies nucleotide vs protein sequences
- **Nucleotide Conversion**: Uses Prodigal for high-quality nucleotide-to-protein conversion
- **Seamless Integration**: Works transparently with existing workflows
- **Multiple Modes**: Supports single genome and metagenome processing

## Classification Phases

1. **Phase 1**: Virus vs Non-virus classification
2. **Phase 2**: (+)ssRNA (Class IV) vs Others
3. **Phase 3**: Coronavirus vs Other ssRNA(+)
4. **Phase 4**: SARS-CoV-2 vs Other Coronaviruses
5. **Phase 5**: SARS-CoV-2 Variant Classification (Omicron, Alpha, Delta, Epsilon, Iota, Gamma, Others)

## Model Architecture

- **Phase 1-4**: Keras-based neural networks with optimized weights
- **Phase 5**: PyTorch Transformer model (ESM-2) for sequence classification
- **Feature Extraction**: CKSAAP (Composition of k-spaced Amino Acid Pairs)
- **Optimization**: Optimized models for reduced memory footprint

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
3. **Memory issues**: The tool automatically uses CPU-only mode and batch processing to avoid memory problems

### Verbose Mode
Use `--verbose` flag to get detailed logging information for debugging.

## Performance Notes

- **Complete Pipeline**: Running all phases takes longer but provides comprehensive analysis
- **Single Phase**: Faster if you only need specific classification
- **Sequence Processing**: Nucleotide sequences are automatically converted to protein sequences using Prodigal
- **Model Loading**: Models are loaded on-demand for each phase
- **Memory Management**: Automatic batch processing prevents out-of-memory errors

## Best Practices

1. **Use `--all-phases`** for comprehensive analysis
2. **Provide output directory** to organize results
3. **Check model files** exist before running
4. **Use verbose mode** for debugging
5. **Start with simple wrapper** if you're new to the tool
6. **For large datasets**: Use smaller batch sizes if memory issues occur

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Citation

If you use DeepCovVar in your research, please cite:

```
@software{deepcovvar2025,
  title={DeepCovVar: Deep Learning-Based COVID-19 Variant Classification Tool},
  author={},
  year={2025},
  url={https://github.com/ckoh04/DeepCovVar}
}
```
