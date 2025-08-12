# DeepCovVar: Deep Learning-Based COVID-19 Variant Classification Tool

A comprehensive bioinformatics package for predicting and classifying COVID-19 variants using state-of-the-art deep learning models.

## Features

- **Five-Phase Classification Pipeline**: Multi-stage classification from virus detection to variant identification
- **Automatic Sequence Processing**: Intelligent detection and conversion of nucleotide sequences to protein sequences
- **Prodigal Integration**: Seamless nucleotide-to-protein conversion using Prodigal gene prediction software
- **Advanced Deep Learning Models**: Keras and PyTorch-based models for optimal performance
- **Quantized Models**: Space-optimized models for efficient deployment
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
pip install -e .
```

## Quick Start

### Command Line Usage

```bash
# Basic usage with automatic nucleotide conversion
python -m DeepCovVar input.fasta --phase 1

# Run complete pipeline (all phases)
python -m DeepCovVar input.fasta --run-all

# Run complete pipeline with custom output directory
python -m DeepCovVar input.fasta --run-all --output-dir results/

# Specify Prodigal path for nucleotide conversion
python -m DeepCovVar input.fasta --prodigal-path /usr/local/bin/prodigal

# Disable automatic conversion
python -m DeepCovVar input.fasta --no-auto-convert

# Force conversion even if sequences appear to be protein
python -m DeepCovVar input.fasta --force-convert
```

### Python API Usage

```python
from DeepCovVar.covid_classifier import COVIDClassifier

# Initialize classifier with automatic nucleotide conversion
classifier = COVIDClassifier(prodigal_path="/usr/local/bin/prodigal")

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
5. **Phase 5**: SARS-CoV-2 Variant Classification

## Model Architecture

- **Phase 1-4**: Keras-based neural networks with quantized weights
- **Phase 5**: PyTorch Transformer model for sequence classification
- **Feature Extraction**: CKSAAP (Composition of k-spaced Amino Acid Pairs)
- **Optimization**: Quantized models for reduced memory footprint

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
