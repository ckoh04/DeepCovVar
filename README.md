# DeepCovVar: Deep Learning-Based COVID-19 Variant Classification Tool

A comprehensive bioinformatics package for predicting and classifying COVID-19 variants using state-of-the-art deep learning models.

## Features

- **Five-Phase Classification Pipeline**: Multi-stage classification from virus detection to variant identification
- **Advanced Deep Learning Models**: Keras and PyTorch-based models for optimal performance
- **Quantized Models**: Space-optimized models for efficient deployment
- **Comprehensive Sequence Analysis**: Feature extraction and preprocessing capabilities
- **Command-Line Interface**: Easy-to-use CLI for batch processing
- **Python API**: Programmatic access for integration into workflows

## Installation

```bash
# Clone the repository
git clone https://github.com/ai-lab/DeepCovVar.git
cd DeepCovVar

# Install dependencies
pip install -r requirements.txt

# Install the package
pip install -e .
```

## Quick Start

### Command Line Usage

```bash
# Basic usage
python -m DeepCovVar -f input.fasta -o results -p 1

# With verbose logging
python -m DeepCovVar -f sequences.fa -o output -p 4 --verbose

# Specify custom model directory
python -m DeepCovVar -f input.fasta -o results -p 2 --model-dir /path/to/models
```

### Python API Usage

```python
from DeepCovVar.covid_classifier import COVIDClassifier

# Initialize classifier
classifier = COVIDClassifier()

# Load model for specific phase
classifier.load_model(phase=1)

# Predict on a sequence
prediction = classifier.predict_phase("ATGCGATCGATCG...", phase=1)
print(f"Prediction: {prediction}")
```

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
@software{deepcovvar2024,
  title={DeepCovVar: Deep Learning-Based COVID-19 Variant Classification Tool},
  author={AI Lab Team},
  year={2024},
  url={https://github.com/ai-lab/DeepCovVar}
}
```
