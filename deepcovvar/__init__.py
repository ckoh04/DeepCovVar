"""
DeepCovVar: Deep Learning-based COVID-19 Variant Classification Tool

A comprehensive tool for classifying COVID-19 variants using deep learning models
across multiple phases of analysis.

Features:
- Multi-phase classification pipeline
- Support for both Keras and PyTorch models
- Automatic sequence type detection and conversion
- Efficient data preprocessing
- deepcovvar_utils: Utility functions (cleaned up)
- sequence_converter: Sequence type detection and conversion
- Optimized model loading and prediction
- Customizable prediction thresholds
- Comprehensive result formatting

Usage:
    >>> from deepcovvar import COVIDClassifier
    >>> classifier = COVIDClassifier()
    >>> results = classifier.predict(1, 'sequences.fasta')

Author: DeepCovVar Team
License: MIT
Version: 1.0.0
"""

from .covid_classifier import COVIDClassifier
from .utils import FEATURE

__version__ = "1.0.0"
__author__ = "DeepCovVar Team"

__all__ = [
    'COVIDClassifier',
    'FEATURE'
]
