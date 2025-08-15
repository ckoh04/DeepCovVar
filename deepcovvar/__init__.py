"""
DeepCovVar: Deep Learning-Based COVID-19 Variant Classification Tool
Author: 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR ANY PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

A comprehensive bioinformatics package for predicting and classifying COVID-19 variants
using state-of-the-art deep learning models. This package provides tools for sequence analysis,
feature extraction, and multi-class prediction of SARS-CoV-2 variants.

Key Features:
    - Five-phase prediction pipeline
    - Advanced deep learning models (Keras and PyTorch)
    - Comprehensive sequence analysis
    - Efficient data preprocessing
    - Detailed output reporting
    - Optimized models for space optimization

Modules:
    - covid_classifier: Main COVID-19 classifier class
        - COVIDClassifier: Main classifier with multi-phase prediction
    - utils: Utility modules
        - features: Feature extraction and processing
        - neminer_utils: Utility functions for sequence processing
        - feature_data: Feature data handling and management
        - sequence_converter: Sequence type detection and conversion

Usage Example:
    >>> from deepcovvar.covid_classifier import COVIDClassifier
>>> from deepcovvar.utils import FEATURE
>>> from deepcovvar.utils import preprocess
    >>> 
    >>> # Preprocess sequences
    >>> data = preprocess('input.fasta', 'DPC', [0, 0], 400)
    >>> 
    >>> # Predict using COVIDClassifier
    >>> classifier = COVIDClassifier()
    >>> results = classifier.predict_phase(data, phase=1)

Dependencies:
    Required:
        - TensorFlow >= 2.0
        - Keras
        - PyTorch
        - NumPy
        - Pandas
        - Biopython
        - Scikit-learn
        - Transformers
    
    Optional:
        - CUDA (for GPU acceleration)

Author: 
Version: 1.0.0
License: GPL-3.0
"""

__version__ = "1.0.0"
__author__ = ""
__license__ = "GPL-3.0"

from .covid_classifier import COVIDClassifier
from .utils import FEATURE, preprocess, preprocessdf
from .utils.feature_data import *

# Define public API
__all__ = [
    'COVIDClassifier',
    'FEATURE',
    'preprocess',
    'preprocessdf'
]
