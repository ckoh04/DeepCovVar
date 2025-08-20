"""
DeepCovVar Utilities Package

This package contains utility modules for the DeepCovVar project:
- features: Feature extraction utilities
- feature_data: Feature data handling
- deepcovvar_utils: Utility functions (cleaned up)
- sequence_converter: Sequence type detection and conversion
"""

from .features import FEATURE
from .feature_data import *
from .sequence_converter import SequenceTypeDetector, SequenceProcessor

__all__ = [
    'FEATURE',
    'SequenceTypeDetector',
    'SequenceProcessor'
]



