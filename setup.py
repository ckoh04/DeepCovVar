#!/usr/bin/env python3
"""
DeepCovVar: Deep Learning-Based COVID-19 Variant Classification Tool
Copyright (C) 2025

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
using state-of-the-art deep learning models.
"""

import os
import setuptools
from pathlib import Path

# Package metadata
NAME = "deepcovvar"
DESCRIPTION = "Deep learning-based COVID-19 variant prediction and classification tool"
AUTHOR = ""
AUTHOR_EMAIL = ""
MAINTAINER = ""
MAINTAINER_EMAIL = ""
# URL = "https://github.com/ckoh04/DeepCovVar"
# DOCUMENTATION = "https://github.com/ckoh04/DeepCovVar/docs"
# REPOSITORY = "https://github.com/ckoh04/DeepCovVar"
LICENSE = "GPL-3.0"
PYTHON_REQUIRES = ">=3.8"

# Version information
VERSION = "1.0.0"  # Following semantic versioning

# Get the long description from README.md
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

def read_requirements(filename="requirements.txt"):
    """Read requirements from file."""
    try:
        with open(filename) as f:
            return [line.strip() for line in f 
                   if line.strip() and not line.startswith("#")]
    except FileNotFoundError:
        return []

# Core package requirements with versions
INSTALL_REQUIRES = [
    'numpy>=1.24.0',
    'pandas',
    'tensorflow>=2.13.0',
    'keras>=2.13.1',
    'torch>=2.0.0',
    'torchvision>=0.15.0',
    'biopython>=1.81',
    'scikit-learn>=1.3.0',
    'scipy',
    'h5py>=3.9.0',
    'matplotlib>=3.7.2',
    'seaborn>=0.12.2',
    'tqdm>=4.66.1',
    'transformers>=4.20.0'
]

# Optional dependencies with specific use cases
EXTRAS_REQUIRE = {
    'gpu': [
        'tensorflow-gpu>=2.13.0',
        'torch>=2.0.0+cu118',
        'cudatoolkit>=11.8.0',
        'cudnn>=8.7.0'
    ],
    'dev': [
        'pytest>=7.4.2',
        'pytest-cov>=4.1.0',
        'black>=23.9.1',
        'flake8>=6.1.0',
        'mypy>=1.5.1',
        'isort>=5.12.0'
    ],
    'docs': [
        'sphinx>=7.1.2',
        'sphinx-rtd-theme>=1.3.0',
        'sphinx-autodoc-typehints>=1.24.0'
    ],
    'viz': [
        'plotly>=5.16.1',
        'dash>=2.13.0'
    ]
}

setuptools.setup(
    # Basic package information
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    maintainer=MAINTAINER,
    maintainer_email=MAINTAINER_EMAIL,
    # url=URL,
    # project_urls={
    #     "Documentation": DOCUMENTATION,
    #     "Source Code": REPOSITORY,
    #     "Bug Tracker": f"{REPOSITORY}/issues",
    # },
    license=LICENSE,
    python_requires=PYTHON_REQUIRES,
    
    # Classifiers for PyPI
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
        "Natural Language :: English",
    ],
    
    # Package discovery
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={
        "deepcovvar": [
            "models/*.keras",
            "models/*.pt",
            "models/*.h5",
            "models/*.json",
            "models/*.txt",
        ],
    },
    
    # Dependencies
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    
    # Entry points for command-line tools
    entry_points={
        "console_scripts": [
            "deepcovvar=deepcovvar.__main__:main",
        ],
    },
    
    # Additional package data
    zip_safe=False,
    
    # Keywords for PyPI search
    keywords=[
        "covid-19",
        "sars-cov-2",
        "variant-classification",
        "deep-learning",
        "bioinformatics",
        "sequence-analysis",
        "machine-learning",
        "neural-networks",
        "tensorflow",
        "pytorch",
        "keras"
    ],
    
    # # Project URLs
    # project_urls={
    #     "Homepage": URL,
    #     "Documentation": DOCUMENTATION,
    #     "Repository": REPOSITORY,
    #     "Bug Reports": f"{REPOSITORY}/issues",
    #     "Source": REPOSITORY,
    # },
)
