.. DeepCovVar documentation master file, created by
   sphinx-quickstart on 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

DeepCovVar: Deep Learning-Based COVID-19 Variant Classification Tool
====================================================================

:Author: DeepCovVar Team
:Date: |today|
:Version: |version|

*************
Introduction
*************

DeepCovVar is a comprehensive deep learning-based tool for COVID-19 variant classification, 
reorganized from the original covid-classifier-tool project. The project has been restructured 
following modern Python packaging best practices and provides a robust, scalable solution for 
SARS-CoV-2 variant analysis.

**Key Features:**
- Five-phase prediction pipeline for comprehensive variant classification
- Advanced deep learning models (Keras and PyTorch)
- Automatic nucleotide-to-protein sequence conversion
- Optimized models for space efficiency
- Professional Python package structure
- Command-line interface and Python API

**Keywords:**
COVID-19, SARS-CoV-2, variant classification, deep learning, bioinformatics, machine learning, 
neural networks, TensorFlow, PyTorch, Keras

*****************************
Classification Phases
*****************************

DeepCovVar implements a five-phase classification pipeline:

1. **Phase 1**: Virus vs Non-virus classification
2. **Phase 2**: (+)ssRNA (Class IV) vs Others
3. **Phase 3**: Coronavirus vs Other ssRNA(+)
4. **Phase 4**: SARS-CoV-2 vs Other Coronaviruses
5. **Phase 5**: SARS-CoV-2 Variant Classification (Omicron, Alpha, Delta, Epsilon, Iota, Gamma, Others)

Each phase uses specialized deep learning models optimized for accuracy and efficiency.

*****************************
Examples and Case-Studies
*****************************

Example usage and case-studies with real data is provided in the examples section.

Contents
=========

.. toctree::
   :maxdepth: 2
   
   installation.rst
   modules.rst
   usage.rst
   api.rst
   examples.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`



