Installation Guide
==================

This section provides detailed instructions for installing DeepCovVar.

Prerequisites
-------------

- Python 3.8 or higher
- pip package manager
- Git (for cloning from repository)

Installation Methods
-------------------

From Source (Development)
^^^^^^^^^^^^^^^^^^^^^^^^^

1. Clone the repository:
   .. code-block:: bash
   
      git clone https://github.com/ckoh04/DeepCovVar.git
      cd DeepCovVar

2. Install dependencies:
   .. code-block:: bash
   
      pip install -r requirements.txt

3. Install in development mode:
   .. code-block:: bash
   
      pip install -e .

From PyPI (Production)
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   pip install deepcovvar

Dependencies
------------

Core Requirements:
- TensorFlow >= 2.13.0
- PyTorch >= 2.0.0
- Keras >= 2.13.1
- Biopython >= 1.81
- NumPy >= 1.24.0
- Pandas >= 2.1.0
- Transformers >= 4.20.0

Optional Dependencies:
- GPU support (CUDA)
- Development tools (pytest, black, flake8)
- Documentation tools (Sphinx)

Verification
------------

After installation, verify that DeepCovVar is working correctly:

.. code-block:: bash

   python -c "import deepcovvar; print('DeepCovVar installed successfully!')"

Or run the test suite:

.. code-block:: bash

   python -m deepcovvar.test_deepcovvar



