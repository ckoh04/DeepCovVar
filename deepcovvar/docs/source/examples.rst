Examples
========

This section provides practical examples of using DeepCovVar.

Basic Usage Examples
-------------------

Simple Classification
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from deepcovvar.covid_classifier import COVIDClassifier

   # Initialize classifier
   classifier = COVIDClassifier()

   # Run single phase
   results = classifier.predict(phase=1, input_file="sequences.fasta")

   print(f"Phase 1 results: {results}")

Complete Pipeline
^^^^^^^^^^^^^^^^

.. code-block:: python

   from deepcovvar.covid_classifier import COVIDClassifier

   # Initialize classifier
   classifier = COVIDClassifier()

   # Run complete pipeline (all phases)
   all_results = classifier.run_all_phases(
       input_file="sequences.fasta",
       output_dir="results/",
       base_filename="my_analysis"
   )

   print("Complete pipeline finished!")
   print(f"Results saved to: results/")

Sequence Processing
^^^^^^^^^^^^^^^^^

.. code-block:: python

   from deepcovvar.covid_classifier import COVIDClassifier

   classifier = COVIDClassifier()

   # Manual sequence processing
   processed_file, was_converted = classifier.process_input_sequences("input.fasta")
   
   if was_converted:
       print("Nucleotide sequences were converted to protein sequences")
   else:
       print("Sequences were already in protein format")

Command Line Examples
--------------------

Basic Classification
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Run all phases (recommended)
   python -m deepcovvar -f input.fasta -o output_dir --all-phases

   # Run specific phase
   python -m deepcovvar -f input.fasta -o output_dir -p 5

   # Use default output directory
   python -m deepcovvar -f input.fasta --all-phases

Shell Script Usage
^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Make script executable
   chmod +x run_deepcovvar.sh

   # Run with shell script
   ./run_deepcovvar.sh input.fasta output_dir 1

Advanced Examples
----------------

Custom Model Directory
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from deepcovvar.covid_classifier import COVIDClassifier

   # Initialize with custom model directory
   classifier = COVIDClassifier(model_dir="/path/to/custom/models")

   # Use as normal
   results = classifier.predict(phase=1, input_file="sequences.fasta")

Batch Processing
^^^^^^^^^^^^^^^

.. code-block:: python

   from deepcovvar.covid_classifier import COVIDClassifier
   import os

   # Initialize classifier
   classifier = COVIDClassifier()

   # Process multiple files
   input_files = ["file1.fasta", "file2.fasta", "file3.fasta"]
   
   for input_file in input_files:
       if os.path.exists(input_file):
           output_dir = f"results_{os.path.splitext(input_file)[0]}"
           results = classifier.run_all_phases(
               input_file=input_file,
               output_dir=output_dir
           )
           print(f"Processed {input_file} -> {output_dir}")

Error Handling
^^^^^^^^^^^^^

.. code-block:: python

   from deepcovvar.covid_classifier import COVIDClassifier

   try:
       classifier = COVIDClassifier()
       results = classifier.predict(phase=1, input_file="sequences.fasta")
   except FileNotFoundError as e:
       print(f"File not found: {e}")
   except Exception as e:
       print(f"Unexpected error: {e}")



