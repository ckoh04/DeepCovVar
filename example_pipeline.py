#!/usr/bin/env python3
"""
Example script demonstrating DeepCovVar pipeline functionality.

This script shows how to:
1. Process sequences with automatic nucleotide detection and conversion
2. Run the complete classification pipeline
3. Handle different input scenarios
"""

import sys
from pathlib import Path

# Add the DeepCovVar package to the path
sys.path.insert(0, str(Path(__file__).parent / "DeepCovVar"))

from DeepCovVar.covid_classifier import COVIDClassifier

def main():
    """Demonstrate the pipeline functionality."""
    print("DeepCovVar Pipeline Example")
    print("=" * 50)
    
    # Initialize classifier
    # You can specify a custom Prodigal path if needed
    # classifier = COVIDClassifier(prodigal_path="/usr/local/bin/prodigal")
    classifier = COVIDClassifier()
    
    # Example input file (should exist in your system)
    input_file = "DeepCovVar/test_sequences.fasta"
    
    if not Path(input_file).exists():
        print(f"Example file {input_file} not found.")
        print("Please create a test FASTA file or modify the input_file variable.")
        return
    
    print(f"Input file: {input_file}")
    print("\nRunning complete pipeline...")
    
    try:
        # Run all phases automatically
        results = classifier.run_all_phases(
            input_file=input_file,
            output_dir="example_results",
            base_filename="example_pipeline"
        )
        
        print(f"\nPipeline completed successfully!")
        print(f"Results from {len(results)} phases:")
        
        for phase, result in results.items():
            if result is not None:
                print(f"  Phase {phase}: {len(result)} sequences processed")
            else:
                print(f"  Phase {phase}: Failed")
        
        print(f"\nAll results saved to 'example_results/' directory")
        print("Check the generated files for detailed results.")
        
    except Exception as e:
        print(f"Error running pipeline: {e}")
        print("This might be due to missing model files or Prodigal installation.")

if __name__ == "__main__":
    main()
