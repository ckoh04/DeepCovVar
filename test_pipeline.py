#!/usr/bin/env python3
"""
Test script for DeepCovVar pipeline functionality.

This script demonstrates how to use the new sequence processing and pipeline features.
"""

import sys
from pathlib import Path

# Add the DeepCovVar package to the path
sys.path.insert(0, str(Path(__file__).parent / "DeepCovVar"))

from DeepCovVar.covid_classifier import COVIDClassifier
from DeepCovVar.sequence_converter import SequenceProcessor

def test_sequence_detection():
    """Test sequence type detection functionality."""
    print("Testing sequence type detection...")
    
    # Test with protein sequences
    protein_file = "DeepCovVar/test_sequences.fasta"
    
    try:
        processor = SequenceProcessor()
        sequence_types = processor.detector.check_fasta_type(protein_file)
        
        print(f"Sequence types for {protein_file}:")
        for seq_id, is_nucleotide in sequence_types:
            seq_type = "Nucleotide" if is_nucleotide else "Protein"
            print(f"  {seq_id}: {seq_type}")
            
    except Exception as e:
        print(f"Error testing sequence detection: {e}")

def test_pipeline():
    """Test the complete pipeline functionality."""
    print("\nTesting complete pipeline functionality...")
    
    try:
        # Initialize classifier
        classifier = COVIDClassifier()
        
        # Test input file
        input_file = "DeepCovVar/test_sequences.fasta"
        
        if not Path(input_file).exists():
            print(f"Test file {input_file} not found. Skipping pipeline test.")
            return
        
        print(f"Running pipeline on {input_file}...")
        
        # Run all phases
        results = classifier.run_all_phases(
            input_file=input_file,
            output_dir="test_results",
            base_filename="test_pipeline"
        )
        
        print(f"Pipeline completed with {len(results)} phases")
        
    except Exception as e:
        print(f"Error testing pipeline: {e}")

def main():
    """Main test function."""
    print("DeepCovVar Pipeline Test")
    print("=" * 40)
    
    # Test sequence detection
    test_sequence_detection()
    
    # Test pipeline
    test_pipeline()
    
    print("\nTest completed!")

if __name__ == "__main__":
    main()
