#!/usr/bin/env python3
"""
Example usage of DeepCovVar package
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

def example_basic_usage():
    """Example of basic usage of DeepCovVar."""
    print("DeepCovVar Basic Usage Example")
    print("=" * 40)
    
    try:
        from deepcovvar.covid_classifier import COVIDClassifier
        
        # Initialize classifier
        print("1. Initializing COVIDClassifier...")
        classifier = COVIDClassifier()
        print("   Classifier initialized")
        
        # Load model for phase 1
        print("\n2. Loading model for Phase 1...")
        classifier.load_model(phase=1)
        print("   Phase 1 model loaded")
        
        # Example sequence (SARS-CoV-2 spike protein fragment)
        example_sequence = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF"
        
        print(f"\n3. Predicting on example sequence (length: {len(example_sequence)})...")
        prediction = classifier.predict_phase(example_sequence, phase=1)
        print(f"   Prediction: {prediction}")
        
        print("\n4. Example completed successfully!")
        
    except Exception as e:
        print(f"   Error: {e}")
        return False
    
    return True

def example_batch_processing():
    """Example of batch processing with DeepCovVar."""
    print("\nDeepCovVar Batch Processing Example")
    print("=" * 40)
    
    try:
        from deepcovvar.covid_classifier import COVIDClassifier
        from Bio import SeqIO
        
        # Initialize classifier
        print("1. Initializing COVIDClassifier...")
        classifier = COVIDClassifier()
        print("   Classifier initialized")
        
        # Load model for phase 1
        print("\n2. Loading model for Phase 1...")
        classifier.load_model(phase=1)
        print("   Phase 1 model loaded")
        
        # Check if test file exists
        test_file = Path(__file__).parent / "test.fasta"
        if not test_file.exists():
            print("   Warning: Test file not found, skipping batch example")
            return True
        
        print(f"\n3. Processing sequences from {test_file.name}...")
        
        # Read sequences
        sequences = list(SeqIO.parse(test_file, "fasta"))
        print(f"   Loaded {len(sequences)} sequences")
        
        # Process first few sequences
        max_sequences = min(3, len(sequences))
        print(f"   Processing first {max_sequences} sequences...")
        
        for i, seq_record in enumerate(sequences[:max_sequences]):
            try:
                prediction = classifier.predict_phase(str(seq_record.seq), phase=1)
                print(f"     Sequence {i+1}: {seq_record.id} → {prediction}")
            except Exception as e:
                print(f"     Sequence {i+1}: {seq_record.id} → ERROR: {e}")
        
        print("\n4. Batch processing example completed!")
        
    except Exception as e:
        print(f"   Error: {e}")
        return False
    
    return True

def main():
    """Run all examples."""
    print("DeepCovVar Examples")
    print("=" * 50)
    
    examples = [
        ("Basic Usage", example_basic_usage),
        ("Batch Processing", example_batch_processing),
    ]
    
    for example_name, example_func in examples:
        print(f"\n{example_name}:")
        print("-" * len(example_name))
        
        if not example_func():
            print(f"{example_name} failed")
            return 1
    
    print("\nAll examples completed successfully!")
    print("\nTo run DeepCovVar from command line:")
    print("  python -m deepcovvar -f input.fasta -o results -p 1")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
