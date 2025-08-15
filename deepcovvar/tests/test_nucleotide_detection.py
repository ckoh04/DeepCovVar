#!/usr/bin/env python3
"""
Test script for nucleotide sequence detection.

This script tests the sequence type detection functionality without requiring Prodigal.
"""

import sys
import os
from pathlib import Path

# Add the parent directory to the path to import DeepCovVar modules
# Handle both direct execution and module import
if __name__ == "__main__":
    # When running directly, use current working directory's parent
    parent_dir = str(Path.cwd().parent)
else:
    # When imported as module, use file location's parent
    parent_dir = str(Path(__file__).parent.parent.absolute())

sys.path.insert(0, parent_dir)

try:
    from sequence_converter import SequenceTypeDetector
except ImportError as e:
    print(f"Error: Could not import sequence_converter module: {e}")
    print(f"Current working directory: {os.getcwd()}")
    print(f"Parent directory: {parent_dir}")
    print(f"Python path: {sys.path[:3]}")
    sys.exit(1)

def test_nucleotide_detection():
    """Test nucleotide sequence detection functionality."""
    print("Testing Nucleotide Sequence Detection")
    print("=" * 40)
    
    # Test file with nucleotide sequences
    test_file = Path(__file__).parent / "test_nucleotide_sequences.fasta"
    
    if not test_file.exists():
        print(f"Test file {test_file} not found!")
        return False
    
    try:
        # Initialize detector
        detector = SequenceTypeDetector()
        
        print(f"Analyzing sequences in: {test_file}")
        print("-" * 40)
        
        # Check sequence types
        sequence_types = detector.check_fasta_type(str(test_file))
        
        print(f"\nResults:")
        print(f"Total sequences: {len(sequence_types)}")
        
        nucleotide_count = sum(1 for _, is_nt in sequence_types if is_nt)
        protein_count = len(sequence_types) - nucleotide_count
        
        print(f"Nucleotide sequences: {nucleotide_count}")
        print(f"Protein sequences: {protein_count}")
        
        print(f"\nDetailed breakdown:")
        for seq_id, is_nucleotide in sequence_types:
            seq_type = "NUCLEOTIDE" if is_nucleotide else "PROTEIN"
            print(f"  {seq_id}: {seq_type}")
        
        if nucleotide_count > 0:
            print(f"\nNucleotide sequences detected successfully!")
            print("These sequences can be converted to protein using Prodigal.")
        else:
            print(f"\nAll sequences appear to be protein.")
        return True
            
    except Exception as e:
        print(f"Error during testing: {e}")
        return False

def test_individual_sequence_detection():
    """Test detection on individual sequences."""
    print("\n" + "="*40)
    print("Testing Individual Sequence Detection")
    print("="*40)
    
    detector = SequenceTypeDetector()
    
    # Test nucleotide sequences
    nucleotide_seqs = [
        "ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACA",
        "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC",
        "ATUGCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUC"
    ]
    
    # Test protein sequences
    protein_seqs = [
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQ",
        "ARNDCQEGHILKMFPSTWYV",
        "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQ"
    ]
    
    print("Nucleotide sequences:")
    for i, seq in enumerate(nucleotide_seqs):
        is_nt = detector.is_nucleotide(seq)
        status = "NUCLEOTIDE" if is_nt else "WRONGLY IDENTIFIED"
        print(f"  {i+1}. {status}")
    
    print("\nProtein sequences:")
    for i, seq in enumerate(protein_seqs):
        is_nt = detector.is_nucleotide(seq)
        status = "PROTEIN" if not is_nt else "WRONGLY IDENTIFIED"
        print(f"  {i+1}. {status}")
    
    return True

def main():
    """Main test function."""
    print("DeepCovVar Nucleotide Detection Tests")
    print("=" * 50)
    
    # Test file-based detection
    test1_passed = test_nucleotide_detection()
    
    # Test individual sequence detection
    test2_passed = test_individual_sequence_detection()
    
    print(f"\n{'='*50}")
    print("Test Summary:")
    print(f"File-based detection: {'PASSED' if test1_passed else 'FAILED'}")
    print(f"Individual detection: {'PASSED' if test2_passed else 'FAILED'}")
    
    if test1_passed and test2_passed:
        print("\nAll tests passed! Nucleotide detection is working correctly.")
        return 0
    else:
        print("\nSome tests failed. Please check the implementation.")
        return 1

if __name__ == "__main__":
    exit(main())
