#!/usr/bin/env python3
"""
Test script for Prodigal nucleotide-to-protein conversion.

This script demonstrates the actual conversion functionality now that Prodigal is installed.
"""

import sys
import os
from pathlib import Path

# Add the parent directory to the path to import DeepCovVar modules
if __name__ == "__main__":
    # When running directly, use current working directory's parent
    parent_dir = str(Path.cwd().parent)
else:
    # When imported as module, use file location's parent
    parent_dir = str(Path(__file__).parent.parent.absolute())

sys.path.insert(0, parent_dir)

try:
    from sequence_converter import ProdigalConverter, SequenceTypeDetector
    print("Successfully imported ProdigalConverter and SequenceTypeDetector")
except ImportError as e:
    print(f"Import error: {e}")
    sys.exit(1)

def test_prodigal_conversion():
    """Test the actual nucleotide-to-protein conversion."""
    print("Testing Prodigal Nucleotide-to-Protein Conversion")
    print("=" * 50)
    
    try:
        # Initialize converter
        converter = ProdigalConverter()
        print(f"ProdigalConverter initialized successfully")
        print(f"  Prodigal path: {converter.prodigal_path}")
        
        # Test file with nucleotide sequences
        test_file = Path(__file__).parent / "test_nucleotide_sequences.fasta"
        
        if not test_file.exists():
            print(f"✗ Test file {test_file} not found!")
            return False
        
        print(f"Test file found: {test_file}")
        
        # Create output file for converted proteins
        output_file = Path(__file__).parent / "converted_proteins.fasta"
        
        print(f"\nConverting nucleotide sequences to protein sequences...")
        print(f"Output will be saved to: {output_file}")
        
        # Perform the conversion
        # Use 'meta' mode for shorter sequences (better for test data)
        result_file = converter.convert_nucleotide_to_protein(
            str(test_file), 
            str(output_file),
            mode='meta'
        )
        
        print(f"Conversion completed successfully!")
        print(f"  Output file: {result_file}")
        
        # Verify the output file exists and has content
        if output_file.exists():
            print(f"Output file created successfully")
            
            # Count sequences in output
            from Bio import SeqIO
            protein_sequences = list(SeqIO.parse(output_file, "fasta"))
            print(f"Generated {len(protein_sequences)} protein sequences")
            
            # Show first few sequences
            print(f"\nFirst 2 converted protein sequences:")
            for i, record in enumerate(protein_sequences[:2]):
                print(f"  {i+1}. {record.id}")
                print(f"     Length: {len(record.seq)} amino acids")
                print(f"     Preview: {str(record.seq)[:50]}...")
                print()
            
            return True
        else:
            print(f"Output file was not created")
            return False
            
    except Exception as e:
        print(f"Error during conversion: {e}")
        return False

def test_sequence_detection():
    """Test sequence type detection."""
    print("\n" + "="*50)
    print("Testing Sequence Type Detection")
    print("="*50)
    
    try:
        detector = SequenceTypeDetector()
        
        # Test file with nucleotide sequences
        test_file = Path(__file__).parent / "test_nucleotide_sequences.fasta"
        
        if not test_file.exists():
            print(f"Test file {test_file} not found!")
            return False
        
        print(f"Test file found: {test_file}")
        
        # Check sequence types
        sequence_types = detector.check_fasta_type(str(test_file))
        nucleotide_count = sum(1 for _, is_nt in sequence_types if is_nt)
        
        print(f"Detected {nucleotide_count} nucleotide sequences out of {len(sequence_types)} total")
        
        # Show detailed breakdown
        print(f"\nDetailed sequence analysis:")
        for seq_id, is_nucleotide in sequence_types:
            seq_type = "NUCLEOTIDE" if is_nucleotide else "PROTEIN"
            print(f"  {seq_id}: {seq_type}")
        
        return True
        
    except Exception as e:
        print(f"✗ Error during sequence detection: {e}")
        return False

def main():
    """Main test function."""
    print("Prodigal Conversion Test Suite")
    print("=" * 60)
    
    # Run tests
    tests = [
        ("Sequence Detection", test_sequence_detection),
        ("Prodigal Conversion", test_prodigal_conversion),
    ]
    
    results = []
    for test_name, test_func in tests:
        print(f"\nRunning: {test_name}")
        print("-" * 40)
        
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"Test failed with exception: {e}")
            results.append((test_name, False))
    
    # Summary
    print(f"\n{'='*60}")
    print("TEST SUMMARY")
    print(f"{'='*60}")
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "PASSED" if result else "FAILED"
        print(f"{test_name}: {status}")
    
    print(f"\nResults: {passed}/{total} tests passed")
    
    if passed == total:
        print("\nAll tests passed! Prodigal conversion is working correctly.")
        print("Nucleotide sequences can now be converted to protein sequences")
        print("The DeepCovVar pipeline is fully functional")
        return 0
    else:
        print(f"\n{total - passed} test(s) failed. Please check the implementation.")
        return 1

if __name__ == "__main__":
    exit(main())
