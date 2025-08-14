#!/usr/bin/env python3
"""
Script to randomly select N sequences from a FASTA file
"""

import random
from Bio import SeqIO
import sys

def select_random_sequences(input_file, output_file, num_sequences):
    """Select random sequences from input FASTA file"""
    
    # Read all sequences
    print(f"Reading sequences from {input_file}...")
    sequences = list(SeqIO.parse(input_file, "fasta"))
    print(f"Found {len(sequences)} sequences")
    
    if num_sequences > len(sequences):
        print(f"Warning: Requested {num_sequences} sequences but only {len(sequences)} available")
        num_sequences = len(sequences)
    
    # Randomly select sequences
    print(f"Randomly selecting {num_sequences} sequences...")
    selected_indices = random.sample(range(len(sequences)), num_sequences)
    selected_sequences = [sequences[i] for i in selected_indices]
    
    # Write selected sequences
    print(f"Writing {len(selected_sequences)} sequences to {output_file}...")
    SeqIO.write(selected_sequences, output_file, "fasta")
    
    print(f"Successfully created {output_file} with {len(selected_sequences)} sequences")
    
    # Print selected sequence IDs
    print("\nSelected sequences:")
    for i, seq in enumerate(selected_sequences, 1):
        print(f"  {i:2d}. {seq.id}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python select_random_sequences.py <input_fasta> <output_fasta> <num_sequences>")
        print("Example: python select_random_sequences.py input.fasta output.fasta 50")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    num_sequences = int(sys.argv[3])
    
    # Set random seed for reproducibility
    random.seed(42)
    
    select_random_sequences(input_file, output_file, num_sequences)
