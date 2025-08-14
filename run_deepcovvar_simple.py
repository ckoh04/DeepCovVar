#!/usr/bin/env python3
"""
Simple wrapper script for DeepCovVar

This script provides a simple interface to run DeepCovVar with minimal arguments.
It automatically runs all phases and uses sensible defaults.

Usage:
    python run_deepcovvar_simple.py input.fasta
    python run_deepcovvar_simple.py input.fasta output_directory
    python run_deepcovvar_simple.py --help

Examples:
    # Run with default output directory (current directory)
    python run_deepcovvar_simple.py sequences.fasta
    
    # Run with custom output directory
    python run_deepcovvar_simple.py sequences.fasta results/
    
    # Run with verbose output
    python run_deepcovvar_simple.py --verbose sequences.fasta results/
"""

import sys
import os
from pathlib import Path

def main():
    if len(sys.argv) < 2 or '--help' in sys.argv or '-h' in sys.argv:
        print(__doc__)
        sys.exit(0)
    
    # Parse arguments
    verbose = '--verbose' in sys.argv
    if verbose:
        sys.argv.remove('--verbose')
    
    if len(sys.argv) == 2:
        # Only input file provided, use current directory as output
        input_file = sys.argv[1]
        output_dir = '.'
    elif len(sys.argv) == 3:
        # Both input file and output directory provided
        input_file = sys.argv[1]
        output_dir = sys.argv[2]
    else:
        print("Error: Invalid number of arguments")
        print(__doc__)
        sys.exit(1)
    
    # Validate input file
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Build the command
    cmd = f"python -m DeepCovVar -f {input_file} -o {output_dir} --all-phases"
    if verbose:
        cmd += " --verbose"
    
    print(f"Running DeepCovVar with command:")
    print(f"  {cmd}")
    print()
    
    # Execute the command
    exit_code = os.system(cmd)
    
    if exit_code == 0:
        print(f"\n✅ DeepCovVar completed successfully!")
        print(f"📁 Results saved to: {output_dir}")
    else:
        print(f"\n❌ DeepCovVar failed with exit code: {exit_code}")
        sys.exit(exit_code)

if __name__ == "__main__":
    main()
