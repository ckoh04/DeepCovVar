"""
Sequence Type Detection and Conversion Module for DeepCovVar

This module provides functionality to:
1. Detect whether input sequences are nucleotide or protein
2. Convert nucleotide sequences to protein sequences using Prodigal
3. Integrate seamlessly with the existing COVID classifier
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Tuple, List, Optional
import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SequenceTypeDetector:
    """Detects whether sequences are nucleotide or protein based on composition."""
    
    def __init__(self, threshold: float = 0.9):
        """
        Initialize the detector.
        
        Args:
            threshold: Minimum fraction of nucleotide bases to classify as nucleotide
        """
        self.threshold = threshold
        self.nt_bases = ['A', 'T', 'G', 'C', 'U', 'N']
        # Amino acid characters excluding those that overlap with DNA bases (A, T, G, C)
        # Only include unique amino acid characters: R, N, D, Q, E, H, I, L, K, M, F, P, S, W, Y, V
        self.aa_chars = 'RNDQEHILKMFPSTWYV'
    
    def is_nucleotide(self, seq: str) -> bool:
        """
        Determine if a sequence is likely nucleotide.
        
        Args:
            seq: Input sequence string
            
        Returns:
            True if sequence is likely nucleotide, False otherwise
        """
        seq = str(seq).upper()
        if not seq:
            return False
            
        # Count nucleotide bases
        nt_count = sum(seq.count(base) for base in self.nt_bases)
        nt_ratio = nt_count / len(seq)
        
        # Additional check: if sequence contains many amino acid characters, it's likely protein
        aa_count = sum(seq.count(aa) for aa in self.aa_chars)
        aa_ratio = aa_count / len(seq)
        
        # If amino acid ratio is high, it's likely protein
        if aa_ratio > 0.7:
            return False
            
        return nt_ratio >= self.threshold
    
    def check_fasta_type(self, fasta_file: str) -> List[Tuple[str, bool]]:
        """
        Check the type of each sequence in a FASTA file.
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            List of tuples (sequence_id, is_nucleotide)
        """
        results = []
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                is_nt = self.is_nucleotide(str(record.seq))
                results.append((record.id, is_nt))
                logger.info(f"{record.id}: {'Nucleotide' if is_nt else 'Protein'}")
        except Exception as e:
            logger.error(f"Error reading FASTA file: {e}")
            raise
            
        return results


class ProdigalConverter:
    """Converts nucleotide sequences to protein sequences using Prodigal."""
    
    def __init__(self, prodigal_path: Optional[str] = None):
        """
        Initialize the converter.
        
        Args:
            prodigal_path: Path to Prodigal executable. If None, will try to find it in PATH.
        """
        self.prodigal_path = prodigal_path or self._find_prodigal()
        if not self.prodigal_path:
            raise RuntimeError(
                "Prodigal not found. Please install Prodigal from "
                "https://github.com/hyattpd/Prodigal or provide the path."
            )
        logger.info(f"Using Prodigal at: {self.prodigal_path}")
    
    def _find_prodigal(self) -> Optional[str]:
        """Find Prodigal executable in system PATH."""
        try:
            result = subprocess.run(['which', 'prodigal'], 
                                  capture_output=True, text=True, check=True)
            return result.stdout.strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            return None
    
    def _check_prodigal_version(self) -> str:
        """Check Prodigal version."""
        try:
            result = subprocess.run([self.prodigal_path, '-v'], 
                                  capture_output=True, text=True, check=True)
            return result.stdout.strip()
        except subprocess.CalledProcessError:
            return "Unknown version"
    
    def convert_nucleotide_to_protein(self, 
                                    input_file: str, 
                                    output_file: str,
                                    mode: str = 'single') -> str:
        """
        Convert nucleotide sequences to protein sequences using Prodigal.
        
        Args:
            input_file: Input FASTA file with nucleotide sequences
            output_file: Output FASTA file for protein sequences
            mode: Prodigal mode ('single' for single genomes, 'meta' for metagenomes)
            
        Returns:
            Path to the output protein file
        """
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        # Prepare Prodigal command
        # Note: Older versions of Prodigal (like 2.6.3) don't support -f fasta for protein output
        # The -a option automatically outputs in FASTA format for protein sequences
        cmd = [
            self.prodigal_path,
            '-i', input_file,
            '-a', output_file,
            '-p', mode,
            '-q'  # Quiet mode
        ]
        
        logger.info(f"Running Prodigal: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info("Prodigal conversion completed successfully")
            return output_file
        except subprocess.CalledProcessError as e:
            logger.error(f"Prodigal failed: {e}")
            logger.error(f"Prodigal stderr: {e.stderr}")
            raise RuntimeError(f"Prodigal conversion failed: {e.stderr}")
    
    def convert_sequences_in_memory(self, 
                                  sequences: List[str], 
                                  seq_ids: List[str],
                                  mode: str = 'single') -> Tuple[List[str], List[str]]:
        """
        Convert nucleotide sequences to protein sequences in memory.
        
        Args:
            sequences: List of nucleotide sequences
            seq_ids: List of sequence IDs
            mode: Prodigal mode
            
        Returns:
            Tuple of (protein_sequences, protein_seq_ids)
        """
        # Create temporary input file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_in:
            for seq_id, seq in zip(seq_ids, sequences):
                temp_in.write(f">{seq_id}\n{seq}\n")
            temp_in_path = temp_in.name
        
        # Create temporary output file
        temp_out_path = temp_in_path.replace('.fasta', '_proteins.fasta')
        
        try:
            # Convert using Prodigal
            self.convert_nucleotide_to_protein(temp_in_path, temp_out_path, mode)
            
            # Read converted proteins
            protein_sequences = []
            protein_seq_ids = []
            
            for record in SeqIO.parse(temp_out_path, "fasta"):
                protein_sequences.append(str(record.seq))
                protein_seq_ids.append(record.id)
                
            return protein_sequences, protein_seq_ids
            
        finally:
            # Clean up temporary files
            for temp_file in [temp_in_path, temp_out_path]:
                if os.path.exists(temp_file):
                    os.unlink(temp_file)


class SequenceProcessor:
    """Main interface for sequence processing that integrates with DeepCovVar."""
    
    def __init__(self, prodigal_path: Optional[str] = None):
        """
        Initialize the sequence processor.
        
        Args:
            prodigal_path: Path to Prodigal executable
        """
        self.detector = SequenceTypeDetector()
        self.converter = ProdigalConverter(prodigal_path)
    
    def process_sequences(self, 
                         input_file: str, 
                         output_file: Optional[str] = None,
                         force_conversion: bool = False) -> Tuple[str, bool]:
        """
        Process sequences: detect type and convert if necessary.
        
        Args:
            input_file: Input FASTA file
            output_file: Output file for converted sequences (if conversion needed)
            force_conversion: Force conversion even if sequences appear to be protein
            
        Returns:
            Tuple of (processed_file_path, was_converted)
        """
        # Check sequence types
        sequence_types = self.detector.check_fasta_type(input_file)
        all_nucleotide = all(is_nt for _, is_nt in sequence_types)
        any_nucleotide = any(is_nt for _, is_nt in sequence_types)
        
        if not any_nucleotide and not force_conversion:
            logger.info("All sequences appear to be protein. No conversion needed.")
            return input_file, False
        
        if all_nucleotide:
            logger.info("All sequences are nucleotide. Converting to protein...")
        else:
            logger.warning("Mixed sequence types detected. Converting nucleotide sequences...")
        
        # Determine output file path
        if output_file is None:
            base_name = Path(input_file).stem
            output_file = f"{base_name}_converted_proteins.fasta"
        
        # Convert sequences
        try:
            self.converter.convert_nucleotide_to_protein(input_file, output_file)
            logger.info(f"Conversion completed. Protein sequences saved to: {output_file}")
            return output_file, True
        except Exception as e:
            logger.error(f"Conversion failed: {e}")
            raise
    
    def process_sequences_in_memory(self, 
                                  sequences: List[str], 
                                  seq_ids: List[str]) -> Tuple[List[str], List[str], bool]:
        """
        Process sequences in memory: detect type and convert if necessary.
        
        Args:
            sequences: List of sequences
            seq_ids: List of sequence IDs
            
        Returns:
            Tuple of (processed_sequences, processed_seq_ids, was_converted)
        """
        # Check if any sequences are nucleotide
        sequence_types = [(seq_id, self.detector.is_nucleotide(seq)) 
                         for seq_id, seq in zip(seq_ids, sequences)]
        
        any_nucleotide = any(is_nt for _, is_nt in sequence_types)
        
        if not any_nucleotide:
            logger.info("All sequences appear to be protein. No conversion needed.")
            return sequences, seq_ids, False
        
        # Convert nucleotide sequences
        logger.info("Converting nucleotide sequences to protein...")
        try:
            protein_sequences, protein_seq_ids = self.converter.convert_sequences_in_memory(
                sequences, seq_ids
            )
            logger.info(f"Conversion completed. {len(protein_sequences)} protein sequences generated.")
            return protein_sequences, protein_seq_ids, True
        except Exception as e:
            logger.error(f"Conversion failed: {e}")
            raise


def main():
    """Command-line interface for testing the sequence processor."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Sequence Type Detection and Conversion Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('--output', '-o', help='Output file for converted sequences')
    parser.add_argument('--prodigal-path', help='Path to Prodigal executable')
    parser.add_argument('--check-only', action='store_true', 
                       help='Only check sequence types, do not convert')
    parser.add_argument('--force', action='store_true',
                       help='Force conversion even if sequences appear to be protein')
    
    args = parser.parse_args()
    
    try:
        processor = SequenceProcessor(args.prodigal_path)
        
        if args.check_only:
            # Just check sequence types
            processor.detector.check_fasta_type(args.input_file)
        else:
            # Process and convert if needed
            output_file, was_converted = processor.process_sequences(
                args.input_file, args.output, args.force
            )
            
            if was_converted:
                print(f"\nConversion completed successfully!")
                print(f"Protein sequences saved to: {output_file}")
            else:
                print(f"\nNo conversion needed. Sequences are already protein.")
                
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
