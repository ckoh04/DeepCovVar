#!/usr/bin/python

"""
DeepCovVar: COVID-19 Variant Classification Tool

This is the main entry point for DeepCovVar, a deep learning-based tool for predicting 
and classifying COVID-19 variants. The tool operates in five phases:

Phase 1: Binary classification of sequences as virus or non-virus
Phase 2: Classification of (+)ssRNA (Class IV) vs others
Phase 3: Further classification of group IV ssRNA(+) (Coronavirus vs others)
Phase 4: SARS-CoV-2 vs other coronaviruses (MERS, SARS, Others)
Phase 5: SARS-CoV-2 variant classification (Omicron, Alpha, Delta, Others)

Author: AI Lab Team
Version: 1.0.0

Usage:
    python -m DeepCovVar [options] -f <fasta_file> -o <output_dir> -p <phase>
    python -m DeepCovVar [options] -f <fasta_file> -o <output_dir> --all-phases

Dependencies:
    - pandas
    - tensorflow
    - keras
    - torch
    - numpy
    - biopython
    - transformers
"""

# Standard library imports
import os
import sys
import time
import logging
import shutil
import argparse
from pathlib import Path
from typing import Tuple, Optional, Dict, Union

# Force CPU-only mode and suppress TensorFlow warnings - must be before importing TensorFlow
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'  # Force CPU usage
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # Suppress all TF logging
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'  # Disable oneDNN optimization warnings
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'  # Prevent CUDA memory errors
os.environ['XLA_FLAGS'] = '--xla_gpu_cuda_data_dir=/usr/local/cuda/'  # Suppress XLA warnings
os.environ['TF_XLA_FLAGS'] = '--tf_xla_enable_xla_devices=false'  # Disable XLA devices
os.environ['TF_ENABLE_EAGER_CLIENT_STREAMING_ENQUEUE'] = 'false'  # Disable eager execution warnings
os.environ['CUDA_DEVICE_ORDER'] = 'PCI_BUS_ID'  # Force PCI bus ID order
os.environ['CUDA_LAUNCH_BLOCKING'] = '1'  # Force synchronous CUDA execution
os.environ['TF_USE_CUDNN'] = '0'  # Disable cuDNN
os.environ['TF_DISABLE_CUDNN_RNN'] = '1'  # Disable cuDNN RNN ops
os.environ['TF_DISABLE_CUDNN_TENSOR_OP_MATH'] = '1'  # Disable cuDNN tensor ops

# Suppress future warnings
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
warnings.filterwarnings('ignore', message='.*experimental feature.*')
warnings.filterwarnings('ignore', message='.*cuda.*')
warnings.filterwarnings('ignore', message='.*gpu.*')

# Configure logging before imports
logging.getLogger('tensorflow').setLevel(logging.ERROR)
logging.getLogger('tensorflow').addHandler(logging.NullHandler())
logging.getLogger('tensorboard').setLevel(logging.ERROR)
logging.getLogger('tensorflow_hub').setLevel(logging.ERROR)
logging.getLogger('h5py').setLevel(logging.ERROR)
logging.getLogger('numexpr').setLevel(logging.ERROR)

# Third-party imports
import pandas as pd
import tensorflow as tf

# Configure TensorFlow to use CPU only - must be done before any TF operations
tf.config.set_visible_devices([], 'GPU')  # Hide all GPUs
tf.config.threading.set_inter_op_parallelism_threads(1)  # Limit to single thread
tf.config.threading.set_intra_op_parallelism_threads(1)  # Limit to single thread

# Additional TensorFlow warning suppression
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
tf.get_logger().setLevel('ERROR')
tf.autograph.set_verbosity(0)
tf.get_logger().addHandler(logging.NullHandler())

# Configure Numpy
import numpy as np
np.seterr(all='ignore')  # Suppress numpy warnings

# Disable TensorFlow debugging and performance logs
import tensorflow.python.util.deprecation as deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False

# Local imports
from deepcovvar import (
    covid_classifier,
    utils,
    __version__
)
from deepcovvar.covid_classifier import COVIDClassifier

# Additional TensorFlow CPU configuration
tf.keras.backend.set_floatx('float32')  # Use float32 for better CPU performance

def setup_logging(log_level: str = "INFO") -> None:
    """Set up logging configuration."""
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {log_level}")
    
    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler('deepcovvar.log')
        ]
    )

def validate_input_file(file_path: str) -> bool:
    """Validate that the input file exists and is a FASTA file."""
    path = Path(file_path)
    if not path.exists():
        print(f"Error: Input file '{file_path}' does not exist.")
        return False
    
    if not path.suffix.lower() in ['.fasta', '.fa', '.fas']:
        print(f"Warning: Input file '{file_path}' may not be a FASTA file.")
    
    return True

def validate_output_dir(output_dir: str) -> bool:
    """Validate and create output directory if it doesn't exist."""
    path = Path(output_dir)
    if not path.exists():
        try:
            path.mkdir(parents=True, exist_ok=True)
            print(f"Created output directory: {output_dir}")
        except Exception as e:
            print(f"Error creating output directory: {e}")
            return False
    return True

def main():
    """Main entry point for DeepCovVar."""
    parser = argparse.ArgumentParser(
        description='DeepCovVar: COVID-19 Variant Classification Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run a specific phase
  python -m deepcovvar -f input.fasta -o output_dir -p 5
  
  # Run all phases (recommended)
  python -m deepcovvar -f input.fasta -o output_dir --all-phases
  
  # Use default output directory (current directory)
  python -m deepcovvar -f input.fasta --all-phases
  
  # Run binary classification with custom thresholds
  python -m deepcovvar -f input.fasta -o output_dir -p 1 --thresholds 40 60
        """
    )
    
    parser.add_argument(
        '-f', '--fasta',
        required=True,
        help='Input FASTA file containing sequences'
    )
    
    parser.add_argument(
        '-o', '--output',
        default='.',
        help='Output directory for results (default: current directory)'
    )
    
    parser.add_argument(
        '-p', '--phase',
        type=int,
        choices=[1, 2, 3, 4, 5],
        help='Specific phase to run (1-5). Use --all-phases to run complete pipeline.'
    )
    
    parser.add_argument(
        '--all-phases',
        action='store_true',
        help='Run all phases from 1 to 5 (recommended)'
    )
    
    parser.add_argument(
        '--model-dir',
        default='models',
        help='Directory containing model files (default: models)'
    )
    
    parser.add_argument(
        '--thresholds',
        nargs=2,
        metavar=('CLASS1_THRESHOLD', 'CLASS2_THRESHOLD'),
        help='Custom thresholds for binary classification (e.g., --thresholds 40 60 for 40%% and 60%%)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version=f'DeepCovVar {__version__}'
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.phase and args.all_phases:
        parser.error("Cannot specify both --phase and --all-phases")
    
    if not args.phase and not args.all_phases:
        parser.error("Must specify either --phase or --all-phases")
    
    # Setup logging
    log_level = "DEBUG" if args.verbose else "INFO"
    setup_logging(log_level)
    
    logger = logging.getLogger(__name__)
    logger.info(f"DeepCovVar {__version__} starting...")
    
    # Validate inputs
    if not validate_input_file(args.fasta):
        sys.exit(1)
    
    if not validate_output_dir(args.output):
        sys.exit(1)
    
    try:
        # Initialize classifier
        model_dir = Path(args.model_dir)
        if not model_dir.is_absolute():
            model_dir = Path(__file__).parent / args.model_dir
        
        logger.info(f"Using model directory: {model_dir}")
        classifier = COVIDClassifier(model_dir=str(model_dir))
        
        if args.all_phases:
            # Run complete pipeline
            logger.info("Running complete COVID classification pipeline (all phases)")
            start_time = time.time()
            
            # Get base filename for output
            base_filename = Path(args.fasta).stem
            
            # Run all phases
            all_results = classifier.run_all_phases(
                input_file=args.fasta,
                output_dir=args.output,
                base_filename=base_filename
            )
            
            elapsed_time = time.time() - start_time
            logger.info(f"Complete pipeline completed in {elapsed_time:.2f} seconds")
            
            # Print summary
            print(f"\nDeepCovVar Complete Pipeline Finished!")
            print(f"Processed: {args.fasta}")
            print(f"Results saved to: {args.output}")
            print(f"Total time: {elapsed_time:.2f} seconds")
            
        else:
            # Run specific phase
            logger.info(f"Loading model for phase {args.phase}")
            classifier.load_model(args.phase)
            
            # Process sequences
            logger.info(f"Processing sequences from {args.fasta}")
            start_time = time.time()
            
            # Read sequences
            from Bio import SeqIO
            sequences = list(SeqIO.parse(args.fasta, "fasta"))
            logger.info(f"Loaded {len(sequences)} sequences")
            
            # Process each sequence
            results = []
            for i, seq_record in enumerate(sequences):
                logger.info(f"Processing sequence {i+1}/{len(sequences)}: {seq_record.id}")
                
                # Extract features and predict
                try:
                    # Use the correct predict method
                    # For binary classification phases, handle thresholds
                    config = classifier.models_config[args.phase]
                    custom_thresholds = None
                    if len(config['classes']) == 2:
                        if args.thresholds:
                            # Use command-line thresholds
                            try:
                                threshold1 = float(args.thresholds[0]) / 100.0
                                threshold2 = float(args.thresholds[1]) / 100.0
                                if 0.0 <= threshold1 <= 1.0 and 0.0 <= threshold2 <= 1.0:
                                    custom_thresholds = {
                                        config['classes'][0]: threshold1,
                                        config['classes'][1]: threshold2
                                    }
                                    print(f"Using command-line thresholds: {config['classes'][0]}={threshold1:.1%}, {config['classes'][1]}={threshold2:.1%}")
                                else:
                                    print("Warning: Thresholds must be between 0 and 100. Using interactive mode.")
                                    custom_thresholds = classifier._get_binary_thresholds(args.phase, config)
                            except ValueError:
                                print("Warning: Invalid threshold values. Using interactive mode.")
                                custom_thresholds = classifier._get_binary_thresholds(args.phase, config)
                        else:
                            # Use interactive mode
                            custom_thresholds = classifier._get_binary_thresholds(args.phase, config)
                    
                    prediction_result = classifier.predict(args.phase, args.fasta, custom_thresholds=custom_thresholds)
                    
                    # Extract prediction for this specific sequence
                    seq_prediction = prediction_result[prediction_result['Sequence_ID'] == seq_record.id]
                    if not seq_prediction.empty:
                        prediction = seq_prediction.iloc[0]['Predicted_Class']
                        confidence = seq_prediction.iloc[0]['Confidence']
                        # Convert percentage string back to float for logging
                        try:
                            confidence_float = float(confidence.strip('%')) / 100.0
                        except:
                            confidence_float = 0.0
                    else:
                        prediction = 'Unknown'
                        confidence = '0.00%'
                        confidence_float = 0.0
                    
                    results.append({
                        'sequence_id': seq_record.id,
                        'sequence': str(seq_record.seq),
                        'phase': args.phase,
                        'prediction': prediction,
                        'confidence': confidence
                    })
                    logger.info(f"Prediction for {seq_record.id}: {prediction} (confidence: {confidence})")
                except Exception as e:
                    logger.error(f"Error processing sequence {seq_record.id}: {e}")
                    results.append({
                        'sequence_id': seq_record.id,
                        'sequence': str(seq_record.seq),
                        'phase': args.phase,
                        'prediction': 'ERROR',
                        'error': str(e)
                    })
            
            # Save results
            output_file = Path(args.output) / f"phase_{args.phase}_results.csv"
            df = pd.DataFrame(results)
            df.to_csv(output_file, index=False)
            
            elapsed_time = time.time() - start_time
            logger.info(f"Processing completed in {elapsed_time:.2f} seconds")
            logger.info(f"Results saved to: {output_file}")
            
            # Print summary
            print(f"\nDeepCovVar Phase {args.phase} Classification Complete!")
            print(f"Processed {len(sequences)} sequences")
            print(f"Results saved to: {output_file}")
            print(f"Total time: {elapsed_time:.2f} seconds")
        
    except Exception as e:
        logger.error(f"Fatal error: {e}")
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
