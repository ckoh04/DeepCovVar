def fix_numpy_compatibility():
    try:
        import numpy as np
        if not hasattr(np, 'typeDict'):
            np.typeDict = np.sctypeDict
    except Exception:
        pass  

fix_numpy_compatibility()

import argparse
import os
import sys
import numpy as np
from pathlib import Path
import tensorflow as tf
import torch
import torch.nn as nn
from Bio import SeqIO
import pandas as pd
import json
import pkg_resources
from transformers import AutoTokenizer, AutoModelForSequenceClassification

from .utils import FEATURE, preprocess, preprocessdf
from .utils import SequenceProcessor

class TransformerModel(nn.Module):
    def __init__(self, vocab_size, d_model=512, nhead=8, num_layers=6, num_classes=3):
        super(TransformerModel, self).__init__()    
        self.embedding = nn.Embedding(vocab_size, d_model)
        self.transformer = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(d_model=d_model, nhead=nhead, batch_first=True),
            num_layers=num_layers
        )
        self.pooler = nn.AdaptiveAvgPool1d(1)
        self.classifier = nn.Linear(d_model, num_classes)
        
    def forward(self, input_ids):
        embedded = self.embedding(input_ids)
        transformer_output = self.transformer(embedded)
        
        pooled = self.pooler(transformer_output.transpose(1, 2)).squeeze(-1)  
        logits = self.classifier(pooled)  
        return logits

class COVIDClassifier:
    def __init__(self, model_dir=None, prodigal_path=None, batch_size=32):
        # Use pkg_resources to get model directory from installed package
        if model_dir is None:
            try:
                self.model_dir = Path(pkg_resources.resource_filename('deepcovvar', 'models'))
            except Exception:
                # Fallback for development when package is not installed
                self.model_dir = Path(__file__).parent / "models"
        else:
            self.model_dir = Path(model_dir)
        self.feature_extractor = FEATURE()
        self.sequence_processor = SequenceProcessor(prodigal_path)
        self.batch_size = batch_size  # Configurable batch size for memory management
        
        self.models_config = {
            1: {
                        'file': 'p1_40_final_model_quantized.keras',
        'description': 'Virus sequences vs Others',
                'classes': ['Virus', 'Non-virus'],
                'feature_size': 2400,  # CKSAAP with gap=5: (5+1)*400 = 2400
                'type': 'keras'
            },
            2: {
                        'file': 'p2_final_model_quantized.keras',
        'description': '(+)ssRNA (Class IV) vs Others',
                'classes': ['Others', '(+)ssRNA'],
                'feature_size': 2400,
                'type': 'keras'
            },
            3: {
                        'file': 'p3_final_model_quantized.keras',
        'description': 'Further classification of group IV ssRNA(+)',
                'classes': ['Other ssRNA(+)', 'Coronavirus'],
                'feature_size': 2400,
                'type': 'keras'
            },
            4: {
                        'file': 'p4_final_model_quantized.keras',
        'description': 'SARS-CoV-2 vs other coronaviruses',
                'classes': ['SARS-CoV-2', 'MERS', 'SARS', 'Others'],
                'feature_size': 2400,
                'type': 'keras'
            },
            5: {
                        'file': 'p5_final_model_quantized.pt',
        'description': 'SARS-CoV-2 variant classification',
                'classes': ['Omicron', 'Alpha', 'Delta', 'Epsilon', 'Iota', 'Gamma', 'Others'], 
                'feature_size': None, 
                'type': 'pytorch_transformer'
            }
        }
        
        self.loaded_models = {}
        self.tokenizer = None
        
        # Get model paths using pkg_resources
        self._update_model_paths()
    
    def _update_model_paths(self):
        """Update model file paths using pkg_resources for installed package."""
        try:
            for phase, config in self.models_config.items():
                if config['type'] == 'keras':
                    # For Keras models, use pkg_resources
                    try:
                        config['_full_path'] = pkg_resources.resource_filename('deepcovvar', f"models/{config['file']}")
                    except Exception:
                        # Fallback to relative path for development
                        config['_full_path'] = str(self.model_dir / config['file'])
                elif config['type'] == 'pytorch_transformer':
                    # For PyTorch models, also use pkg_resources
                    try:
                        config['_full_path'] = pkg_resources.resource_filename('deepcovvar', f"models/{config['file']}")
                    except Exception:
                        # Fallback to relative path for development
                        config['_full_path'] = str(self.model_dir / config['file'])
        except Exception as e:
            print(f"Warning: Could not update model paths with pkg_resources: {e}")
            # Fallback: use relative paths
            for phase, config in self.models_config.items():
                config['_full_path'] = str(self.model_dir / config['file'])
    
    def load_transformer_model(self, model_path):
        model_dir = Path(model_path)
        
        # Try to get config file path using pkg_resources first
        try:
            config_path = Path(pkg_resources.resource_filename('deepcovvar', 'models/config.pt'))
        except Exception:
            config_path = model_dir / "config.pt"
        
        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")
        
        config = torch.load(config_path, map_location='cpu')
        print(f"Loaded model config: {config}")
        
        
        base_model_name = config.get('base_model_name', 'facebook/esm2_t33_650M_UR50D')
        num_classes = config.get('num_classes', 7)
        
        print(f"Loading ESM-2 model: {base_model_name}")
        

        try:
            self.tokenizer = AutoTokenizer.from_pretrained(base_model_name)
            print("Tokenizer loaded successfully")
        except Exception as e:
            print(f"Warning: Could not load tokenizer from {base_model_name}, using fallback")
            self.tokenizer = None
        
        try:
            
            print("Loading base ESM-2 model...")
            model = AutoModelForSequenceClassification.from_pretrained(
                base_model_name,
                num_labels=num_classes,
                ignore_mismatched_sizes=True
            )
            
            # Try to load state dict first, fall back to original if not found
            try:
                state_dict_path = Path(pkg_resources.resource_filename('deepcovvar', 'models/model_state_dict_quantized.pt'))
            except Exception:
                state_dict_path = model_dir / "model_state_dict_quantized.pt"
                
            if not state_dict_path.exists():
                try:
                    state_dict_path = Path(pkg_resources.resource_filename('deepcovvar', 'models/model_state_dict.pt'))
                except Exception:
                    state_dict_path = model_dir / "model_state_dict.pt"
                print("Model not found, using original model")
            else:
                print("Loading model weights...")
            
            if not state_dict_path.exists():
                raise FileNotFoundError(f"Model state dict not found: {state_dict_path}")
            
            print("Loading trained weights...")
            state_dict = torch.load(state_dict_path, map_location='cpu')
            
            missing_keys, unexpected_keys = model.load_state_dict(state_dict, strict=False)
            
            if missing_keys:
                print(f"Warning: Missing keys: {missing_keys}")
            if unexpected_keys:
                print(f"Warning: Unexpected keys: {unexpected_keys}")
            
            model.eval()
            print("ESM-2 model loaded successfully")
            return model
            
        except Exception as e:
            print(f"Error loading ESM-2 model: {e}")
            raise
    
    def tokenize_sequences(self, sequences, max_length=512):
        if self.tokenizer is None:
            raise ValueError("Tokenizer not initialized. Load model first.")
        
        encoded = self.tokenizer(
            sequences,
            padding=True,
            truncation=True,
            max_length=max_length,
            return_tensors='pt'
        )
        
        return encoded['input_ids'], encoded['attention_mask']
    
    def load_model(self, phase):
        if phase in self.loaded_models:
            return self.loaded_models[phase]
        
        config = self.models_config[phase]
        
        if config['type'] == 'pytorch_transformer':
            
            # For Phase 5, we need to load from the models directory where config.pt and other files are located
            model_path = self.model_dir
            
            if not model_path.exists():
                raise FileNotFoundError(f"Model directory not found: {model_path}")
            
            print(f"Loading transformer model for Phase {phase}: {config['description']}")
            model = self.load_transformer_model(model_path)
        else:
            # Use the full path from pkg_resources if available
            if '_full_path' in config:
                model_path = Path(config['_full_path'])
            else:
                model_path = self.model_dir / config['file']
            
            # If model doesn't exist, fall back to original
            if not model_path.exists():
                original_file = config['file'].replace('_quantized.keras', '.keras')
                if '_full_path' in config:
                    # Try to get original file path from pkg_resources
                    try:
                        original_path = pkg_resources.resource_filename('deepcovvar', f"models/{original_file}")
                        model_path = Path(original_path)
                    except Exception:
                        # Fallback to relative path
                        model_path = self.model_dir / original_file
                else:
                    model_path = self.model_dir / original_file
                    
                if model_path.exists():
                    print(f"Model not found, using original: {original_file}")
                    config['description'] = config['description'].replace(' (Quantized)', '')
                else:
                    raise FileNotFoundError(f"Model file not found: {model_path}")
            else:
                print(f"Using model: {config['file']}")
            
            print(f"Loading model for Phase {phase}: {config['description']}")
                
            if config['type'] == 'keras':
                model = tf.keras.models.load_model(str(model_path))
            elif config['type'] == 'pytorch':
                model = torch.load(str(model_path), map_location='cpu')
                model.eval()
        
        self.loaded_models[phase] = model
        return model
    
    def read_sequences(self, input_file, auto_convert=True):
        """
        Read sequences from FASTA file with automatic nucleotide detection and conversion.
        
        Args:
            input_file: Path to input FASTA file
            auto_convert: Whether to automatically convert nucleotide sequences to protein
            
        Returns:
            Tuple of (sequences, sequence_ids)
        """
        sequences = []
        seq_ids = []
        
        try:
            for record in SeqIO.parse(input_file, "fasta"):
                sequences.append(str(record.seq).upper())
                seq_ids.append(record.id)
            
            if not sequences:
                raise ValueError("No sequences found in the input file")
            
            print(f"Successfully read {len(sequences)} sequences from {input_file}")
            
            # Check if conversion is needed
            if auto_convert:
                try:
                    converted_sequences, converted_seq_ids, was_converted = \
                        self.sequence_processor.process_sequences_in_memory(sequences, seq_ids)
                    
                    if was_converted:
                        print(f"Converted {len(converted_sequences)} nucleotide sequences to protein sequences")
                        return converted_sequences, converted_seq_ids
                    else:
                        print("All sequences are already protein. No conversion needed.")
                        return sequences, seq_ids
                        
                except Exception as e:
                    print(f"Warning: Sequence conversion failed: {e}")
                    print("Proceeding with original sequences...")
                    return sequences, seq_ids
            
            return sequences, seq_ids
            
        except Exception as e:
            raise ValueError(f"Error reading sequences from {input_file}: {str(e)}")
    
    def extract_features(self, sequences, feature_size=2400):
        print("Extracting CKSAAP features...")
        
        features = []
        for i, seq in enumerate(sequences):
            try:
                clean_seq = ''.join([aa for aa in seq if aa in 'ARNDCQEGHILKMFPSTWYV'])
                
                if len(clean_seq) < 2:
                    raise ValueError(f"Sequence {i+1} too short after cleaning")
                

                cksaap_features = self.feature_extractor.CKSAAP(clean_seq, gap=5)
                
                if len(cksaap_features) != feature_size:
                    if len(cksaap_features) < feature_size:
                        cksaap_features.extend([0.0] * (feature_size - len(cksaap_features)))
                    else:
                        cksaap_features = cksaap_features[:feature_size]
                
                features.append(cksaap_features)
                
            except Exception as e:
                print(f"Warning: Error processing sequence {i+1}: {str(e)}")
                features.append([0.0] * feature_size)
        
        return np.array(features, dtype=np.float32)
    
    def process_input_sequences(self, input_file, output_file=None, force_conversion=False):
        """
        Process input sequences: detect type and convert if necessary.
        
        Args:
            input_file: Input FASTA file
            output_file: Output file for converted sequences
            force_conversion: Force conversion even if sequences appear to be protein
            
        Returns:
            Tuple of (processed_file_path, was_converted)
        """
        try:
            return self.sequence_processor.process_sequences(
                input_file, output_file, force_conversion
            )
        except Exception as e:
            print(f"Error processing sequences: {e}")
            raise
    
    def predict(self, phase, input_file, output_file=None, custom_thresholds=None):
        model = self.load_model(phase)
        config = self.models_config[phase]
        
        # Handle custom thresholds for binary classification
        if custom_thresholds is None and len(config['classes']) == 2:
            custom_thresholds = self._get_binary_thresholds(phase, config)
        
        sequences, seq_ids = self.read_sequences(input_file)

        if config['type'] == 'pytorch_transformer':
            features = None 
        else:
            features = self.extract_features(sequences, config['feature_size'])
        
        print(f"Making predictions using Phase {phase} model...")
        
        if config['type'] == 'keras':
            # Check if model expects 4D input (e.g., (None, 1, 2400, 1))
            expected_shape = None
            try:
                expected_shape = model.input_shape
            except Exception:
                pass
            if expected_shape and len(expected_shape) == 4:
                # Reshape features to (batch, 1, 2400, 1)
                features = features.reshape((features.shape[0], 1, features.shape[1], 1))
            predictions = model.predict(features, verbose=0)
            
            # Convert probabilities to class predictions
            if predictions.shape[1] > 1:
                predicted_classes = np.argmax(predictions, axis=1)
                confidence_scores = np.max(predictions, axis=1)
            else:
                # Binary classification with custom thresholds
                if custom_thresholds:
                    # Use custom threshold for first class (index 0)
                    threshold = custom_thresholds[config['classes'][0]]
                    predicted_classes = (predictions >= threshold).astype(int).flatten()
                    confidence_scores = np.where(predicted_classes == 1, 
                                               predictions.flatten(), 
                                               1 - predictions.flatten())
                else:
                    # Default 50% threshold
                    predicted_classes = (predictions > 0.5).astype(int).flatten()
                    confidence_scores = np.where(predicted_classes == 1, 
                                               predictions.flatten(), 
                                               1 - predictions.flatten())
                
        elif config['type'] == 'pytorch':
            with torch.no_grad():
                features_tensor = torch.FloatTensor(features)
                outputs = model(features_tensor)
                
                if hasattr(torch, 'softmax'):
                    predictions = torch.softmax(outputs, dim=1).numpy()
                else:
                    predictions = torch.nn.functional.softmax(outputs, dim=1).numpy()
                
                predicted_classes = np.argmax(predictions, axis=1)
                confidence_scores = np.max(predictions, axis=1)
                
        elif config['type'] == 'pytorch_transformer':
            # For ESM-2 transformer model, we need to tokenize sequences instead of using CKSAAP features
            print("Tokenizing sequences for ESM-2 transformer model...")
            
            # Process in batches to avoid memory issues
            batch_size = self.batch_size  # Use configurable batch size
            predicted_classes = []
            confidence_scores = []
            all_predictions = []  # Store all predictions for probability display
            
            print(f"Processing {len(sequences)} sequences in batches of {batch_size}")
            
            for i in range(0, len(sequences), batch_size):
                batch_end = min(i + batch_size, len(sequences))
                batch_sequences = sequences[i:batch_end]
                batch_seq_ids = seq_ids[i:batch_end]
                
                print(f"Processing batch {i//batch_size + 1}/{(len(sequences) + batch_size - 1)//batch_size} "
                      f"(sequences {i+1}-{batch_end})")
                
                try:
                    # Tokenize batch
                    input_ids, attention_mask = self.tokenize_sequences(batch_sequences)
                    
                    with torch.no_grad():
                        outputs = model(input_ids=input_ids, attention_mask=attention_mask)
                        logits = outputs.logits if hasattr(outputs, 'logits') else outputs
                        batch_predictions = torch.softmax(logits, dim=1).numpy()
                        
                        batch_pred_classes = np.argmax(batch_predictions, axis=1)
                        batch_conf_scores = np.max(batch_predictions, axis=1)
                        
                        predicted_classes.extend(batch_pred_classes)
                        confidence_scores.extend(batch_conf_scores)
                        all_predictions.extend(batch_predictions)  # Store all predictions
                        
                        # Clear GPU memory if available
                        if torch.cuda.is_available():
                            torch.cuda.empty_cache()
                            
                except Exception as e:
                    print(f"Error processing batch {i//batch_size + 1}: {e}")
                    # Fill with default values for failed batch
                    batch_size_actual = batch_end - i
                    predicted_classes.extend([0] * batch_size_actual)
                    confidence_scores.extend([0.0] * batch_size_actual)
                    # Create dummy predictions for failed batch
                    dummy_predictions = np.zeros((batch_size_actual, len(config['classes'])))
                    dummy_predictions[:, 0] = 1.0  # Set first class to 100%
                    all_predictions.extend(dummy_predictions)
            
            predicted_classes = np.array(predicted_classes)
            confidence_scores = np.array(confidence_scores)
            all_predictions = np.array(all_predictions)
    
        results = []
        for i, (seq_id, pred_class, confidence) in enumerate(zip(seq_ids, predicted_classes, confidence_scores)):
            class_name = config['classes'][pred_class] if pred_class < len(config['classes']) else 'Unknown'
            
            # Get probabilities for all classes
            if config['type'] == 'keras':
                if predictions.shape[1] > 1:
                    # Multi-class: get probabilities for all classes
                    all_probs = predictions[i]
                else:
                    # Binary: create [1-pred, pred] for [class0, class1]
                    pred_val = predictions[i][0] if len(predictions[i]) == 1 else predictions[i][0]
                    all_probs = [1 - pred_val, pred_val]
            elif config['type'] == 'pytorch':
                all_probs = predictions[i]
            elif config['type'] == 'pytorch_transformer':
                all_probs = all_predictions[i]
            
            # Create result row with all class probabilities
            result_row = {
                'Sequence_ID': seq_id,
                'Predicted_Class': class_name,
                'Confidence': f"{confidence:.2%}"
            }
            
            # Add probabilities for each class
            for j, class_name_j in enumerate(config['classes']):
                if j < len(all_probs):
                    prob_value = all_probs[j]
                    result_row[f'{class_name_j}_Probability'] = f"{prob_value:.2%}"
                else:
                    result_row[f'{class_name_j}_Probability'] = "0.00%"
            
            results.append(result_row)
        

        results_df = pd.DataFrame(results)
        

        print(f"\n{'='*60}")
        print(f"PHASE {phase} RESULTS: {config['description']}")
        print(f"{'='*60}")
        print(results_df.to_string(index=False))
        

        if output_file:
            results_df.to_csv(output_file, index=False)
            print(f"\nResults saved to: {output_file}")
        
        return results_df
    
    def _get_binary_thresholds(self, phase, config):
        """
        Get custom thresholds for binary classification phases from user input.
        
        Args:
            phase: Phase number
            config: Model configuration
            
        Returns:
            Dictionary mapping class names to thresholds, or None if user chooses default
        """
        if len(config['classes']) != 2:
            return None
            
        print(f"\n{'='*60}")
        print(f"BINARY CLASSIFICATION THRESHOLD SETTING - PHASE {phase}")
        print(f"{'='*60}")
        print(f"Phase: {config['description']}")
        print(f"Classes: {config['classes'][0]} vs {config['classes'][1]}")
        print(f"\nCurrent default threshold: 50% (0.5)")
        print(f"Example: If you set {config['classes'][0]} threshold to 40%, then")
        print(f"any sequence with ≥40% {config['classes'][0]} probability will be classified as {config['classes'][0]}")
        print(f"regardless of the {config['classes'][1]} probability.")
        
        while True:
            choice = input(f"\nDo you want to set custom thresholds? (y/n): ").lower().strip()
            if choice in ['y', 'yes']:
                break
            elif choice in ['n', 'no']:
                print(f"Using default 50% threshold for Phase {phase}")
                return None
            else:
                print("Please enter 'y' or 'n'")
        
        thresholds = {}
        for i, class_name in enumerate(config['classes']):
            while True:
                try:
                    threshold_input = input(f"Enter threshold for {class_name} (0-100%): ").strip()
                    if threshold_input.endswith('%'):
                        threshold_input = threshold_input[:-1]
                    
                    threshold = float(threshold_input) / 100.0
                    if 0.0 <= threshold <= 1.0:
                        thresholds[class_name] = threshold
                        print(f"Set {class_name} threshold to {threshold:.1%}")
                        break
                    else:
                        print("Threshold must be between 0% and 100%")
                except ValueError:
                    print("Please enter a valid number")
        
        print(f"\nCustom thresholds set for Phase {phase}:")
        for class_name, threshold in thresholds.items():
            print(f"  {class_name}: {threshold:.1%}")
        
        return thresholds
    
    def run_all_phases(self, input_file, output_dir=None, base_filename=None):
        """
        Run all phases of the COVID classifier pipeline and save results separately.
        
        Args:
            input_file: Input FASTA file
            output_dir: Directory to save results (default: current directory)
            base_filename: Base filename for output files (default: input filename without extension)
            
        Returns:
            Dictionary containing results from all phases
        """
        if output_dir is None:
            output_dir = Path.cwd()
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(exist_ok=True)
        
        if base_filename is None:
            base_filename = Path(input_file).stem
        
        all_results = {}
        
        print(f"\n{'='*80}")
        print("RUNNING COMPLETE COVID CLASSIFICATION PIPELINE")
        print(f"{'='*80}")
        print(f"Input file: {input_file}")
        print(f"Output directory: {output_dir}")
        print(f"Base filename: {base_filename}")
        print(f"{'='*80}\n")
        
        # Process sequences once (nucleotide conversion if needed)
        try:
            processed_file, was_converted = self.process_input_sequences(input_file)
            if was_converted:
                print(f"Using converted protein sequences from: {processed_file}")
                working_file = processed_file
            else:
                working_file = input_file
        except Exception as e:
            print(f"Warning: Sequence processing failed: {e}")
            print("Proceeding with original input file...")
            working_file = input_file
        
        # Run through all phases
        for phase in sorted(self.models_config.keys()):
            try:
                print(f"\n{'='*60}")
                print(f"RUNNING PHASE {phase}: {self.models_config[phase]['description']}")
                print(f"{'='*60}")
                
                # Generate output filename for this phase
                phase_output_file = output_dir / f"{base_filename}_phase_{phase}_results.csv"
                
                # Run prediction for this phase
                results_df = self.predict(phase, working_file, str(phase_output_file))
                all_results[phase] = results_df
                
                print(f"Phase {phase} completed successfully!")
                print(f"Results saved to: {phase_output_file}")
                
            except Exception as e:
                print(f"Error in Phase {phase}: {e}")
                all_results[phase] = None
                continue
        
        # Generate summary report
        summary_file = output_dir / f"{base_filename}_pipeline_summary.txt"
        self._generate_pipeline_summary(all_results, summary_file, input_file)
        
        print(f"\n{'='*80}")
        print("PIPELINE COMPLETED")
        print(f"{'='*80}")
        print(f"Summary report saved to: {summary_file}")
        print(f"Individual phase results saved to: {output_dir}")
        print(f"{'='*80}")
        
        return all_results
    
    def _generate_pipeline_summary(self, all_results, summary_file, input_file):
        """Generate a summary report of all pipeline phases."""
        with open(summary_file, 'w') as f:
            f.write("COVID CLASSIFICATION PIPELINE SUMMARY REPORT\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Input file: {input_file}\n")
            f.write(f"Generated: {pd.Timestamp.now()}\n\n")
            
            f.write("PHASE RESULTS SUMMARY:\n")
            f.write("-" * 30 + "\n")
            
            for phase, results in all_results.items():
                if results is not None:
                    config = self.models_config[phase]
                    f.write(f"\nPhase {phase}: {config['description']}\n")
                    f.write(f"Status: Completed\n")
                    f.write(f"Classes: {', '.join(config['classes'])}\n")
                    f.write(f"Sequences processed: {len(results)}\n")
                    
                    # Count predictions for each class
                    class_counts = results['Predicted_Class'].value_counts()
                    f.write("Predictions:\n")
                    for class_name, count in class_counts.items():
                        f.write(f"  {class_name}: {count}\n")
                else:
                    f.write(f"\nPhase {phase}: FAILED\n")
                    f.write(f"Status: Error occurred during processing\n")
            
            f.write(f"\n{'='*50}\n")
            f.write("End of Report\n")
    
    def interactive_mode(self, input_file):

        print("\n" + "="*60)
        print("DeepCovVar - INTERACTIVE MODE")
        print("="*60)
        
        print("\nAvailable options:")
        print("  0. Run ALL phases (complete pipeline)")
        for phase, config in self.models_config.items():
            print(f"  {phase}. {config['description']}")
        
        while True:
            try:
                choice = int(input(f"\nSelect option (0-{len(self.models_config)}): "))
                if choice == 0 or choice in self.models_config:
                    break
                else:
                    print(f"Invalid choice. Please select from 0-{len(self.models_config)}")
            except ValueError:
                print("Please enter a valid number")
        
        if choice == 0:
            # Run all phases
            print("\nRunning complete pipeline...")
            output_dir = input("Output directory (press Enter for current directory): ").strip()
            if not output_dir:
                output_dir = None
            
            base_filename = input("Base filename for results (press Enter for auto): ").strip()
            if not base_filename:
                base_filename = None
            
            return self.run_all_phases(input_file, output_dir, base_filename)
        else:
            # Run single phase
            phase = choice
            save_output = input("\nSave results to file? (y/n): ").lower().startswith('y')
            output_file = None
            if save_output:
                default_output = f"phase_{phase}_results.csv"
                output_file = input(f"Output filename [{default_output}]: ").strip()
                if not output_file:
                    output_file = default_output
            
            # For binary classification phases, ask about thresholds
            config = self.models_config[phase]
            custom_thresholds = None
            if len(config['classes']) == 2:
                custom_thresholds = self._get_binary_thresholds(phase, config)
            
            return self.predict(phase, input_file, output_file, custom_thresholds)


def main():
    parser = argparse.ArgumentParser(
        description="DeepCovVar - COVID-19 Variant Classifier with Automatic Nucleotide Conversion",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog= """
                    Examples:
                    # Interactive mode with automatic nucleotide conversion
                    python covid_classifier.py input_sequences.fasta
                    
                    # Direct phase selection with automatic conversion
                    python covid_classifier.py input_sequences.fasta --phase 1
                    
                    # With output file and automatic conversion
                    python covid_classifier.py input_sequences.fasta --phase 1 --output results.csv
                    
                    # Specify Prodigal path for nucleotide conversion
                    python covid_classifier.py input_sequences.fasta --prodigal-path /usr/local/bin/prodigal
                    
                    # Disable automatic conversion
                    python covid_classifier.py input_sequences.fasta --no-auto-convert
                    
                    # Force conversion even if sequences appear to be protein
                    python covid_classifier.py input_sequences.fasta --force-convert
                    
                    # Run complete pipeline (all phases)
                    python covid_classifier.py input_sequences.fasta --run-all
                    
                    # Run complete pipeline with custom output directory
                    python covid_classifier.py input_sequences.fasta --run-all --output-dir results/
                """
    )
    
    parser.add_argument('input_file', 
                       help='Input FASTA file containing protein sequences (.fa, .faa, .fasta)')
    
    parser.add_argument('--phase', '-p', 
                       type=int, 
                       choices=[1, 2, 3, 4, 5],
                       help='Phase/model to use (1-5)')
    
    parser.add_argument('--output', '-o',
                       help='Output CSV file for results')
    
    parser.add_argument('--model-dir', '-m',
                       default="models",
                       help='Directory containing model files')
    
    parser.add_argument('--prodigal-path', 
                       help='Path to Prodigal executable for nucleotide conversion')
    
    parser.add_argument('--no-auto-convert', action='store_true',
                       help='Disable automatic nucleotide to protein conversion')
    
    parser.add_argument('--force-convert', action='store_true',
                       help='Force conversion even if sequences appear to be protein')
    
    parser.add_argument('--run-all', action='store_true',
                       help='Run all phases of the classification pipeline')
    
    parser.add_argument('--output-dir', '-d',
                       help='Output directory for pipeline results (when using --run-all)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found")
        sys.exit(1)
    
    valid_extensions = ['.fa', '.faa', '.fasta']
    if not any(args.input_file.lower().endswith(ext) for ext in valid_extensions):
        print(f"Warning: Input file should have one of these extensions: {valid_extensions}")
    
    try:
        classifier = COVIDClassifier(args.model_dir, args.prodigal_path)
        
        # Process sequences if requested
        if args.force_convert or not args.no_auto_convert:
            try:
                processed_file, was_converted = classifier.process_input_sequences(
                    args.input_file, force_conversion=args.force_convert
                )
                if was_converted:
                    print(f"Using converted protein sequences from: {processed_file}")
                    input_file = processed_file
                else:
                    input_file = args.input_file
            except Exception as e:
                print(f"Warning: Sequence processing failed: {e}")
                print("Proceeding with original input file...")
                input_file = args.input_file
        else:
            input_file = args.input_file
        
        if args.run_all:
            # Run complete pipeline
            classifier.run_all_phases(input_file, args.output_dir)
        elif args.phase:
            classifier.predict(args.phase, input_file, args.output)
        else:
            classifier.interactive_mode(input_file)
            
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()