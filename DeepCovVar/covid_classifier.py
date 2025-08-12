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
from transformers import AutoTokenizer, AutoModelForSequenceClassification

from .features import FEATURE
from .neminer_utils import preprocess, preprocessdf

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
    def __init__(self, model_dir="models"):
        self.model_dir = Path(model_dir)
        self.feature_extractor = FEATURE()
        
        self.models_config = {
            1: {
                'file': 'p1_40_final_model_quantized.keras',
                'description': 'Virus sequences vs Others (Quantized)',
                'classes': ['Virus', 'Non-virus'],
                'feature_size': 2400,  # CKSAAP with gap=5: (5+1)*400 = 2400
                'type': 'keras'
            },
            2: {
                'file': 'p2_final_model_quantized.keras',
                'description': '(+)ssRNA (Class IV) vs Others (Quantized)',
                'classes': ['Others', '(+)ssRNA'],
                'feature_size': 2400,
                'type': 'keras'
            },
            3: {
                'file': 'p3_final_model_quantized.keras',
                'description': 'Further classification of group IV ssRNA(+) (Quantized)',
                'classes': ['Other ssRNA(+)', 'Coronavirus'],
                'feature_size': 2400,
                'type': 'keras'
            },
            4: {
                'file': 'p4_final_model_quantized.keras',
                'description': 'SARS-CoV-2 vs other coronaviruses (Quantized)',
                'classes': ['SARS-CoV-2', 'MERS', 'SARS', 'Others'],
                'feature_size': 2400,
                'type': 'keras'
            },
            5: {
                'file': 'p5_final_model_quantized.pt',
                'description': 'SARS-CoV-2 variant classification (Quantized)',
                'classes':  ['Omicron', 'Alpha', 'Delta', 'Others'], 
                'feature_size': None, 
                'type': 'pytorch_transformer'
            }
        }
        
        self.loaded_models = {}
        self.tokenizer = None
    
    def load_transformer_model(self, model_path):
        model_dir = Path(model_path)
        
        config_path = model_dir / "config.pt"
        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")
        
        config = torch.load(config_path, map_location='gpu')
        print(f"Loaded model config: {config}")
        
        
        base_model_name = config.get('base_model_name', 'facebook/esm2_t33_650M_UR50D')
        num_classes = config.get('num_classes', 4)
        
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
            
            # Try to load quantized state dict first, fall back to original if not found
            state_dict_path = model_dir / "model_state_dict_quantized.pt"
            if not state_dict_path.exists():
                state_dict_path = model_dir / "model_state_dict.pt"
                print("Quantized model not found, using original model")
            else:
                print("Loading quantized model weights...")
            
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
            model_path = self.model_dir / config['file']
            
            # If quantized model doesn't exist, fall back to original
            if not model_path.exists():
                original_file = config['file'].replace('_quantized.keras', '.keras')
                model_path = self.model_dir / original_file
                if model_path.exists():
                    print(f"Quantized model not found, using original: {original_file}")
                    config['description'] = config['description'].replace(' (Quantized)', '')
                else:
                    raise FileNotFoundError(f"Model file not found: {model_path}")
            else:
                print(f"Using quantized model: {config['file']}")
            
            print(f"Loading model for Phase {phase}: {config['description']}")
            
            if config['type'] == 'keras':
                model = tf.keras.models.load_model(str(model_path))
            elif config['type'] == 'pytorch':
                model = torch.load(str(model_path), map_location='gpu')
                model.eval()
        
        self.loaded_models[phase] = model
        return model
    
    def read_sequences(self, input_file):
        sequences = []
        seq_ids = []
        
        try:
            for record in SeqIO.parse(input_file, "fasta"):
                sequences.append(str(record.seq).upper())
                seq_ids.append(record.id)
            
            if not sequences:
                raise ValueError("No sequences found in the input file")
                
            print(f"Successfully read {len(sequences)} sequences from {input_file}")
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
    
    def predict(self, phase, input_file, output_file=None):
        model = self.load_model(phase)
        config = self.models_config[phase]
        
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
                # Binary classification
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
            input_ids, attention_mask = self.tokenize_sequences(sequences)
            
            with torch.no_grad():
                outputs = model(input_ids=input_ids, attention_mask=attention_mask)
                logits = outputs.logits if hasattr(outputs, 'logits') else outputs
                predictions = torch.softmax(logits, dim=1).numpy()
                
                predicted_classes = np.argmax(predictions, axis=1)
                confidence_scores = np.max(predictions, axis=1)
    
        results = []
        for i, (seq_id, pred_class, confidence) in enumerate(zip(seq_ids, predicted_classes, confidence_scores)):
            class_name = config['classes'][pred_class] if pred_class < len(config['classes']) else 'Unknown'
            results.append({
                'Sequence_ID': seq_id,
                'Predicted_Class': class_name,
                'Confidence': confidence,
                'Class_Index': pred_class
            })
        

        results_df = pd.DataFrame(results)
        

        print(f"\n{'='*60}")
        print(f"PHASE {phase} RESULTS: {config['description']}")
        print(f"{'='*60}")
        print(results_df.to_string(index=False))
        

        if output_file:
            results_df.to_csv(output_file, index=False)
            print(f"\nResults saved to: {output_file}")
        
        return results_df
    
    def interactive_mode(self, input_file):

        print("\n" + "="*60)
        print("SARSuite - INTERACTIVE MODE")
        print("="*60)
        
        print("\nAvailable phases:")
        for phase, config in self.models_config.items():
            print(f"  {phase}. {config['description']}")
        
        while True:
            try:
                phase = int(input(f"\nSelect phase (1-{len(self.models_config)}): "))
                if phase in self.models_config:
                    break
                else:
                    print(f"Invalid phase. Please select from 1-{len(self.models_config)}")
            except ValueError:
                print("Please enter a valid number")
        
        save_output = input("\nSave results to file? (y/n): ").lower().startswith('y')
        output_file = None
        if save_output:
            default_output = f"phase_{phase}_results.csv"
            output_file = input(f"Output filename [{default_output}]: ").strip()
            if not output_file:
                output_file = default_output
        
        return self.predict(phase, input_file, output_file)


def main():
    parser = argparse.ArgumentParser(
        description="SARSuite - COVID-19 Variant Classifier (Using Quantized Models)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog= """
                    Examples:
                    # Interactive mode
                    python covid_classifier.py input_sequences.fasta
                    
                    # Direct phase selection
                    python covid_classifier.py input_sequences.fasta --phase 1
                    
                    # With output file
                    python covid_classifier.py input_sequences.fasta --phase 1 --output results.csv
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
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found")
        sys.exit(1)
    
    valid_extensions = ['.fa', '.faa', '.fasta']
    if not any(args.input_file.lower().endswith(ext) for ext in valid_extensions):
        print(f"Warning: Input file should have one of these extensions: {valid_extensions}")
    
    try:
        classifier = COVIDClassifier(args.model_dir)
        
        if args.phase:
            classifier.predict(args.phase, args.input_file, args.output)
        else:
            classifier.interactive_mode(args.input_file)
            
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()