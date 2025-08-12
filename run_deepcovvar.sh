#!/bin/bash

# DeepCovVar: COVID-19 Variant Classification Tool
# Usage: ./run_deepcovvar.sh <fasta_file> <output_dir> <phase>

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to show usage
show_usage() {
    echo "DeepCovVar: COVID-19 Variant Classification Tool"
    echo "Usage: $0 <fasta_file> <output_dir> <phase>"
    echo ""
    echo "Arguments:"
    echo "  fasta_file  Input FASTA file containing sequences"
    echo "  output_dir  Output directory for results"
    echo "  phase       Classification phase (1-5)"
    echo ""
    echo "Phases:"
    echo "  1: Virus vs Non-virus classification"
    echo "  2: (+)ssRNA (Class IV) vs Others"
    echo "  3: Coronavirus vs Other ssRNA(+)"
    echo "  4: SARS-CoV-2 vs Other Coronaviruses"
    echo "  5: SARS-CoV-2 Variant Classification"
    echo ""
    echo "Examples:"
    echo "  $0 input.fasta results 1"
    echo "  $0 sequences.fa output 4"
    echo ""
    echo "For more options, use: python -m DeepCovVar --help"
}

# Check if help is requested
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_usage
    exit 0
fi

# Check number of arguments
if [[ $# -ne 3 ]]; then
    print_error "Incorrect number of arguments. Expected 3, got $#"
    show_usage
    exit 1
fi

FASTA_FILE="$1"
OUTPUT_DIR="$2"
PHASE="$3"

# Validate arguments
print_info "Validating arguments..."

# Check if FASTA file exists
if [[ ! -f "$FASTA_FILE" ]]; then
    print_error "FASTA file not found: $FASTA_FILE"
    exit 1
fi

# Check if phase is valid
if [[ ! "$PHASE" =~ ^[1-5]$ ]]; then
    print_error "Invalid phase: $PHASE. Must be 1, 2, 3, 4, or 5"
    exit 1
fi

# Create output directory if it doesn't exist
if [[ ! -d "$OUTPUT_DIR" ]]; then
    print_info "Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

print_success "Arguments validated successfully"

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEEPVAR_DIR="$SCRIPT_DIR/DeepCovVar"

# Check if DeepCovVar is installed
if [[ ! -d "$DEEPVAR_DIR" ]]; then
    print_error "DeepCovVar package not found in $DEEPVAR_DIR"
    print_info "Please ensure the package is properly installed"
    exit 1
fi

# Check if models directory exists
MODELS_DIR="$DEEPVAR_DIR/models"
if [[ ! -d "$MODELS_DIR" ]]; then
    print_error "Models directory not found: $MODELS_DIR"
    exit 1
fi

# Check for quantized models
QUANTIZED_MODELS=$(find "$MODELS_DIR" -name "*quantized.keras" | wc -l)
if [[ $QUANTIZED_MODELS -eq 0 ]]; then
    print_warning "No quantized models found in $MODELS_DIR"
    print_info "This may cause issues with model loading"
else
    print_info "Found $QUANTIZED_MODELS quantized models"
fi

# Run DeepCovVar
print_info "Starting DeepCovVar classification..."
print_info "Input file: $FASTA_FILE"
print_info "Output directory: $OUTPUT_DIR"
print_info "Phase: $PHASE"

# Change to DeepCovVar directory and run
cd "$DEEPVAR_DIR"

# Run the classification
if python -m DeepCovVar -f "$FASTA_FILE" -o "$OUTPUT_DIR" -p "$PHASE"; then
    print_success "Classification completed successfully!"
    print_info "Results saved to: $OUTPUT_DIR/phase_${PHASE}_results.csv"
else
    print_error "Classification failed!"
    exit 1
fi

print_info "DeepCovVar run completed"
print_info "Check the output directory for results: $OUTPUT_DIR"
