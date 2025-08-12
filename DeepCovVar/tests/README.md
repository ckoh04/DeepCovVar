# DeepCovVar Test Suite

This directory contains comprehensive tests for the DeepCovVar COVID classifier with nucleotide sequence processing capabilities.

## Test Files

### Test Scripts
1. **`test_nucleotide_detection.py`**
   - Tests the sequence type detection functionality
   - File-based sequence type detection
   - Individual sequence detection
   - Nucleotide vs protein classification accuracy

2. **`test_pipeline_integration.py`**
   - Tests the complete pipeline integration
   - Sequence processor functionality
   - COVID classifier integration
   - Pipeline method availability
   - Prodigal availability check

3. **`test_deepcovvar.py`**
   - Original DeepCovVar test script
   - Tests basic functionality and model loading

4. **`run_tests.py`**
   - Test runner script that executes all available tests

### Test Data Files
1. **`test_nucleotide_sequences.fasta`**
   - Contains dummy nucleotide sequences for testing
   - **SARS_CoV_2_spike_gene**: Coronavirus spike protein gene
   - **MERS_spike_gene**: MERS coronavirus spike protein gene
   - **Influenza_H1N1_HA_gene**: Influenza hemagglutinin gene
   - **Human_ribosomal_RNA**: Human ribosomal RNA sequence

2. **`p3_test.fasta`**
   - Test sequences for Phase 3 classification
   - Coronavirus vs Other ssRNA(+) classification

3. **`p5_test.fasta`**
   - Test sequences for Phase 5 classification
   - SARS-CoV-2 variant classification

4. **`test.fasta`**
   - General test sequences for various phases

5. **`test_sequences.fasta`**
   - Basic test sequences for functionality testing

## Running Tests

### Run All Tests
```bash
cd DeepCovVar/tests
python run_tests.py
```

### Run Individual Tests
```bash
# Test nucleotide detection only
python test_nucleotide_detection.py

# Test pipeline integration only
python test_pipeline_integration.py
```

### Run Tests from Project Root
```bash
cd CovProj/DeepCovVar
python DeepCovVar/tests/run_tests.py
```

## Test Requirements

### Required Dependencies
- Python 3.7+
- BioPython
- All DeepCovVar dependencies

### Optional Dependencies
- **Prodigal**: For nucleotide-to-protein conversion testing
  ```bash
  conda install -c bioconda prodigal
  ```

## Expected Test Results

### Without Prodigal
- ✅ Sequence detection tests will pass
- ✅ Pipeline integration tests will pass
- Warning: Prodigal availability test will show warning (normal)

### With Prodigal
- ✅ All tests should pass
- ✅ Nucleotide conversion functionality will be tested
- ✅ Complete pipeline can be tested end-to-end

## Test Output

Tests provide detailed output including:
- Success indicators
- Failure indicators
- Warning messages
- Detailed error information
- Summary of passed/failed tests

## Troubleshooting

### Common Issues

1. **Import Errors**: Ensure you're running from the correct directory
2. **File Not Found**: Check that test files exist in the tests directory
3. **Prodigal Not Found**: Install Prodigal or ignore the warning (conversion won't work)

### Debug Mode

For detailed debugging, run individual test scripts directly:
```bash
python -v test_nucleotide_detection.py
```

## Adding New Tests

To add new tests:

1. Create test file in this directory
2. Follow the naming convention: `test_*.py`
3. Add the test to `run_tests.py` if it should be part of the test suite
4. Ensure tests return `True` for success, `False` for failure

## Test Coverage

Current tests cover:
- ✅ Sequence type detection
- ✅ Sequence processor integration
- ✅ COVID classifier integration
- ✅ Pipeline method availability
- ✅ Prodigal integration (when available)
- ✅ Basic DeepCovVar functionality
- ✅ Phase-specific test data (P3, P5)
- ✅ General test sequences

Test data includes:
- ✅ Nucleotide sequences for conversion testing
- ✅ Protein sequences for classification testing
- ✅ Phase-specific sequences for targeted testing
- ✅ Mixed sequence types for comprehensive testing

Future tests could cover:
- End-to-end pipeline execution
- Model loading and prediction
- Error handling scenarios
- Performance benchmarks
- Cross-phase consistency validation
