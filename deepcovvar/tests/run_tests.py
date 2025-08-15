#!/usr/bin/env python3
import sys
import subprocess
from pathlib import Path

def run_test_script(script_path):
    try:
        result = subprocess.run([sys.executable, str(script_path)], 
                              capture_output=True, text=True, timeout=60)
        
        print(f"\n{'='*60}")
        print(f"Running: {script_path.name}")
        print(f"{'='*60}")
        
        if result.stdout:
            print(result.stdout)
        
        if result.stderr:
            print("STDERR:")
            print(result.stderr)
        
        success = result.returncode == 0
        status = "PASSED" if success else "FAILED"
        print(f"\nResult: {status}")
        
        return success
        
    except subprocess.TimeoutExpired:
        print(f"Test {script_path.name} timed out after 60 seconds")
        return False
    except Exception as e:
        print(f"Error running {script_path.name}: {e}")
        return False

def main():
    print("DeepCovVar Test Suite")
    print("=" * 50)
    
    # Find all test scripts
    tests_dir = Path(__file__).parent
    test_scripts = [
        tests_dir / "test_nucleotide_detection.py",
        tests_dir / "test_pipeline_integration.py",
        tests_dir / "test_deepcovvar.py"
    ]
    
    # Filter to only existing scripts
    available_tests = [script for script in test_scripts if script.exists()]
    
    if not available_tests:
        print("No test scripts found!")
        return 1
    
    print(f"Found {len(available_tests)} test script(s):")
    for script in available_tests:
        print(f"  - {script.name}")
    
    print(f"\nRunning tests...")
    
    # Run all tests
    results = []
    for script in available_tests:
        success = run_test_script(script)
        results.append((script.name, success))
    
    # Summary
    print(f"\n{'='*60}")
    print("TEST SUMMARY")
    print(f"{'='*60}")
    
    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    for script_name, success in results:
        status = "PASSED" if success else "FAILED"
        print(f"{script_name}: {status}")
    
    print(f"\nOverall: {passed}/{total} test suites passed")
    
    if passed == total:
        print("\nAll tests passed!")
        return 0
    else:
        print(f"\n{total - passed} test suite(s) failed.")
        return 1

if __name__ == "__main__":
    exit(main())
