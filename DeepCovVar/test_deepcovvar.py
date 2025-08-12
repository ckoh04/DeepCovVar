#!/usr/bin/env python3
"""
Test script for DeepCovVar package
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

def test_imports():
    """Test that all modules can be imported correctly."""
    try:
        from DeepCovVar.covid_classifier import COVIDClassifier
        print("✓ COVIDClassifier imported successfully")
        
        from DeepCovVar.features import FEATURE
        print("✓ FEATURE imported successfully")
        
        from DeepCovVar.neminer_utils import preprocess, preprocessdf
        print("✓ neminer_utils imported successfully")
        
        from DeepCovVar import feature_data
        print("✓ feature_data imported successfully")
        
        return True
    except ImportError as e:
        print(f"✗ Import error: {e}")
        return False

def test_classifier_initialization():
    """Test that the COVIDClassifier can be initialized."""
    try:
        from DeepCovVar.covid_classifier import COVIDClassifier
        
        # Test with default model directory
        classifier = COVIDClassifier()
        print("✓ COVIDClassifier initialized with default model directory")
        
        # Test with custom model directory
        custom_dir = Path(__file__).parent / "models"
        classifier = COVIDClassifier(model_dir=str(custom_dir))
        print("✓ COVIDClassifier initialized with custom model directory")
        
        return True
    except Exception as e:
        print(f"✗ Initialization error: {e}")
        return False

def test_model_files():
    """Test that model files exist and are accessible."""
    try:
        model_dir = Path(__file__).parent / "models"
        
        if not model_dir.exists():
            print(f"✗ Model directory not found: {model_dir}")
            return False
        
        # Check for quantized models
        quantized_models = list(model_dir.glob("*quantized.keras"))
        if not quantized_models:
            print("✗ No quantized models found")
            return False
        
        print(f"✓ Found {len(quantized_models)} quantized models:")
        for model in quantized_models:
            print(f"  - {model.name}")
        
        return True
    except Exception as e:
        print(f"✗ Model file test error: {e}")
        return False

def main():
    """Run all tests."""
    print("DeepCovVar Package Tests")
    print("=" * 40)
    
    tests = [
        ("Import Tests", test_imports),
        ("Classifier Initialization", test_classifier_initialization),
        ("Model Files", test_model_files),
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        print(f"\n{test_name}:")
        print("-" * len(test_name))
        
        if test_func():
            passed += 1
            print(f"✓ {test_name} PASSED")
        else:
            print(f"✗ {test_name} FAILED")
    
    print(f"\nTest Summary: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! DeepCovVar is ready to use.")
        return 0
    else:
        print("❌ Some tests failed. Please check the errors above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
