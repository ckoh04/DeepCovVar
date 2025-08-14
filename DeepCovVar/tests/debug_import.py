#!/usr/bin/env python3
import sys
import os
from pathlib import Path

print("Current working directory:", os.getcwd())
print("Script location:", __file__)
print("Script parent:", Path(__file__).parent)
print("Script parent parent:", Path(__file__).parent.parent)
print("Script parent parent absolute:", Path(__file__).parent.parent.absolute())

# Add the parent directory to the path
parent_dir = str(Path(__file__).parent.parent.absolute())
print("Parent dir to add:", parent_dir)

sys.path.insert(0, parent_dir)
print("Python path after insertion:", sys.path[:3])

# Try to list files in parent directory
try:
    parent_files = os.listdir(parent_dir)
    print("Files in parent directory:", [f for f in parent_files if f.endswith('.py')])
except Exception as e:
    print("Error listing parent directory:", e)

# Try to import sequence_converter
try:
    import sequence_converter
    print("✓ sequence_converter imported successfully")
except ImportError as e:
    print("✗ sequence_converter import failed:", e)
