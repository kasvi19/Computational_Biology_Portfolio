import sys
import os
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

# Debugging the path
print("sys.path:", sys.path)

import pytest
from melting_temp import calculate_gc_content, calculate_at_content, calculating_melting_temp
SEQUENCE = "ATGC"

def test_gc_content():
    assert calculate_gc_content(SEQUENCE) == 50.0, "GC content calculation failed."

def test_at_content():
    assert calculate_at_content(SEQUENCE) == 50.0, "AT content calculation failed."

def test_melting_temperature():
    assert calculating_melting_temp(SEQUENCE) == 12.0, "Melting temperature calculation failed."

def test_empty_sequence():
    assert calculate_gc_content('') == 0.0, "GC content for empty sequence failed."
    assert calculate_at_content('') == 0.0, "AT content for empty sequence failed."
    assert calculating_melting_temp('') == 0.0, "Melting temperature for empty sequence failed."

        
        
