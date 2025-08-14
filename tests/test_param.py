import pytest
import numpy as np
import sys
import os

# Add the parent rep to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from bin import param as p

def test_file_present():
    # Test if file present
    filename_OK = "src/PackMem_prot.py"
    p.file_present(str(filename_OK))

    # Test if file absent
    filename_NO = "piou.txt"
    with pytest.raises(FileNotFoundError):
        p.file_present(str(filename_NO))