import pytest
import numpy as np

import sys
import os

# Add the parent rep to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))


from  bin import listes as l

def test_get_glyc_lipids():
    glycerol = {'DOP': 'C2', 'DOE': 'C2', 'DPP': 'C2'}
    lipids = ['DOP', 'DPP']
    tested_output = l.get_glyc_lipids(lipids, glycerol)
    wanted_output = 'C2'
    assert tested_output == wanted_output
