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

def test_min_max_mean():
    data = [1, 0.4, 5.6, 0.89]
    out_min, out_max, out_mean = l.min_max_mean(data)
    wanted_min = 0.4
    wanted_max = 5.6
    wanted_mean = 1.9725
    assert out_min == wanted_min
    assert out_max == wanted_max
    assert pytest.approx(out_mean, 0.001) == wanted_mean

def test_create_array():
    tested_output = l.create_array(1, 2.4, 0.5)
    wanted_output = np.array([1, 1.5, 2])
    np.testing.assert_array_equal(tested_output, wanted_output)

def test_create_arrayZ():
    residues = [] # pb ici car objet MDAnalysis
    list_resids = [1, 2]
    atom = 'C2'
    tested_output = l.create_arrayZ(residues, list_resids, atom, 1.0, 13.0)
    # supposons z_coord = 7 puis 10
    wanted_output = {1: [13., 12., 11., 10., 9., 8., 7., 6.], 2: [13., 12., 11., 10., 9.]}
    assert tested_output == wanted_output
