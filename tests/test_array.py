import pytest
import numpy as np
import MDAnalysis as mda

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
    n_atoms = 2
    u = mda.Universe.empty(n_atoms=n_atoms,
                           n_residues=n_atoms,
                           atom_resindex=list(range(n_atoms)),
                           trajectory=True)

    # Add atom names
    atom_name = "C2"
    u.add_TopologyAttr("name", [atom_name] * n_atoms)

    # Add residue IDs
    list_resids = [1, 2]
    u.add_TopologyAttr("resid", list_resids)

    # Add coordinates
    coords = np.zeros((n_atoms, 3))
    # Specify the z values
    coords[:, 2] = [7.0, 10.0]
    u.atoms.positions = coords

    tested_output = l.create_arrayZ(u.residues, list_resids, atom_name, 1.0, 13.0, up=True)
    wanted_output = {1: [13., 12., 11., 10., 9., 8., 7., 6.], 2: [13., 12., 11., 10., 9.]}
    # Compare each element
    for resid in wanted_output:
        assert np.allclose(tested_output[resid], wanted_output[resid])
