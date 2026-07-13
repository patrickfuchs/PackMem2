import pytest
import numpy as np
import MDAnalysis as mda

import packmem2.core.arrays as a

def test_get_glyc_lipids():
    glycerol = {'DOP': 'C2', 'DOE': 'C2', 'DPP': 'C2'}
    lipids = ['DOP', 'DPP']
    tested_output = a.get_glyc_lipids(lipids, glycerol)
    wanted_output = 'C2'
    assert tested_output == wanted_output

def test_min_max_mean():
    data = [1, 0.4, 5.6, 0.89]
    out_min, out_max, out_mean = a.min_max_mean(data)
    wanted_min = 0.4
    wanted_max = 5.6
    wanted_mean = 1.9725
    assert out_min == wanted_min
    assert out_max == wanted_max
    assert pytest.approx(out_mean, 0.001) == wanted_mean

def test_create_array():
    tested_output = a.create_array(1, 2.4, 0.5)
    wanted_output = np.array([1, 1.5, 2])
    np.testing.assert_array_equal(tested_output, wanted_output)

def test_create_arrayZ():
    dict_glyc ={"DMPC": "C2", "POPC" : "C2", "ARG": "N"}
    n_atoms = 2
    u = mda.Universe.empty(n_atoms=n_atoms,
                           n_residues=n_atoms,
                           atom_resindex=list(range(n_atoms)),
                           trajectory=True)

    # Add atom names
    atom_name = "C2"
    u.add_TopologyAttr("name", [atom_name] * n_atoms)
    
    # Add residue names
    res_name = "DMPC"
    u.add_TopologyAttr("resname", [res_name] * n_atoms)

    # Add residue IDs
    list_resids = [1, 2]
    u.add_TopologyAttr("resid", list_resids)

    # Add coordinates
    coords = np.zeros((n_atoms, 3))
    # Specify the z values
    coords[:, 2] = [7.0, 10.0]
    u.atoms.positions = coords

    # Test Upper leaflet
    tested_output = a.create_arrayZ(u.residues, list_resids, dict_glyc, 1.0, 13.0, up=True)
    wanted_output = {1: np.array([13., 12., 11., 10., 9., 8., 7., 6.]),
                     2: np.array([13., 12., 11., 10., 9.])}
    # Compare each element
    for resid in wanted_output:
        assert np.allclose(tested_output[resid], wanted_output[resid])
    
    # Test Lower leaflet
    # Add coordinates
    coords = np.zeros((n_atoms, 3))
    # Specify the z values
    coords[:, 2] = [7.34002, 10.7]
    u.atoms.positions = coords

    tested_output = a.create_arrayZ(u.residues, list_resids, dict_glyc, 1.0, 3.24, up=False)
    wanted_output = {1: np.array([2.34, 3.34, 4.34, 5.34, 6.34, 7.34, 8.34]),
                     2: np.array([2.7, 3.7, 4.7, 5.7, 6.7, 7.7, 8.7, 9.7,
                                  10.7, 11.7])}
    # Compare each element
    for resid in wanted_output:
        assert np.allclose(tested_output[resid], wanted_output[resid])
