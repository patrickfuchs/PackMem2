import sys
import os
import numpy as np
import MDAnalysis as mda

# Add the parent rep to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))


from bin import protdist as prot

def test_find_protein():
    n_atoms = 2
    u = mda.Universe.empty(n_atoms=n_atoms,
                           n_residues=n_atoms,
                           atom_resindex=list(range(n_atoms)),
                           trajectory=True)

    # Add residue IDs
    list_resids = [1, 2]
    u.add_TopologyAttr("resid", list_resids)

    # Add coordinates
    coords = np.array([[2.140001, 0.930004, 0.370003],
                [0.360003, 1.550003, 1.860004]])
    u.atoms.positions = coords

    matrix = np.zeros((3,3))
    defect = 'up'
    protein = u.atoms
    arrayX = np.array([0., 1., 2.])
    arrayY = np.array([0., 1., 2.])
    zmean =  1.0

    tested_out = prot.find_protein(matrix, defect, protein, arrayX, arrayY, zmean)
    wanted_out = np.array([[0, 1, 0],
                            [0, 0, 0],
                            [0, 0, 0]])
    np.testing.assert_array_equal(tested_out, wanted_out)

def test_find_edges():
    mat = np.array([[2, 2, 0, 3],
                   [2, 0, 0, 0],
                   [0, 0, 4, 4]])
    
    expected_out = prot.find_edges(mat)
    wanted_out = {2:[[0, 0], [0, 1], [1, 0]], 3:[[0, 3]], 4:[[2, 2], [2, 3]]}
    assert expected_out == wanted_out

