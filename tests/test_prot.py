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
    
    tested_out = prot.find_edges(mat)
    wanted_out = {2:np.array([np.array([0, 0]), np.array([0, 1]), np.array([1, 0])]),
                3:np.array([np.array([0, 3])]),
                4:np.array([np.array([2, 2]), np.array([2, 3])])}

    # Compare each element
    for label in wanted_out:
        assert np.allclose(tested_out[label], wanted_out[label])

def test_find_shortest_sqdist():
    coor_A = np.array([np.array([0, 1]), np.array([2, 5]), np.array([3, 4])])
    coor_B = np.array([np.array([6, 7]), np.array([5, 1]), np.array([9, 0]), np.array([4, 8])])

    tested_out = prot.find_shortest_sqdist(coor_A, coor_B)
    wanted_out = 13
    assert tested_out == wanted_out

def test_assign_dist_group():
    dict_edge_defect = {1: np.array([np.array([6, 7]), np.array([5, 1]), np.array([9, 0]), np.array([4, 8])]),
                        3: np.array([np.array([8, 3]), np.array([7, 4]), np.array([6, 1]), np.array([8, 8])]),
                        4: np.array([np.array([0, 4]), np.array([6, 8]), np.array([4, 1]), np.array([5, 9])])}
    prot_edge = np.array([np.array([0, 1]), np.array([2, 5]), np.array([3, 4])])
    dist_threshold = 3

    tested_out = prot.assign_dist_group(prot_edge, dict_edge_defect, dist_threshold)
    wanted_out = {1: 'far', 3: 'far', 4: 'close'}

    # Compare each element
    for label in wanted_out:
        assert tested_out[label] == wanted_out[label]
