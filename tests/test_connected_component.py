import sys
import os
import numpy as np
import MDAnalysis as mda

# Add the parent rep to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))


from bin import connected_component as cc

def test_intersect():
    array_1 = np.array([1, 2, 3, 4, 5])
    array_2 = np.array([2, 4, 5, 5, 6])

    tested_out = cc.intersect(array_1, array_2)
    wanted_out = np.array([2, 4, 5])
    assert np.allclose(tested_out, wanted_out)
