import pytest
import io
import numpy as np

import sys
import os

# Add the parent rep to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))


from bin import pdb as pdb

def test_write_a_pdb_line():
    # Create fake file
    fake_file = io.StringIO()
    # Write
    pdb.write_a_pdb_line(fake_file, 1, 'EDG', 1, [-8.0, -9.0, 7.310], -1.0)

    # Read the fake file
    fake_file.seek(0)  # Seek the start of the file
    content = fake_file.read()

    wanted_line = f"ATOM      1   H1 EDG     1      -8.000  -9.000   7.310  1.00 -1.00\n"
    assert content == wanted_line

def  test_outputPDB_Total_matrix():
    arrayX = [-1., 0.]
    arrayY = [-3., -2.]
    z_extr = 4.430
    Matrix_final = np.array([[1, np.nan],
                             [1, 1]])

    pdb.outputPDB_Total_matrix('output', 1, arrayX, arrayY, z_extr, Matrix_final)

    with open("output.pdb", 'r') as file_test:
        content = file_test.readlines()
    wanted_lines = ["MODEL        1\n",
                    "ATOM      1   H1 MAT     1      -1.000  -3.000   4.430  1.00  1.00\n",
                    "ATOM      1   H1 EDG     1      -1.000  -2.000   4.430  1.00 -1.00\n",
                    "ATOM      2   H1 MAT     2       0.000  -3.000   4.430  1.00  1.00\n",
                    "ATOM      2   H1 MAT     2       0.000  -2.000   4.430  1.00  1.00\n",
                    "ENDMDL\n"]
    assert content == wanted_lines
