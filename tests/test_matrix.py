import pytest
import numpy as np

import sys
import os

# Add the parent rep to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))


from  bin import matrix as m

def test_initialize_matrix2D():
    input = m.initialize_matrix2D(4, 4, np.nan)
    wanted_output = np.array([[np.nan, np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan, np.nan]])
    np.testing.assert_array_equal(input, wanted_output)


def test_fill_matrix():
    Matrix = np.full((4, 4), np.nan)
    listX = np.arange(-1., 6., 1)
    listY = np.arange(-1., 6., 1)
    input = m.fill_matrix(Matrix, [0,4,2], listX, listY, np.array([4, 3, 2, 1, 0]), 1.34, 'all', 'a')
    wanted_output = np.array([[0.0, 0.003, 0.003, 0.003], [0.001, 0.003, 0.005, 0.003], [0.0, 0.003, 0.003, 0.003], [0.0, 0.0, 0.001, 0.0]])
    np.testing.assert_array_equal(input, wanted_output)

def test_binarize_matrix_without0():
    Matrix_bin = np.full((4, 4), 0.)
    Matrix = np.array([[0.0, 0.003, 5.0, 0.5], [0.001, 0.003, 1.005, 0.003], [0.99, 0.003, 0.003, 0.003], [0.0, 0.0, 0.001, 0.0]])
    input = Matrix_bin = m.binarize_matrix_without0(Matrix, Matrix_bin , -0.01, 0.99)
    wanted_output = np.array([[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 1.0, 0.0], [1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]])
    np.testing.assert_array_equal(input, wanted_output)
