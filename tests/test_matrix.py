import pytest
import numpy as np

import sys
import os

# Add the parent rep to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))


from  bin import matrix as m

def test_initialize_matrix2D():
    tested_ouput = m.initialize_matrix2D(4, 4, np.nan)
    wanted_output = np.array([[np.nan, np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan, np.nan]])
    np.testing.assert_array_equal(tested_ouput, wanted_output)

def test_diff_Z():
    listZ = np.array([63.54, 62.54, 61.54, 60.54, 59.54, 58.54, 57.54, 56.54, 55.54, 54.54, 53.54, 52.54, 51.54, 50.54])
    tested_ouput = round(m.diff_Z(listZ, 52.18),2)
    wanted_output = -1.64
    assert(tested_ouput == wanted_output)

def test_find_X_Y():
    coord = np.array([44.79, 32.78, 52.18])
    listX = np.arange(-9.0, 95.0)
    listY =  np.arange(-11.0, 96.0)
    tested_ouputX, tested_ouputY = m.find_X_Y(coord, listX, listY)
    wanted_outputX, wanted_outputY = 53, 43
    assert(tested_ouputX == wanted_outputX)
    assert(tested_ouputY == wanted_outputY)

def test_check_edges():
    tested_ouput_T = m.check_edges(43, 5, 107)
    tested_ouput_F = m.check_edges(2, 5, 104)
    wanted_output_T = 43
    wanted_output_F = 5
    assert(tested_ouput_T == wanted_output_T)
    assert(tested_ouput_F ==  wanted_output_F)

def test_setDefects():
    tested_ouput_polar = m.setDefects('n', 0.0)
    polar_output = 1.0
    tested_ouput_apolar = m.setDefects('a', 0.0)
    apolar_output = 0.001
    assert(tested_ouput_polar == polar_output)
    assert(tested_ouput_apolar == apolar_output)

def test_my_fill_matrix():
    listX = np.arange(-1., 6., 1)
    listY = np.arange(-1., 6., 1)
    Matrix = np.full((len(listX), len(listY)), np.nan)
    tested_ouput = m.my_fill_matrix(Matrix, 1.34, 'a', [0, 4, 2], listX, listY, np.array([4, 3, 2, 1, 0]))
    wanted_output = np.array([[np.nan, np.nan, np.nan, np.nan, 0.5, 0.5, 0.5],
                              [np.nan, np.nan, np.nan, 0.5, 0.5, 0.5, 0.5],
                              [np.nan, np.nan, np.nan, np.nan, 0.5, 0.5, 0.5],
                              [np.nan, np.nan, np.nan, np.nan, np.nan, 0.5, np.nan],
                              [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                              [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                              [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]])
    np.testing.assert_array_equal(tested_ouput, wanted_output)

def test_fill_matrix():
    listX = np.arange(-1., 6., 1)
    listY = np.arange(-1., 6., 1)
    Matrix = np.full((len(listX), len(listY)), np.nan)
    tested_ouput = m.fill_matrix(Matrix, [0,4,2], listX, listY, np.array([4, 3, 2, 1, 0]), 1.34, 'all', 'a')
    wanted_output = np.array([[np.nan, np.nan, np.nan, 0.0, 0.003, 0.003, 0.003],
                              [np.nan, np.nan, np.nan, 0.001, 0.003, 0.005, 0.003],
                              [np.nan, np.nan, np.nan, 0.0, 0.003, 0.003, 0.003],
                              [np.nan, np.nan, np.nan, 0.0, 0.0, 0.001, 0.0],
                              [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                              [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
                              [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]])
    #wanted_output = np.array([[0.0, 3., 3., 3.],
    #                           [1., 3., 5., 3.],
    #                           [0.0, 3., 3., 3.],
    #                           [0.0, 0.0, 1., 0.0]])
    np.testing.assert_array_equal(tested_ouput, wanted_output)

def test_binarize_matrix_without0():
    Matrix_bin = np.full((4, 4), 0.)
    Matrix = np.array([[0.0, 0.003, 5.0, 0.5], [0.001, 0.003, 1.005, 0.003], [0.99, 0.003, 0.003, 0.003], [0.0, 0.0, 0.001, 0.0]])
    tested_ouput = Matrix_bin = m.binarize_matrix_without0(Matrix, Matrix_bin , -0.01, 0.99)
    wanted_output = np.array([[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 1.0, 0.0], [1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]])
    np.testing.assert_array_equal(tested_ouput, wanted_output)

def test_modify_matrix():
    list_lab_edge = np.array([1, 5])
    mat_lab = np.array([[1, 1, 1, 1, 0],
                [1, 3, 3, 0, 0],
                [1, 0, 0, 4, 5],
                [2, 2, 0, 4, 5],
                [2, 0, 5, 5, 5]])
    mat = np.array([[1, 1, 1, 1, 0],
                [1, 1, 1, 0, 0],
                [1, 0, 0, 1, 1],
                [1, 1, 0, 1, 1],
                [1, 0, 1, 1, 1]])
    tested_ouput = m.modify_matrix(mat_lab, mat, list_lab_edge)
    wanted_output = np.array([[0, 0, 0, 0, 0],
                [0, 1, 1, 0, 0],
                [0, 0, 0, 1, 0],
                [1, 1, 0, 1, 0],
                [1, 0, 0, 0, 0]])
    np.testing.assert_array_equal(tested_ouput, wanted_output)

def test_count_edge_area():
    cluster_edge = np.array([1, 5])
    area_cluster = {1: 6, 2 : 3, 3 : 2, 4 : 2, 5 : 5}
    tested_ouput = m.count_edge_area(area_cluster, cluster_edge)
    wanted_output = 11
    assert(tested_ouput == wanted_output)

def test_clean_NA_inside():
    list_lab_edge = np.array([1, 5])
    mat_lab = np.array([[1, 1, 1, 1, 6],
                [1, 3, 3, 0, 0],
                [1, 0, 0, 4, 5],
                [2, 2, 0, 4, 5],
                [2, 7, 5, 5, 5]])
    mat = np.array([[1, 1, 1, 1, np.nan],
                [1, 1, 1, 0, 0],
                [1, 0, 0, 1, 1],
                [1, 1, 0, 1, 1],
                [1, np.nan, 1, 1, 1]])
    total_edge = 11
    tested_ouput_mat_lab, tested_ouput_tot_edge, tested_ouput_clustPb = m.clean_NA_inside(mat_lab, list_lab_edge, mat, total_edge)
    output_mat_lab = np.array([[1, 1, 1, 1, 1],
                [1, 3, 3, 0, 0],
                [1, 0, 0, 4, 5],
                [2, 2, 0, 4, 5],
                [2, 1, 5, 5, 5]])
    output_tot_edge = 13
    output_clustPb = {6:1, 7:1}
    np.testing.assert_array_equal(tested_ouput_mat_lab, output_mat_lab)
    assert(tested_ouput_tot_edge == output_tot_edge)
    assert(tested_ouput_clustPb == output_clustPb)
