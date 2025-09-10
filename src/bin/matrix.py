#-*- coding: utf-8 -*-
"""Funtions for matrix."""
# R. Gautier A. Bacle 2015
# M. Zygadlo 2025

import math
import numpy as np

# matrix size  (square size = 1A)
SIZE=1.0
# calculate limit size (half diagonal) 3D
SIZE_SIDE=0.5*(math.sqrt((math.pow(SIZE, 2))* 3.))


def initialize_matrix2D(val1, val2, default):
    """
    Initialise a 2D matrix with a default value.
    
    --------------------
    INPUT
    val1: int
        Dimension in X/Y for the array
    val2: int
        Dimension in X/Y for the array
    default : int / float / string
        Default value in the matrix cells
        
    --------------------
    OUTPUT
    numpy matrix
        Matrix initialised with default value
    """
    matrix = np.full((val1, val2), default)
    return matrix

def diff_Z(arrayZ, coordZ):
    """
    Compute the difference between the last value in listZ and the Z for one atom.

    --------------------
    INPUT
    arrayZ: numpy array
        List from zmax+1 to z_C2_coord-1 by step of 1.0  OR
        List from z_C2_coord+1 to zmin-1 by step of 1.0
    coordZ: float
        The z coordinate for one atom

    --------------------
    OUPUT
    float
        The difference between the last value in listZ and the Z for one atom
    """
    return arrayZ[-1] - coordZ

def find_X_Y(coord, matX, matY):
    """
    Determine from x,y coordinates the index in matX and matY.
    
    Has been said to determine the central value in matrix

    --------------------
    INPUT
    coord: numpy array
        Contains the coordinates x,y,z of an atom
    matX: numpy array
        Contains integer from xmin-1 to xmax+1 by step of 1.0
    matY: numpy array
        Contains integer from ymin-1 to ymax+1 by step of 1.0
    
    --------------------
    OUTPUT
    tuple of int
        the index of the coord x and y of the atom in matX and matY
    """
    # Convert to integer
    tmpX = int(coord[0])
    tmpY = int(coord[1])
    # Find the index of the value given (x_atom or y_atom)
    iX = np.where(matX == tmpX)[0][0]
    iY = np.where(matY == tmpY)[0][0]
    return iX,iY

def check_edges(val, val_lim1, val_lim2):
    """
    Check if the cell index is lower than val_lim1 or higher than val_lim2. Changes the value of the cell index if so.

    --------------------
    INPUT
    val: int
        Index of the coordinate of the atom in listX or listY
    val_lim1: int
        The number of cell to work around the atom
    val_lim2: int
        The length of listX or listY

    --------------------
    OUTPUT
    int
        The cell index that represents the center of the search
    """
    if val < val_lim1:
        return val_lim1
    elif val >= val_lim2-val_lim1:
        return val_lim2-(val_lim1+1)
    else:
        return val

def dist_oneAxis(coord1, coord2):
    """
    Compute the distance for one axis without square root.

    --------------------
    INPUT
    coord1: float
        The position on one axis
    coord2:
        The position on one axis
    --------------------
    OUTPUT
    float
        The distance between the two input points
    """
    return (coord1 - coord2)**2

def dist(coord1, coord2):
    """
    Compute the euclidean distance in x,y,z without square root.

    --------------------
    INPUT
    coord1: numpy array
        The position in x,y,z
    coord2: numpy array
        The position in x,y,z
    --------------------
    OUTPUT
    float
        The distance between the two input points
    """
    return ((coord1[0]-coord2[0])**2 +
           (coord1[1]-coord2[1])**2 +
           (coord1[2]-coord2[2])**2)

def setDefects(type_aliphatic, val_mat):
    """
    Fill matrix cell depending on the defect type (Deep, Shallow).

    --------------------
    INPUT
    type_aliphatic: string
        The nature of the atom. a : aliphatique / n : polar
    val_mat: float
        The value in the matrix cell corresponding to the atom

    --------------------
    OUTPUT
    float
        The new value of the matrix cell
    """
    if type_aliphatic == "a":
        val_mat += 0.001
    else:
        val_mat += 1.
    return val_mat

def fill_matrix(mat, coordtmp, arrayX, arrayY, arrayZ,
                radius_res, FlagPDtype, aliph_atoms):
    """
    Fill the matrix for each atom.

    ---------------
    INPUT
    mat: numpy array 2D
        Matrix to be filled where there are atoms
    coordtmp: numpy array
        Contains the coordinates x,y,z of the atom
    arrayX: numpy array
        Array from xmin-1 to xmax+1 by step of 1.0
    arrayY: numpy array
        Array from ymin-1 to ymax+1 by step of 1.0
    arrayZ: numpy array
        Array from zmax+1 to z_C2_coord-1 by step of 1.0  OR
        Array from z_C2_coord+1 to zmin-1 by step of 1.0
    radius_res: float
        the radius of the atom type for the residue
    FlagPDtype: str
        The type of defect to analyse. all / deep / shallow
    aliph_atoms: string
        The nature of the atom. a : aliphatique / n : polar

    -----------------
    OUTPUT
    numpy array 2D
        Contains a value 0 < a < 1 for aliphatic atom
        > 1 if polar OR deep
        Defects = 0
    """
    # Number of cells to work around
    v=5
    # Find the index of x_atom and y_atom in listX and listY
    # this corresponds to the location of the atom in the matrix
    iX,iY = find_X_Y(coordtmp, arrayX, arrayY)
    # Change the indexes (cell index in matrix) if < v or > len(matrix)-5
    iX = check_edges(iX, v, len(arrayX))
    iY = check_edges(iY, v, len(arrayY))
    # Select the cells to work in at i+-v
    listXM = arrayX[iX - v : iX + (v + 1)]
    listYM = arrayY[iY - v : iY + (v + 1)]
    # Limit distance in a radius
    dist_lim = (SIZE + radius_res) ** 2
    # Limit distance in diagonal
    dist_meet = (SIZE_SIDE + radius_res) ** 2

    # Select valid positions to search in the radius of the distance
    # from the position of the atom to the slice in Z / X / Y
    validZ = [z for z in arrayZ if dist_oneAxis(coordtmp[2], z) <= dist_lim]
    validX = [(ix, x) for ix, x in enumerate(listXM) if dist_oneAxis(coordtmp[0], x) <= dist_lim]
    validY = [(iy, y) for iy, y in enumerate(listYM) if dist_oneAxis(coordtmp[1], y) <= dist_lim]
    # Loop on the different z positions of the upper OR lower leaflet
    for sliceZ in validZ:
        # Loop on the different cells in x dimension
        for indX, sliceX in validX:
            # Loop on the different cells in y dimension
            for indY, sliceY in validY:
                X = indX + (iX - v)
                Y = indY + (iY - v)
                # Create a list with the coordinates x,y,z of a cell
                # in a 5 cell radius of the atom
                coordCenter = np.array([sliceX, sliceY, sliceZ], dtype=float)
                # Compute the distance between the atom and this position of the matrix
                distance = dist(coordCenter, coordtmp)
                # If the matrix cell was empty, put 0.0
                if np.isnan(mat[X, Y]):
                    mat[X, Y] = 0.
                # Get the value in this cell of the matrix
                if distance <= dist_meet:
                    # If shallow or all defect
                    if FlagPDtype == "shallow" or FlagPDtype == "all":
                        mat[X, Y] = setDefects(aliph_atoms, mat[X, Y])
                    # If deep defect
                    elif FlagPDtype == "deep":
                        mat[X, Y] += 1.
    return mat

def binarize_matrix_without0(mat, mat_ini, val1=0, val2=0.99):
    """
    Binarise the presence / absence of aliphatic atom (and packing defects) in matrix.

    --------------------
    INPUT
    mat: numpy matrix
        Contains where there are atoms in the simulation box and their type
    mat_ini: numpy matrix
        Matrix to be binarised
    val1: float
        lower limit value
    val2: float
        higher limit value

    --------------------
    OUTPUT
    numpy matrix
        Contains the position of the aliphatic / polar atoms (+ packing defects) (1)
        in the simulation box    
    """
    # Get the index of the polar atoms / deep defects if val2 = 0.99
    # Get the index of the apolar atoms if val 2 = 0.001
    index = np.argwhere((mat >= val2) | (mat <= val1))
    mat_ini[index[:,0], index[:,1]] = 1. 
    return mat_ini

def count_edge_area(area_defects, edge_labels):
    """
    Count the total area of the clusters on the edge.

    --------------------
    INPUT
    area_defects: dictionary
        Contains the labels of each defects and their areas
    edge_labels: list
        Contains the labels of the defects on the edge
    
    --------------------
    OUTPUT
    int
        The sum of the defect areas on the edge
    """
    total_edge_area = 0
    for key in edge_labels:
        total_edge_area += area_defects[key]
    return  total_edge_area

def modify_matrix(mat_labels, mat, edge_labels):
    """
    Modify the binary matrix to take account edges (determined by all packing defects).
    
    The edges are put to 0.0

    --------------------
    INPUT
    mat_labels: numpy matrix
        Contains the labels of the packing defects + aliphatic atoms
    mat: numpy matrix
        Contains the position of the aliphatic atoms (0)
    edge_labels: list
        Contains the labels of the clusters on the edges
    
    --------------------
    OUTPUT
    matrix
        Contains mat where the clusters on the edge were put to 0.0
    """
    # Get the index in mat1 where the labels are in listval
    mask = np.isin(mat_labels, edge_labels)
    index = np.argwhere(mask)
    mat[index[:,0], index[:,1]] = 0. 
    return mat

def clean_NA_inside(mat_labels, edge_labels, mat_ini, total_edge):
    """
    Clean up the label matrix, the total edge area if there were NaN that weren't changed in the first fill matrix.
    
    Then gives a dictionary of the wrongly labeled defects.

    --------------------
    INPUT
    mat_labels : numpy matrix
        Contains the defects' labels
    edge_labels : list
        Contains the labels on each edge of the matrix
    mat_ini : numpy matrix
        Contains the positions of the polar atoms (int)
        aliphatic atoms (0.0 < float < 0.99)
        and packing defects (0.0)
    total_edge : int
        The area taken by the packing defects on the edge of the matrix

    --------------------
    OUTPUT
    numpy matrix
        Contains the labels without nan
    int
        The updated area taken by the packing defects on the edge of the matrix
    dictionnary
        Contains the labels and their area concerned by the NaN problem
    """
    # Select the NaN that are inside the memb and not registered as edges
    index_nan_inside = np.argwhere((np.isnan(mat_ini)) & (~np.isin(mat_labels, edge_labels)))
    labels_Pb = mat_labels[index_nan_inside[:,0], index_nan_inside[:,1]]
    # Get the unique labels and count their occurrence => dict
    unique, counts = np.unique(labels_Pb, return_counts=True)
    clustPb = dict(zip(unique, counts))
    # Correct the nan by giving it the first label of clust_edge
    mat_labels[index_nan_inside[:,0], index_nan_inside[:,1]] = edge_labels[0]
    # Correcy the total edge area
    total_edge += len(labels_Pb)
    return mat_labels, total_edge, clustPb
