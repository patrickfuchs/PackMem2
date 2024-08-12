#-*- coding: utf-8 -*-
# Funtions for matrix
# R. Gautier A. Bacle 2015

import math
import sys
import numpy as np
from bin import BasicFunctions as bfrg

# matrix size  (square size = 1A)
SIZE=1.0
# calculate limit size (half diagonal) 3D
SIZE_SIDE=0.5*(math.sqrt((math.pow(SIZE, 2))* 3.))


def initialize_matrix2D(val1, val2, default):
    """
    Initialise a 2D matrix with a default value
    
    --------------------
    INPUT
    val1 : int
        Dimension in X/Y for the array
    val2 : int
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

def get_radius(radius, res_name, atom_name):
    """
    Get radius from [lipid_name atom_name]

    --------------------
    INPUT
    radius : dictionnary
        Contains the atom's radiuses of lipids and amino acids
    res_name : string
        Name of the residue (lipid / amino acid)
    atom_name : string
        Name of the atom
    
    --------------------
    OUTPUT
    float
        The radius of an lipid's / amino acid's atom
    """
    key = res_name+" "+atom_name
    if key not in radius.keys():
        print(f"ERROR : Association of lipid {res_name} and atom {atom_name} not found in parameter file")
        sys.exit()
    return radius[key]

def get_aliphatic(aliphatic, res_name, atom_name):
    """
    Get aliphatic flag [Lipid_name Atom_name]

    n for polar atoms
    a for aliphatic atoms

    --------------------
    INPUT
    aliphatic : dictionnary
        Contains the nature of the atoms
    res_name : string
        Name of the residue (lipid / amino acid)
    atom_name : string
        Name of the atom
    
    --------------------
    OUTPUT
    float
        The nature of an lipid's / amino acid's atom
    """
    key = res_name+" "+atom_name
    if key not in aliphatic.keys():
        print(f"ERROR : Association of lipid {res_name} and atom {atom_name} not found in parameter file")
        sys.exit()
    return aliphatic[key]

def find_Z(coordZ, listZ):
    """
    Compute the difference between the last value in listZ
    and the Z for one atom

    --------------------
    INPUT
    coordZ : float
        The z coordinate for one atom
    listZ : list
        List from zmax+1 to z_C2_coord-1 by step of 1.0  OR
        List from z_C2_coord+1 to zmin-1 by step of 1.0

    --------------------
    OUPUT
    float
        The difference between the last value in listZ and the Z for one atom
    """
    return listZ[-1]-coordZ

def find_X_Y(coord, matX, matY):
    """
    Determine from x,y coordinates the index in matX and matY
    Has  been said to determine the central value in matrix

    --------------------
    INPUT
    coord : list
        Contains the coordinates x,y,z of an atom
    matX : list
        Contains integer from xmin-1 to xmax+1 by step of 1.0
    matY : list
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
    iX = matX.index(tmpX)
    iY = matY.index(tmpY)
    return iX,iY

def check_edges(val, val_lim1, val_lim2):
    """
    Check if the cell index is lower than 5 or higher than dimension-5
    Changes the value of the cell index if so

    --------------------
    INPUT
    val : int
        Index of the coordinate of the atom in listX or listY
    val_lim1 : int
        The number of cell to work around the atom
    val_lim2 : int
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

def setDefects(type_aliphatic, val_mat):
    """
    Fill matrix cell depending on the defect type (Deep, Shallow)

    --------------------
    INPUT
    type_aliphatic : string
        The nature of the atom. a : aliphatique / n : polar
    val_mat : float
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

def fill_matrix(matrix, coordtmp, listX, listY, listZ,
                radius_res, FlagPDtype, aliph_atoms):
    """
    fill the matrix for each atom
    ---------------
    INPUT
    matrix : numpy matrix
        Matrix to be filled where there are atoms
    coordtmp : list
        Contains the coordinates x,y,z of the atom
    listX : list
        List from xmin-1 to xmax+1 by step of 1.0
    listY : list
        List from ymin-1 to ymax+1 by step of 1.0
    listZ : list
        List from zmax+1 to z_C2_coord-1 by step of 1.0  OR
        List from z_C2_coord+1 to zmin-1 by step of 1.0
    radius_res : float
        the radius of the atom type for the residue
    FlagPDtype : int
        The type of defect to analyse. 0 : all / 1 : deep / 2 : shallow
    aliph_atoms : string
        The nature of the atom. a : aliphatique / n : polar

    -----------------
    OUTPUT
    numpy matrix
        Matrix filled where there are atoms
    """
    # Number of cells to work around
    v=5
    # Find the index of x_atom and y_atom in listX and listY
    # this corresponds to the location of the atom in the matrix
    iX,iY = find_X_Y(coordtmp, listX, listY)
    # Change the indexes (cell index in matrix) if < v or > len(matrix)-5
    iX = check_edges(iX, v, len(listX))
    iY = check_edges(iY, v, len(listY))
    # Select the cells to work in at i+-v
    listXM = listX[iX - v : iX + (v + 1)]
    listYM = listY[iY - v : iY + (v + 1)]
    # Limit distance in a radius
    dist_lim = (SIZE+radius_res) ** 2
    # Limit distance in diagonal
    dist_meet = (SIZE_SIDE + radius_res) ** 2
    # Loop on the different z positions of the upper OR lower leaflet
    for sliceZ in listZ:
        # Distance from the position of the atom to the slice in Z
        distZ_to_sliceZ = bfrg.dist_oneAxis(coordtmp[2],sliceZ)
        # Check that we can search in the all the radius
        if distZ_to_sliceZ > dist_lim:
            continue
        # Loop on the different cells in x dimension
        for indX,sliceX in enumerate(listXM):
            # Distance from the position of the atom to a cell in X
            distX_to_sliceX = bfrg.dist_oneAxis(coordtmp[0], sliceX)
            # Check that we can search in the all the radius
            if distX_to_sliceX > dist_lim:
                continue
            # Loop on the different cells in y dimension
            for indY,sliceY in enumerate(listYM):
                # Distance from the position of the atom to a cell in X
                distY_to_sliceY = bfrg.dist_oneAxis(coordtmp[1], sliceY)
                # Check that we can search in the all the radius
                if distY_to_sliceY > dist_lim:
                    continue
                # 
                X = indX + (iX - v)
                Y = indY + (iY - v)
                # Create a list with the coordinates x,y,z of a cell
                # in a 5 cell radius of the atom
                coordCenter = [float(sliceX), float(sliceY), float(sliceZ)]
                # Compute the distance between the atom and this position of the matrix
                distance = bfrg.dist(coordCenter, coordtmp)
                # If the matrix cell was empty, put 0.0
                if np.isnan(matrix[X][Y]):
                    matrix[X][Y]=0.
                # Get the value in this cell of the matrix
                val_mat = matrix[X][Y]
                # 
                if distance <= dist_meet:
                    # If shallow (2) or all (0) defect
                    if FlagPDtype == 2 or FlagPDtype == 0:
                        # 
                        matrix[X][Y] = setDefects(aliph_atoms, val_mat)
                    # If deep (1) defect
                    else:
                        matrix[X][Y] +=1.
    return matrix

def binarize_matrix_without0(matrix, matrix_ini, val1=0, val2=0.99):
    """
    Binarise the presence of aliphatic atom (and packing defects) in matrix

    --------------------
    INPUT
    matrix : numpy matrix
        Contains where there are atoms in the simulation box and their type
    matrix_ini : numpy matrix
        Matrix to be binarised
    val1 : float
        lower limit value
    val2 : float
        higher limit value

    --------------------
    OUTPUT
    numpy matrix
        Contains the position of the aliphatic atoms (+ packing defects) (0) in the simulationb box    
    """
    for i in range(0, matrix.shape[0]):
        for j in range(0, matrix.shape[1]):
            if np.isnan(matrix[i][j]):
                matrix_ini[i][j]=0.
            # Check Matrix, if there is an aliphatic atom (or pacing defect) : 0
            # Otherwise : 1
            else:
                if (matrix[i][j] > val1 and  matrix[i][j] < val2):
                    matrix_ini[i][j]=0.
                else:
                    matrix_ini[i][j]=1.
    return matrix_ini

def binarize_matrix(matrix, matrix_ini, val=0.):
    """
    Binarise the presence of atom in matrix

    --------------------
    INPUT
    matrix : numpy matrix
        Contains where there are atoms in the simulation box
    matrix_ini : numpy matrix
        Matrix to be binarised
    val : float
        lower limit value

    --------------------
    OUTPUT
    numpy matrix
        Contains the position of the atoms (1) in the simulation box
    """
    for i in range(0, matrix.shape[0]):
        for j in range(0, matrix.shape[1]):
            #print(matrix[i][j])
            if np.isnan(matrix[i][j]):
                matrix_ini[i][j]=0.
            else:
                if matrix[i][j] > val:
                    matrix_ini[i][j]=1.
                else:
                    matrix_ini[i][j]=0.
    return matrix_ini

def modify_matrix(mat1, mat2, listval):
    """
    Modify the binary matrix to take account edges (determined by all packing defects)
    The edges are put to 0.0

    --------------------
    INPUT
    mat1 : numpy matrix
        Contains the labels of the packing defects + aliphatic atoms
    mat2 : numpy matrix
        Contains the position of the aliphatic atoms (0)
    listval : list
        Contains the labels of the clusters on the edges
    
    --------------------
    OUTPUT
    matrix
        Contains mat2 where the clusters on the edge were put to 0.0
    """
    for i in range(0, mat1.shape[0]):
        for j in range(0, mat1.shape[1]):
            if mat1[i][j] in listval:
                mat2[i][j]= 0.0
    return mat2

def initialize_matrix3D(val1, val2, val3, defaut):
    Matrix=list(range(0,val1))
    for i in range(0, val1):
        Matrix[i]=list(range(0,val2))
        for j in range(0, val2):
            Matrix[i][j]=list(range(0,val3))
            for k in range(0, val3):
                Matrix[i][j][k]=defaut
    return Matrix

def binarize_matrix_whithout0_bis(matrix, matrix_ini, val=0.99):
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix[0])):
            #print(matrix[i][j])
            if matrix[i][j] == "NA":
                matrix_ini[i][j]=0.
            else:
                if (matrix[i][j] > 0. and  matrix[i][j] < val):
                    matrix_ini[i][j]=0.
                else:
                    matrix_ini[i][j]=1.
    return matrix_ini

def clean_NA_inside(Matrix_labels, clust_edge, Matrix_ini, total_edge):
    """

    --------------------
    INPUT
    Matrix_labels : numpy matrix
        Contains the labels
    clust_edge : list
        Contains the labels on each edge of the matrix
    Matrix_ini : numpy matrix
        Contains the positions of the polar atoms (int)
        aliphatic atoms (float > 0.0)
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
        Contains the labels and the number of cells concerned that and nan in Matrix_ini
    """
    clustPb={}
    # Loop on the Matrix_ini to find potential nan
    for i in range(0, Matrix_ini.shape[0]):
        for j in range(0, Matrix_ini.shape[1]):
            if np.isnan(Matrix_ini[i][j]) and Matrix_labels[i][j] not in clust_edge:
                if Matrix_labels[i][j] in clustPb:
                    clustPb[Matrix_labels[i][j]]+=1
                else:
                    clustPb[Matrix_labels[i][j]]=1
                # Correct the nan by giving it the first label of clust_edge
                Matrix_labels[i][j]=clust_edge[0]
                total_edge+=1
    return Matrix_labels, total_edge, clustPb
            
def determinelastNA(matrix):
    for i in range(len(matrix)-1, -1, -1):
        print(i, matrix[i][0][0])

def find_indexfromvalue(value, liste):
    mini=100.
    mini_i=len(liste)
    for i, val in enumerate(liste):
        if math.fabs(val-value) < mini:
            mini=math.fabs(val-value)
            mini_i=i
    return mini_i
    
