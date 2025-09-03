#-*- coding: utf-8 -*-
# Pg assign a distance group for each packing defect
# R. Gerard 2024
# M. Zygadlo 2025

import numpy as np

from bin import matrix as m

def find_protein(matrix, defect, protein, arrayX, arrayY, zmean):
    """
    Locate the protein in the matrix.

    --------------------
    INPUT
    matrix : numpy array 2D
        Empty matrix (filled with 0)
    defect : string
        Take 'up' or 'lo' value. If we are working with the up or lo matrix
    protein : MDAnalysis group
        The protein(s)
    arrayX : numpy array
        List from xmin-1 to xmax+1 by step of 1.0
    arrayY : numpy array
        List from ymin-1 to ymax+1 by step of 1.0
    zmean : float
        The mean z position of the membrane

    --------------------
    OUTPUT
    numpy array 2D
        Contains the position of the protein(s) marked by 1
    """
    for i in range(len(protein.resids)) :
        # Get the coordinates X Y Z
        coordtmp = [protein.positions[i,0].round(2), protein.positions[i,1].round(2), protein.positions[i,2].round(2)]
        iX,iY = m.find_X_Y(coordtmp, arrayX, arrayY)
        # Upper leaflet
        if coordtmp[2] > zmean and defect == "up":
            matrix[iX, iY] = 1
        # Lower leaflet
        if coordtmp[2] < zmean and defect == "lo":
            matrix[iX, iY] = 1
    return matrix

def find_edges(mat):
    """
    Find the border coordinates of packing defects based on their labels.

    --------------------
    INPUT
    mat : numpy array 2D
        Contains the packing defect labels or the value 0
        OR
        Contains the position of the protein(s) marked by 1 (0 otherwise)
    --------------------
    OUTPUT
    dictionary
        Contains labels as keys and values are arrays of coordinates of the edge cells for each defect (label)
    """
    # Get the indexes of where the defects are
    defect_X = np.where(mat != 0)[0]
    defect_Y = np.where(mat != 0)[1]
    # Dictionary to store the coordinates of the packing defects' edges based on their label
    defects_edges_coors = {}

    for i in range(len(defect_X)):
        indX = defect_X[i]
        indY = defect_Y[i]
        # If this cell is a defect edge
        if (    mat[(indX-1) % mat.shape[0], indY] == 0 or  # top
                mat[(indX+1) % mat.shape[0], indY] == 0 or  # bottom
                mat[indX, (indY-1) % mat.shape[1]] == 0 or  # left
                mat[indX, (indY+1) % mat.shape[1]] == 0     # right
            ):
            # Get the label
            label = mat[indX, indY]
            # Add the coordinates of the edge cell to the dictionary with the corresponding label
            if label not in defects_edges_coors:
                defects_edges_coors[label] = np.array([np.array([indX, indY])])
            else:
                defects_edges_coors[label] = np.vstack((defects_edges_coors[label], np.array([indX, indY])))
    return defects_edges_coors

def find_shortest_sqdist(coor_prot_edge, coor_defect_edge):
    """
    Find the shortest squared distance between the protein edge coordinates and packing defects one.

    --------------------
    INPUT
    coor_prot_edge: numpy array
        Contains the coordinates of the protein edge cells
    coor_defect_edge: numpy array
        Contains the coordinates of the packing defect edge cells

    --------------------
    OUTPUT
    float
        The shortest squared distance between the edge coordinates of the protein and packing defects
    """
    # Compute all the square distances between tthe two arrays
    square_dist = np.sum((coor_prot_edge[:, None] - coor_defect_edge[None, :])**2, axis=2)
    # Get the minimum one
    min_sqdist = np.min(square_dist)
    return min_sqdist

def assign_dist_group(coor_prot_edge, dict_coor_defect_edge, sqdist_thres):
    """
    Assign distance groups to packing defects based on their proximity to the protein.

    --------------------
    INPUT
    coor_prot_edge: numpy array
        Contains the coordinates of the protein edge cells
    dict_coor_defect_edge: dictionary
        Contains labels as keys and arrays of defects edge cells coordinates as values
    sqdist_thres: int
        Squared distance threshold to determine if a packing defect is 'close' or 'far' from the protein

    --------------------
    INPUT
    dictionary
        Contains labels as keys and group ('close'/'far') based on the distance threshold as values
    """
    # Create a dictionnary to associate each label with their distance group
    defects_labels_group = {}

    for label in dict_coor_defect_edge:
        # Find the shortest squared distance to the protein
        dist2_tmp = find_shortest_sqdist(coor_prot_edge, dict_coor_defect_edge[label])
        # Assign a distance group to the label
        if dist2_tmp < sqdist_thres:
            defects_labels_group[label] = 'close'
        else :
            defects_labels_group[label] = 'far'
    return defects_labels_group

def outputTXT_defects_prot(out_name, dict_labels_group, dict_defect_area):
    """
    Create output file with information about the packing defects and their distance groups.

    --------------------
    INPUT
    out_name: string
        The output file name
    dict_labels_group: dictionary
        Contains labels as keys and group ('close' or 'far') as values
    dict_defect_area: dictionary
        Contains labels as keys and areas of each defect as values
    """
    with open(f"{out_name}.txt","w") as f:
        for label in dict_labels_group:
            if label in dict_defect_area.keys():
                f.write(f"{label},{dict_labels_group[label]},{dict_defect_area[label]}\n")
