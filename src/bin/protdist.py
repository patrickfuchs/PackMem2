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
        Contains labels as keys and values are lists of coordinates of the edge cells for each defect (label)
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
                defects_edges_coors[label] = []
            defects_edges_coors[label].append([indX, indY])                
    return defects_edges_coors


def find_short_dist2(prot_bord_coor, pd_bord_coor):

    """Find the shortest squared distance between the border coordinates of the protein and packing defects.

    Args:
        prot_bord_coor (list): List of coordinates (tuples) of the protein border cells.
        pd_bord_coor (list): List of coordinates (tuples) of the packing defect border cells.

    Returns:
        float: The shortest squared distance between the border coordinates of the protein and packing defects.
    """


    # Initialize a minimum squared distance with the first coordinates, arbitrarily using the first coordinates from each of the two lists.
    dist2 = (prot_bord_coor[0][0]-pd_bord_coor[0][0])**2 + (prot_bord_coor[0][1]-pd_bord_coor[0][1])**2

    # Calculate all possible squared distances between the coordinates of the two lists and return the smallest one.
    for prot_coor in prot_bord_coor :
        for pd_coor in pd_bord_coor :
            dist2_tmp = (prot_coor[0]-pd_coor[0])**2 + (prot_coor[1]-pd_coor[1])**2

            if dist2_tmp < dist2 :
                dist2 = dist2_tmp

    return dist2


def assign_dist_group(prot_bord_coor, dico_pd_bord_coor, dist_thr):

    """Assign distance groups to packing defects based on their proximity to the protein.

    Args:
        prot_bord_coor (list): List of coordinates (tuples) of the protein border cells.
        dico_pd_bord_coor (dict): Dictionary where keys are labels and values are lists of coordinates (tuples) of the packing defect border cells.
        dist_thr (int): Distance threshold to determine if a packing defect is 'close' or 'far' from the protein.

    Returns:
        dict: A dictionary where keys are labels and values are 'close' or 'far' based on the distance threshold.
    """

    if np.sum(prot_bord_coor) == 0:
        print(f"assign_dist_group : No protein found") # If the protein is not in the given leaflet.
        return 0

    # Create a dictionnary to associate each label with their distance group. {key=label : value=dist_group}
    pd_labels_group = {}

    # Set squared distance threshold
    dist2_thr = dist_thr**2

    # For each label, find the shortest squared distance to the protein. Given this result and the threshold, assign a distance group to the label.
    for label in dico_pd_bord_coor:
        dist2_tmp = find_short_dist2(prot_bord_coor, dico_pd_bord_coor[label])

        if dist2_tmp < dist2_thr:
            pd_labels_group[label] = 'close'
        else :
            pd_labels_group[label] = 'far'

    return pd_labels_group

def outputTXT_defects_prot(outputname, FlagPDtype, leaflet, dico_labels_group, dico_def_area):
    """Output a text file with information about the packing defects and their distance groups.

    Args:
        outputname (str): The base name of the output file.
        FlagPDtype (str): The type of packing defects.
        leaflet (str): The leaflet ('Up' or 'Lo') being processed.
        dico_labels_group (dict): Dictionary where keys are labels and values are 'close' or 'far' based on the distance threshold.
        dico_def_area (dict): Dictionary where keys are labels and values are the areas of the defects.
    """
    flag_to_defect = ["All", "Deep", "Shallow"]
    outputname = outputname + "_" + leaflet+ "_"+ FlagPDtype +"_prot.txt"

    if dico_labels_group == 0:
        with open(outputname,"w") as f:
            f.write(f"0,NA,0\n")

    with open(outputname,"w") as f:
        for label in dico_labels_group:
            if label in dico_labels_group.keys() and label in dico_def_area.keys():
                f.write("%s,%s,%s\n" % (label, dico_labels_group[label], dico_def_area[label]))
