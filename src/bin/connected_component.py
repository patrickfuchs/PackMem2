# -*- coding: utf-8 -*-
##### Here is implemented the Connected-component_labeling algorithm
##### Source: http://en.wikipedia.org/wiki/Connected-component_labeling
# P. Fuchs 08/2014
# M. Zygadlo 2025

import numpy as np
import sys
from bin import matrix as m

def find_neighbours(matrix, ind_i, ind_j):
    """
    Find the defect neighbour cells in the upper left corner of a defect cell.

    --------------------
    INPUT
    matrix: numpy array 2D
        Contains the positions of the defects (0), 1 otherwise
    ind_i : int
        The index i in the matrix
    ind j : int
        The index j in the matrix
    
    --------------------
    OUTPUT
    list
        Contains the index of the defect neighbours of the given cell
    """
    defect_neighbours = []
    # Look for connectivity of the current cell
    if ind_i > 0:
        if ind_j >0 and matrix[(ind_i-1) % matrix.shape[0], (ind_j-1) % matrix.shape[1]] == 0: # North-West
            defect_neighbours.append([(ind_i-1) % matrix.shape[0], (ind_j-1) % matrix.shape[1]])                
        if matrix[(ind_i-1) % matrix.shape[0], ind_j] == 0: # North
            defect_neighbours.append([(ind_i-1) % matrix.shape[0], ind_j])
        if matrix[(ind_i-1) % matrix.shape[0], (ind_j+1) % matrix.shape[1]] == 0: # North-East
            defect_neighbours.append([(ind_i-1) % matrix.shape[0], (ind_j+1) % matrix.shape[1]])
    if ind_j > 0 and matrix[ind_i, (ind_j-1) % matrix.shape[1]] == 0: # West
        defect_neighbours.append([ind_i, (ind_j-1) % matrix.shape[1]])
    return defect_neighbours

def find_label_neighbours(matrix, defect_neighbours):
    """
    Find the label of the neighbours coordinates given.
    
    --------------------
    INPUT
    matrix: numpy array 2D
        Contains the labels of the defects
    defect_neighbours: list of lists
        Contains the indexes of the neighbouring cells
    
    --------------------
    OUTPUT
    list
        Contains the sorted labels of the neighbouring cells without duplicates
    """
    label_neighbours = []
    # Loop over the neighbours
    for neighbour in defect_neighbours:
        neighbour_i = neighbour[0]
        neighbour_j = neighbour[1]
        # Store their label in a list
        label_neighbours.append(matrix[neighbour_i, neighbour_j])
    # Eliminate duplicate & sort
    label_neighbours = list(set(label_neighbours))
    label_neighbours.sort()
    return label_neighbours

def intersect(l1, l2):
    """
    Return the intersection of 2 lists.

    --------------------
    INPUT
    l1: list
    l2: list
    
    --------------------
    OUTPUT
    list
        Contains the set of the intersection of the two lists
    """
    return list(set(l1) & set(l2))

def sort_without_duplicate_merged_lists(l1, l2):
    """
    Return the intersection of 2 lists.

    --------------------
    INPUT
    l1: list
    l2: list
    
    --------------------
    OUTPUT
    list
        Contains a sorted set of the two lists added
    """
    new_list = list(set(l1 + l2))
    new_list.sort()
    return new_list
 
def merge_list_of_lists(list_of_lists):
    """
    Merge a list of lists.

    merge the sublists that intersect in a list_of_lists, e.g.
    init list_of_lists:  [[1, 2, 3, 9], [4, 5, 8], [1, 7, 9], [8, 9, 10], [13, 16], [12, 44], [16, 54]]
    final list_of_lists: [[1, 2, 3, 4, 5, 7, 8, 9, 10], [13, 16, 54], [12, 44]]
    
    --------------------
    INPUT
    list_of_lists: list of lists
    
    --------------------
    OUTPUT
    list of lists
        Contains the sublists merged with the ones that intersectted with them
    """
    for i in range(len(list_of_lists)-1, -1, -1):
        for j in range(i-1, -1, -1):
            if intersect(list_of_lists[j], list_of_lists[i]):
                list_of_lists[j] = list(set(list_of_lists[i] + list_of_lists[j]))
                list_of_lists[j].sort()
                del list_of_lists[i]
                break
    return list_of_lists

def get_uniq_labels(nb_labels, equiv_labels):
    """
    Build a dictionnary for assigning a uniq label.

    If equiv_labels = [[1,2], [3,8,9], [6,7]]
    At the end we want to get something like this:
    {1: 1, 2: 1, 3: 3, 4: 4, 5: 5, 6: 6, 7: 6, 8:3, 9:3}
    [1, 3, 4, 5, 6]

    --------------------
    INPUT
    nb_labels: int
        The number of labels
    equiv_labels: list of lists
        Contains the equivalent labels in each sublists
    
    --------------------
    OUTPUT
    dictionary
        Contains the labels as keys and the unique minimal label for each equivalent labels as values
    list
        Contains the unique minimal label for each equivalent labels
    """
    dict_connected_labels = {}
    uniq_labels = []
    for label in range(1, nb_labels+1):
        for sublist in equiv_labels:
            if label in sublist:
                dict_connected_labels[label] = sublist[0]
                break
            dict_connected_labels[label] = label
    uniq_labels = list(set(dict_connected_labels.values()))
    return dict_connected_labels, uniq_labels

def connect_labels_in_matrix(mat_labels, dict_connected_labels):
    """
    Change the labels in mat_labels to the smallest one they are connected to.

    If tthe input is:
    mat = np.array([[0, 0, 0, 0, 0, 0, 0],
                    [0, 1, 2, 0, 0, 8, 0],
                    [0, 1, 0, 0, 3, 9, 0],
                    [0, 0, 0, 0, 0, 9, 0],
                    [0, 6, 0, 4, 0, 0, 0],
                    [0, 7, 0, 0, 0, 5, 0],
                    [0, 0, 0, 0, 0, 0, 0]])
    dict_connected_labels = {1: 1, 2: 1, 3: 3, 4: 4, 5: 5, 6: 6, 7: 6, 8: 3, 9: 3}
    The output is:
    np.array([[0, 0, 0, 0, 0, 0, 0],
            [0, 1, 1, 0, 0, 3, 0],
            [0, 1, 0, 0, 3, 3, 0],
            [0, 0, 0, 0, 0, 3, 0],
            [0, 6, 0, 4, 0, 0, 0],
            [0, 6, 0, 0, 0, 5, 0],
            [0, 0, 0, 0, 0, 0, 0]])
    --------------------
    INPUT
    mat_labels: numpy array 2D
        Contains the raw labels of each defect
    dict_connected_labels: dictionary
        Contains the labels as keys and the unique minimal label for each equivalent labels as values

    --------------------
    OUTPUT
    numpy array 2D
        Contains the smallest label the previous label was connected to for each defect
    """
    for label in dict_connected_labels:
        mat_labels = np.where(mat_labels == label, dict_connected_labels[label], mat_labels)
    return mat_labels

def get_area_first_coor_defects(matrix, uniq_labels):
    """
    Get the area of each defect and their first coordinates.

    If the input is:
    mat = np.array([[0, 0, 0, 0, 0, 0, 0],
                    [0, 1, 1, 0, 0, 3, 0],
                    [0, 1, 0, 0, 3, 3, 0],
                    [0, 0, 0, 0, 0, 3, 0],
                    [0, 6, 0, 4, 0, 0, 0],
                    [0, 6, 0, 0, 0, 5, 0],
                    [0, 0, 0, 0, 0, 0, 0]])
    uniq_labels = [1, 3, 4, 5, 6]
    The output will be:
    area_defects = {1:3, 3:4, 4:1, 5:1, 6:2}
    coor_defects = {1:[1,1], 3:[1,5], 4:[4,3], 5:[5,5], 6:[4,1]}

    --------------------
    INPUT
    mat: numpy array 2D
        Contains the labels of the defects
    uniq_labels: list
        Contains the unique minimal labels
    
    --------------------
    OUTPUT
    dictionary
        Contains the area of each defect
    dictionary
        Contains the first coordinates of each defect
    """
    area_defects = {}
    first_coor_defects = {}

    for label in uniq_labels:
        # Get the positions of this label
        label_X = np.where(matrix == label)[0]
        label_Y = np.where(matrix == label)[1]
        # Add the surface they cover to area_defects
        if label not in area_defects:
            area_defects[label] = 0
            # Take the first occurence as reference coordinates
            first_coor_defects[label] = [label_X[0], label_Y[0]]
        area_defects[label] += len(label_X)
    return area_defects, first_coor_defects

def get_connected_components(matrix):
    """
    Connect and label the defects. Get a list of the labels, tehir area and ttheir first coordinates.

    --------------------
    INPUT
    mat: numpy array 2D
        Contains the position of the aliphatic atoms + packing defects (marked by 0), 1 otherwise

    --------------------
    OUPUT
    numpy array 2D
        Contains the connected labels
    list
        Contains the set of labels in this matrix
    dictionnary
        Contains the area of each label
    dictionnary
        Contains the first appearance of the label in the matrix
    """
    mat_labels = m.initialize_matrix2D(matrix.shape[0], matrix.shape[1], 0)

    nb_labels = 0 
    equiv_labels = []

    # Find where there are defects
    defect_X = np.where(matrix == 0)[0]
    defect_Y = np.where(matrix == 0)[1]

    for i in range(len(defect_X)):
        indX = defect_X[i]
        indY = defect_Y[i]
        # Create a list that contains the coordinates of the defect neighbours
        defect_neighbours = find_neighbours(matrix, indX, indY)
        
        # If the cell is an upper left defect edge corner, this cell is a new label
        if len(defect_neighbours) == 0:
            nb_labels += 1
            mat_labels[indX, indY] = nb_labels
        # If the cell is part of the lower right defect edge
        if len(defect_neighbours) == 1:
            neighbour_i = defect_neighbours[0][0]
            neighbour_j = defect_neighbours[0][1]
            mat_labels[indX, indY] = mat_labels[neighbour_i, neighbour_j]
        # If the cell is inside a defect
        if len(defect_neighbours) > 1:
            # Find all the labels of the neighbours
            label_neighbours = find_label_neighbours(mat_labels, defect_neighbours)
            # assign the smallest label to that cell
            mat_labels[indX][indY] = min(label_neighbours)
            
            # In case of multiple labels
            if len(label_neighbours) > 1:
                if len(equiv_labels) == 0:
                    equiv_labels.append(label_neighbours)
                else:
                    is_present_in_equiv_labels = False
                    for k,tmplist in enumerate(equiv_labels):
                        # If one the labels already is in equiv_labels
                        if intersect(tmplist, label_neighbours):
                            # Fuse the lists, avoid duplicates and sort
                            equiv_labels[k] = sort_without_duplicate_merged_lists(equiv_labels[k], label_neighbours)
                            is_present_in_equiv_labels = True
                    # If the labels don't exist in equiv_labels
                    if not is_present_in_equiv_labels:
                        equiv_labels.append(label_neighbours)

    # If there are redundancies, clean them up
    equiv_labels = merge_list_of_lists(equiv_labels)

    # Get the labels to be linked to their smallest equivalent label
    # Get the list of the unique values in the dictionary
    dict_connected_labels, uniq_labels = get_uniq_labels(nb_labels, equiv_labels)

    # Change the labels in the matrix to the smallest one
    mat_labels = connect_labels_in_matrix(mat_labels, dict_connected_labels)

    # Get the area of each defect and their first coordinates
    area_defects, first_coor_defects = get_area_first_coor_defects(mat_labels, uniq_labels)
                
    return mat_labels, uniq_labels, area_defects, first_coor_defects

def get_clusters_on_the_edge(M_labels):
    """
    Get the labels that are on the edge of the matrix

    --------------------
    INPUT
    M_labels : numpy matrix
        Contains the labels of the packing defects (+ aliphatic atoms)
    
    --------------------
    OUTPUT
    list
        Contains the labels of each edge of the matrix
    """
    nrows = M_labels.shape[0]
    ncols = M_labels.shape[1]
    clusters_on_the_edge = []
    # Get set of labels on 1st & last row
    clusters_on_the_edge += list(set([i for i in M_labels[0] if i != 0]))
    clusters_on_the_edge += list(set([i for i in M_labels[nrows-1] if i != 0]))
    # check clusters on 1st & last col
    clusters_on_the_edge += list(set([i for i in M_labels[:, 0] if i != 0]))
    clusters_on_the_edge += list(set([i for i in M_labels[:, ncols-1] if i != 0]))
    #for i in range(nrows):
    #    if M_labels[i][0] not in clusters_on_the_edge and M_labels[i][0] != 0:
    #        clusters_on_the_edge.append(M_labels[i][0])
    #    if M_labels[i][ncols-1] not in clusters_on_the_edge and M_labels[i][ncols-1] != 0:
    #        clusters_on_the_edge.append(M_labels[i][ncols-1])
    # get rid of duplicates & sort
    clusters_on_the_edge = sorted(list(set(clusters_on_the_edge)))
    return clusters_on_the_edge

def delete_NApoints_inside(clustPb, Matrix_labels, root_labels, area_clusters):
    """
    Get rid of the defect comprised enterely of nan

    --------------------
    INPUT
    clustPb : dictionnary
        Contains the labels and the number of cells concerned that has nan in Matrix_ini
    Matrix_labels : numpy matrix
        Contains the labels
    root_labels : list
        Contains the set of labels in this matrix
    area_clusters : dictionnary
        Contains the area of each label

    --------------------
    OUTPUT
    list
        Contains the set of labels in this matrix that are attributed to a non-nan cluster
    dictionnary
        Contains the area of each non-nan label
    dictionnary
        Contains the first appearance of the non-nan labels in the matrix
    """
    clust_2_delete=[]
    # Loop on all the clusters that had nan in them
    for key in clustPb:
        #  If one of the clusters is enterely composed of nan
        if area_clusters[key] == clustPb[key]:
            # Add this label to the delete list
            clust_2_delete.append(key)
    # if there are some labels to be deleted
    if len(clust_2_delete) != 0:
        for num in clust_2_delete:
            # Get their indexes in root_label
            tmp=root_labels.index(num)
            # Delete them from root_label
            del root_labels[tmp]
    # to store the area of each cluster
    area_clusters = {}
    # to store the coor (in the matrix) of the 1st element of each cluster
    coor_clusters = {}
    # Loop on the packing labels in the matrix
    for root_label in root_labels:
        area_clusters[root_label] = 0
    # Loop on label matrix to compute the area of starting cell of cluster
    for i in range(Matrix_labels.shape[0]):
        for j in range(Matrix_labels.shape[1]):
            # if there is a label
            if Matrix_labels[i][j]:
                # Add this cell to the area of the label
                area_clusters[Matrix_labels[i][j]] += 1
                # Store the coordinates of the first appearance of the label
                if Matrix_labels[i][j] not in coor_clusters:
                    coor_clusters[Matrix_labels[i][j]]=[i,j]
    return root_labels, area_clusters, coor_clusters
    
