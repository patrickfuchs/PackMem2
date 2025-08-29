# -*- coding: utf-8 -*-

##### Here is implemented the Connected-component_labeling algorithm
##### Source: http://en.wikipedia.org/wiki/Connected-component_labeling
#####
##### P. Fuchs 08/2014
#####

import numpy as np
import sys
from bin import matrix as m

def get_uniq_labels(nb_labels,list_equiv_labels):
    # list_equiv_labels looks like that: [[1,2], [3,4,5,8], [6,7]...]
    # At the end we want to get something like this:
    # {1: 1, 2: 1, 3: 3, 4: 3, 5: 3, 6: 6, 7: 6, 8:3}
    dico_uniq_labels = {} ; root_labels = []
    for label in range(1,nb_labels+1):
        # is label in list_equiv_labels?
        is_label_in_list_equiv_labels = 0
        for tmp_list_labels in list_equiv_labels:
            if label in tmp_list_labels:
                # assign the root label to the current label
                dico_uniq_labels[label] = tmp_list_labels[0]
                is_label_in_list_equiv_labels = 1
                # is it a root label (thus at the first position of the list)?
                if label == tmp_list_labels[0]:
                    if label not in root_labels: # in case of big cluster, it avoids writing multiple times the same label!
                        root_labels.append(label)
        if not is_label_in_list_equiv_labels:
            # if we end up overhere, the label is a root label but unique
            dico_uniq_labels[label] = label
            root_labels.append(label)
    DEBUG=0
    if DEBUG:
        print("list_equiv_labels:" , list_equiv_labels)
        print("nb labels=%i" % nb_labels)
        print("Dictionnary of equivalence:")
        print(dico_uniq_labels)
        print("nb of root labels=%i" % len(root_labels))
        print("root labels are:"); print(root_labels)
        #exit()
    return dico_uniq_labels, root_labels

# returns the intersection of 2 lists
def intersect(l1, l2):
    return list(set(l1) & set(l2))

def merge_avoid_duplicate_sort_list(l1,l2):
    new_l = list(set(l1 + l2))
    new_l.sort()
    return new_l

def is_duplicate_in_list_of_lists(list_of_lists, starting_index):
    L = len(list_of_lists)
    for i in range(starting_index+1,L):
        if len(intersect(list_of_lists[starting_index],list_of_lists[i])) > 0:
            return True
    return False

def is_empty_sublist(list_of_lists):
    L = len(list_of_lists)
    for i in range(L):
        if len(list_of_lists[i]) == 0:
            return True
    return False
    
# merge a list of list
def merge_list_of_lists(list_of_lists):
    # we want to merge the sublists that intersect in a list_of_lists, e.g.
    # init list_of_lists:  [[1, 2, 3, 9], [4, 5, 8], [1, 7, 9], [8, 9, 10], [13, 16], [12, 44], [16, 54]]
    # final list_of_lists: [[1, 2, 3, 4, 5, 7, 8, 9, 10], [13, 16, 54], [12, 44]]
    #print "Intial list_of_lists", list_of_lists
    L = len(list_of_lists)
    starting_index = 0
    while starting_index <= L - 1:
        # we try to merge sublist starting_index with all other sublists
        while is_duplicate_in_list_of_lists(list_of_lists,starting_index):
            for i in range(starting_index+1,L):
                if len(intersect(list_of_lists[starting_index],list_of_lists[i])) > 0:
                    # OK we merge sublist starting_index with sublist i
                    list_of_lists[starting_index] = merge_avoid_duplicate_sort_list(list_of_lists[starting_index], list_of_lists[i])
                    # we delete sublist i
                    list_of_lists[i] = []
        # we increment starting index
        starting_index += 1
    # cleanup empty sublists
    while is_empty_sublist(list_of_lists):
        L = len(list_of_lists)
        for i in range(L):
            if len(list_of_lists[i]) == 0:
                list_of_lists.pop(i)
                break
    #print "Final list_of_lists" , list_of_lists
    return list_of_lists


# this is the important fct 
def get_connected_components(M, val_bin=0):
    """
    Connect the shallow correspond cells together

    --------------------
    INPUT
    M : numpy matrix
        Contains the position of the aliphatic atoms (+ packing defects) (0) in the simulationb box
    val_bin : int / float
        The value to differenciate between pore or atom

    --------------------
    OUPUT
    numpy matrix
        Contains the labels
    list
        Contains the set of labels in this matrix
    dictionnary
        Contains the area of each label
    dictionnary
        Contains the first appearance of the label in the matrix
    """
    # Initialize a new matrix (for labels)
    # Set each cell to 0 (integer)
    M_labels = m.initialize_matrix2D(M.shape[0], M.shape[1], 0)
    # counter for the number of labels
    nb_labels = 0 
    # list for storing equivalent labels (see below)
    list_equiv_labels = []

    # First pass
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            # If there is a positive cell (== 0)
            if M[i][j] == val_bin:
                # create a list that contains the positive neighbors
                # e.g.: [[1, 2]] means the neighboring cell [1,2] (of the current cell) is positive
                # e.g.: [[6, 6], [6, 4]] means the neighboring cells (of the current cell)
                #       [6,6] and [6,4] are positive
                positive_neighbors = []
                # look for connectivity of the current positive cell
                if i > 0:
                    if j > 0 and M[i-1][j-1] == val_bin: # North-West
                        positive_neighbors.append([i-1,j-1])                
                    if M[i-1][j] == val_bin: # North
                        positive_neighbors.append([i-1,j])
                    if j < M.shape[1] - 1 and M[i-1][j+1] == val_bin: # North-East
                        positive_neighbors.append([i-1,j+1])
                if j > 0 and M[i][j-1] == val_bin: # West
                    positive_neighbors.append([i,j-1])
                # if no neighbor is positive, this cell is a new label
                if len(positive_neighbors) == 0:
                    nb_labels += 1
                    M_labels[i][j] = nb_labels
                # if only 1 neighbor is positive, assign the neighbor label to current cell
                if len(positive_neighbors) == 1:
                    row_uniq_neighbor = positive_neighbors[0][0]
                    col_uniq_neighbor = positive_neighbors[0][1]
                    M_labels[i][j] = M_labels[row_uniq_neighbor][col_uniq_neighbor]
                # if multiple neighbors are positive, assign one of their label
                #  --> it doesn't matter which one: arbitrarily, we choose the first one
                # in the list positive_neighbors
                if len(positive_neighbors) > 1:
                    # loop over positive_neighbors and store their label in a list called label_neighbors
                    label_neighbors = []
                    for neighbor in positive_neighbors: # neighbor is a list (e.g. [11,81]) containing coor of the neighbor in M
                        row_neighbor = neighbor[0]
                        col_neighbor = neighbor[1]
                        label_neighbors.append(M_labels[row_neighbor][col_neighbor])
                    # eliminate duplicate & sort
                    label_neighbors = list(set(label_neighbors))
                    label_neighbors.sort()
                    # assign the label (min of label_neighbors) to that cell
                    # I tested, any label label_neighbors work as well!
                    M_labels[i][j] = min(label_neighbors)
                    
                    # in case of multiple labels, keep track they are equivalent in list_equiv_labels:
                    # e.g. [[1,2],[3,4,5,7],[6,8], ...] (each sublist is ordered)
                    #  -> means 1&2 are equiv, 3,..,7 are equiv, etc
                    if len(label_neighbors) > 1:
                        # is list_equiv_labels empty? -> easy to fill in
                        if len(list_equiv_labels) == 0:
                            list_equiv_labels.append(label_neighbors)
                        # otherwise, it's not empty! -> and we have to fill it in
                        else:
                            # is one of the label_neighbors present in list_equiv_labels?
                            is_present_in_list_equiv_labels = 0
                            for k,tmplist in enumerate(list_equiv_labels):
                                if len(intersect(tmplist,label_neighbors)) > 0:
                                    # found it! So we fuse the lists, avoid duplicates and sort!
                                    list_equiv_labels[k] = merge_avoid_duplicate_sort_list(list_equiv_labels[k], label_neighbors)
                                    is_present_in_list_equiv_labels = 1
                            # label is not present, so we create a new sublist
                            if not is_present_in_list_equiv_labels:
                                list_equiv_labels.append(label_neighbors)

    # There are redundancies in list_equiv_labels -> clean them up
    list_equiv_labels = merge_list_of_lists(list_equiv_labels)

    # Second pass
    # Build a dictionnary for assigning a uniq label
    # should look like this: {1: 1, 2: 1, 3: 3, 4: 4}
    # --> label 1 & 2 are equiv, thus we set label 2 to label 1
    dico_of_uniq_labels, root_labels = get_uniq_labels(nb_labels,list_equiv_labels)
    area_clusters = {} # to store the area of each cluster
    coor_clusters = {} # to store the coor (in the matrix) of the 1st element of each cluster
    # Loop on the labels that are not equivalent to any other
    for root_label in root_labels:
        area_clusters[root_label] = 0
    # Loop on label matrix to get rid of equivalent labels
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            # if there is a label
            if M_labels[i][j]:
                # Change label to its equivalent and kept one
                M_labels[i][j] = dico_of_uniq_labels[M_labels[i][j]]
                area_clusters[M_labels[i][j]] += 1
                # Store the coordinates of the first appearance of the label
                if M_labels[i][j] not in coor_clusters:
                    coor_clusters[M_labels[i][j]]=[i,j]
                
    return M_labels, root_labels, area_clusters, coor_clusters

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
    
