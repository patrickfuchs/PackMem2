#-*- coding: utf-8 -*-
# Functions about dictionaries
# R. Gautier 2015
# M. Zygadlo 2025

import sys

def get_value(dict_info, resname, atom_name):
    """
    Get value from a dictionary of the form "lipid_name atom_name" : value.

    --------------------
    INPUT
    dict_info: dictionnary
        Contains the atom radius or nature (n for polar / a for apolar)
    resname: string
        Name of the residue (lipid / amino acid)
    atom_name: string
        Name of the atom
    
    --------------------
    OUTPUT
    float
        The atom radius or nature
    """
    key = resname+" "+atom_name
    if key not in dict_info.keys():
        print(f"ERROR : Association of lipid {resname} and atom {atom_name} not found in parameter file")
        sys.exit()
    return dict_info[key]

def del_key_dict(dict_info, list_key):
    """
    Remove in a dictionary keys from list_key.

    --------------------
    INPUT
    dict_info: dictionary
        Contains information on the label clusters
    list_key: list
        Contains the labels on the edge of the matrix

    --------------------
    OUTPUT
    dictionary
        Contains information on the label clusters minus the edges ones
    """
    if len(list_key) != 0:
        for key in list_key:
            del dict_info[key]
    return dict_info

# return the maximum value from dico[key]=list
def max_value_dico(dico):
    return max([max(data2) for data2 in [dico[data] for data in dico.keys()]])

# return the minimum value from dico[key]=list
def min_value_dico(dico):
    return min([min(data2) for data2 in [dico[data] for data in dico.keys()]])

# from dico[liste] determine max length
def determine_lenMax2(dicoTot):
    return max([len(dicoTot[key]) for key in dicoTot])
