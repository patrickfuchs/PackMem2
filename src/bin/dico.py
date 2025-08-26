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

def max_value_dict(dict_arrayZ):
    """
    Return the maximum value from the dictionary.

    --------------------
    INPUT
    dict_arrayZ:  dictionary
        For each resid, contains array of floats ranging from z_coord to zmin
        or from zmax to z_coord by 1.0 steps
    
    --------------------
    OUTPUT
    float
        The maximum value from the dictionary
    """
    return max([max(arrayZ) for arrayZ in dict_arrayZ.values()])

def min_value_dict(dict_arrayZ):
    """
    Return the minimum value from the dictionary.

    --------------------
    INPUT
    dict_arrayZ:  dictionary
        For each resid, contains array of floats ranging from z_coord to zmin
        or from zmax to z_coord by 1.0 steps
    
    --------------------
    OUTPUT
    float
        The minimum value from the dictionary
    """
    return min([min(arrayZ) for arrayZ in dict_arrayZ.values()])
