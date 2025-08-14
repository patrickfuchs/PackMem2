#-*- coding: utf-8 -*-
# R. Gautier A. Bacle 2015
# M. Zygadlo 2025

import argparse
import os
import numpy as np

def file_present(filename):
    """
    Check if the file exists.

    --------------------
    INPUT
    filename: str
        The name if the file to test the presence of
    """
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"ERROR : file '{filename}' not found.")

def get_args():
    """
    Get the arguments for the script and check that the inputfiles are valid.

    --------------------
    OUTPUT
    parser.parse_args
        Contains all the arguments for the script
    """
    parser = argparse.ArgumentParser(description = 'Arguments for PackMem_prot.py')
    parser.add_argument('-f', action='store', dest='traj', required=True,
                        help = 'Trajectory file (.xtc)')
    parser.add_argument('-s', action='store', dest='topo', required=True,
                        help = 'Topology file (.gro)')
    parser.add_argument('-l', action='store', dest='lipid', required=True,
                        help = 'Lipid name in the .gro file')
    parser.add_argument('-b', action='store', dest='start', type=int,
                        default = 0,
                        help = 'Frame to start the analysis (default: 0)')
    parser.add_argument('-e', action='store', dest='end', type=int,
                        default = None,
                        help = 'Frame to end the analysis (default: None)')
    parser.add_argument('-r', action='store', dest='filesrad',
                        default = 'vdw_radii_Charmm.txt',
                        help = 'File for the atom radius (default: vdw_radii_Charmm.txt)')
    parser.add_argument('-p', action='store', dest='paramFile',
                        default = 'param_Charmm.txt',
                        help = 'File for lipid parameters (default: param_Charmm.txt)')
    parser.add_argument('-o', action='store', dest='outputname',
                        default = 'output',
                        help = 'Name for output file (default: output)')
    parser.add_argument('-d', action='store', dest='dist_suppl_Z', type=float,
                        default = 1.0, 
                        help = 'Distance to differenciate Deep from Shallow defects (default: 1.0)')
    parser.add_argument('-n', action='store', dest='indexFile',
                        help = 'Index file (Gromacs ndx style) of only Lower/Upper group')
    parser.add_argument('-pdb', dest = 'pdbout', action = 'store_true',
                        help = 'Get .pdb outputs of the packing defects')
    parser.add_argument('-prot', dest = 'protein', action = 'store_true',
                        help = 'Analyse the packing defects close/far of the protein')
                        
    args = parser.parse_args()

    # Check that the files exist
    file_present(args.traj)
    file_present(args.topo)
    file_present(args.paramFile)
    file_present(args.filesrad)    

    # Check that the variable d is positive
    if args.dist_suppl_Z < 0.0 :
        raise Exception("ERROR : The distance for the -d option must be > 0.0")
    
    return args

# load file return list[line(string)]
def read_file(filename):
    try: 
        with open(filename) as f:
            data = f.readlines()
    except : 
        print("ERROR : Something went wrong with the file %s" % filename)
    return data

def dict_2columns(list_str):
    """
    Transform list of strings with 2 columns into a dictionary.

    --------------------
    INPUT
    list_str: list
        Contains strings as values
    
    --------------------
    OUTPUT
    dictionary
        Contains for each string the key as first word and value as second word
    """
    dic = {}
    for line in list_str:
        # Use space a separator
        data = line.strip().split(' ')
        dic[data[0]] = data[1]
    return dic

def set_params(filename):
    """
    Read the parameter file and return a dictionary of the information.

    --------------------
    INPUT
    filename: str
        Name of the parameter file
    
    --------------------
    OUTPUT
    dictionary
        Contains the name of the lipids as key and their central atoms as value
    """
    # Read the file
    lines = read_file(filename)
    # Create a dictionary from the file
    resname_glyc = dict_2columns(lines)
    return resname_glyc

# open VdW radius -> dict["RES(3L) ATOMname"]=radius in A 
def read_radius(filename):
    dico={}
    data = read_file(filename)
    for i in range(0, len(data)):
        data[i]=data[i].strip()
        if data[i] != "":
            data[i]=data[i].split()
            cle=data[i][0]+" "+data[i][1]
            dico[cle]=float(data[i][2])
    return dico

#new 2016-08-18
# open aliphatic table-> dict["RES(3L) ATOMname"]= p (polar) a (aliphatic) 
def read_aliphatic(filename):
    dico={}
    data = read_file(filename)
    for i in range(0, len(data)):
        data[i]=data[i].strip()
        if data[i] != "":
            data[i]=data[i].split()
            cle=data[i][0]+" "+data[i][1]
            dico[cle]=data[i][3]
    return dico

def read_ndx(ndx_file):
    """
    Read a .ndx file to set the lower/upper residue number lists.

    --------------------
    INPUT
    ndx_file: str
        Name if the index file
    --------------------
    OUTPUT
    array numpy
        Contains the residue number of each lipid in the upper/lower leaflet
    """
    # Read the lines in the index file
    lines = bfrg.read_file(ndx_file)
    # Get the 2nd and last last into a numpy array
    # Because these lines are where the residue numbers are
    index1 = np.array(lines[1].strip().split(' '), dtype=int)
    index2 = np.array(lines[3].strip().split(' '), dtype=int)

    # If the first line : "[ Upper leaflet ]" or "[ upper leaflet ]"
    if "upper" in lines[0] or "Upper" in lines[0]:
        return index1, index2
    # If the first  line : "[ Lower leaflet ]" or "[ lower leaflet ]"
    else:
        return index2, index1
