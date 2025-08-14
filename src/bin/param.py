#-*- coding: utf-8 -*-
# R. Gautier  A. Bacle 2015

import argparse
import os
import re
from bin import BasicFunctions as bfrg

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

# transform 2 columns of file to a dictionary
def dict_2columns(tab1):
    lipid_3L  = {}
    for i in tab1:
        data = i.strip().split(' ')
        lipid_3L[data[0]] = data[1]
    return lipid_3L

#RG 2016 01 11
# transform  a lipid dictionary [LIPIDFILE_ATOMFILE]=LIPID3L to a list of LIPIDFILE
def dict_lipid(lipid_3L):
    lipid  = []
    for key in lipid_3L:
        data = key.strip()
        lipid.append(data)
    return lipid

# read the parameters file and return 3 tables
#lis le fichier de paramètres en renvoir 3 tableaux
#tab2 : Glycerol atoms by lipid type
#tab3 : Dictionary of aliphatic atoms for each lipid type
def set_params(filename):
    lines = bfrg.read_file(filename)
    resname_glyc = dict_2columns(lines[103:])
    lipid=dict_lipid(resname_glyc)
    return(resname_glyc, lipid)

#return the lines number for each lipid
def limits_lip(tab):
    list_limit_lip = []
    for i, lip in enumerate(tab):
        line = len(lip.strip().split(' '))
        if line > 1:
            list_limit_lip.append(i)
    return list_limit_lip

#build dictionary with key= lipd name, values = aliphatic atoms list
def dic_aliph_atoms(tab, list_limits):
    dic_aliph_atoms = {}
    for limit in list_limits :
        name_lip = tab[limit].strip().split(' ')[0]
        nb_aliph = int(tab[limit].strip().split(' ')[1])
        tab_aliph = []
        for i in range(nb_aliph):
            tab_aliph.append(tab[limit+i+1].strip())
        dic_aliph_atoms[name_lip] = tab_aliph
    return(dic_aliph_atoms)

#Read a ndx file to set the lower/upper residue number lists
def read_ndx(ndx_file):
    lines = bfrg.read_file(ndx_file)
    lower_leaflet = []
    upper_leaflet = []
    index1 = lines[1].split(' ')
    index2 = lines[3].split(' ')
    for number in index1:
        res_num = int(number.strip())
        if "Lower" in lines[0] or "lower" in lines[0]:
            lower_leaflet.append(res_num)
        else: 
            upper_leaflet.append(res_num)
    for number in index2:
        res_num = int(number.strip())
        if "Lower" in lines[0] or "lower" in lines[0]:
            upper_leaflet.append(res_num)
        else:
            lower_leaflet.append(res_num)
    return(lower_leaflet, upper_leaflet)
    







