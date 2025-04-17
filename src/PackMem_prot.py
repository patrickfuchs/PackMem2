#!/usr/bin/python3
#-*- coding: utf-8 -*-
# Pg compute packing defects on membranes simulations
# R. Gautier A. Bacle april 2016
# R. Gerard may 2024
# M. Zygadlo august 2024

import sys
import argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder 

from bin import listes as l
from bin import matrix as m
from bin import pdb as pdb
from bin import connected_component as cc
from bin import dico as d
from bin import BasicFunctions as bfrg
from bin import param as p
from bin import protdist as pdist

def create_listZ(residues, list_resids, atom_mb, dist_suppl_Z, z_extr, up=True):
    """
    Create the listZ array for upper and lower leaflet
    ---------------------------------------------------------------------------
    INPUT:
    residues : MDAnalysis ResidueGroup
        Contains the names of all the residues selected
    list_resids : numpy array
        Contains all the residue numbers selected
    atom_mb : str
        Name of the reference atom for the lipids (C2)
    dist_suppl_Z : float
        Supplementary distance from the z coord
    z_extr : float
        Either maximum or minimum value of z coord
    up : boolean
        If we are on the upper leaflet
    ---------------------------------------------------------------------------
    OUTPUT:
    numpy array
        Contains floats ranging from z_coord to zmin
        or from zmax to z_coord by 1.0 steps
    """
    leaflet_listZ = {}
    for resid in list_resids:
        # Get the residue at a given residue number (resid)
        residue = residues[residues.resids == resid][0]
        # Get the atom given in the parameter file
        atom = residue.atoms[residue.atoms.names == atom_mb][0]
        # Get the z of this atom
        z_coord = atom.position[2]
        if up:
            tmp = l.create_array(round(z_coord - dist_suppl_Z, 2),
                                        round(z_extr +1.0, 2), m.SIZE)
        else:
            tmp = l.create_array(round(z_coord + dist_suppl_Z, 2),
                                            round(z_extr - 1.0, 2), -m.SIZE)
        # Reverse it
        tmp = np.flip(tmp)
        leaflet_listZ[resid] = tmp
    return leaflet_listZ

##########################################################################################
##### main
if __name__ == '__main__':
    
    """
    Python
    Script to compute Packing defect in flat bilayers, Use pdb file (from for example editconf) 
    be careful to PBC to create pdb file
    Lipids parameters adapted to Berger lipid FF (with corrections). Be careful to the atoms name if you use other lipids
    Lipids parameters from CHarmm36 FF (with Klauda corrections). Be careful to the atoms name if you use other lipids
    Lipids parameters from Martini FF. Be careful to the atoms name if you use other lipids
    fileRadius.txt example:
    DOPC  C02 1.875 a (aliphatic)
    DOPC  O8  1.48 n (non aliphatic)
    DOPC  C25 1.98 n (non aliphatic)
    Output files:
    outputname_Up/Lo_Shallow/Deep/All_result.txt
    ## MatrixSize  9646  9801         # Membrane matrix size, Total matrix size
    ## Total   51   582 11.41 6.034   # number of packing defects, total area of packing defects, average size, pourcent of membrane (Membrane matrix size)
    1    1    77.00     2.00    # for each packing defect num size x_position (first pixel) y_mean (first pixel)
    2    1    21.00     3.00 
    3    8    59.12     7.00 

    outputname_TotalUp/Lo_Shallow/Deep/All.pdb # PDB format only if verbose
    Matrix x,y with, for each cell, the value of "packing defect" (in B_factor column, the last column) 
    if 0 = Deep packing defects The value increases with the number of atoms in the cell.
    if >0 and < 1 = Shallow defect
    if -1 = Edge
    ATOM      3   H1 EDG     3       4.000  54.000  56.940  1.00 -1.00
    ATOM      3   H1 MAT     3       4.000  54.000  56.940  1.00  0.00
    ATOM      3   H1 MAT     3       4.000  54.000  56.940  1.00  0.00

    outputname+"_DefectUp/Lo_Shallow/Deep/All.pdb # PDB format only if verbose
    Packing defects in pdb format:
    the residue number corresponds to the different packing defects.
    ATOM      5   H1 DEF     2       6.000  22.000  56.940  1.00  2.00
    ATOM      5   H1 DEF     3       6.000  85.000  56.940  1.00  3.00
    ATOM      5   H1 DEF     3       6.000  86.000  56.940  1.00  3.00
    ATOM      6   H1 DEF     3       7.000  86.000  56.940  1.00  3.00
    """
    #######PARAMETRES et INPUT#####
    
    try:
        outputname = "output"
        pdbout = 0
        parser = argparse.ArgumentParser(description = 'Arguments for the app')
        parser.add_argument('-f', action = 'store', dest = 'traj', help = 'trajectory file (.xtc)')
        parser.add_argument('-s', action='store', dest='topo', help='Topology file (.gro)')
        parser.add_argument('-l', action='store', dest='lipid', help='Lipid name in the .gro file')
        parser.add_argument('-r', action = 'store', dest = 'filesrad',
                            help = 'File for radius (default vdw_radiiFinal2014.txt)',
                            default = 'vdw_radiiFinal2014.txt')
        parser.add_argument('-o', action = 'store', dest = 'outputname',
                            default = 'output', help = 'Name for output file (default output)')
        parser.add_argument('-d', action = 'store', dest = 'dist_suppl_Z', default = 1.0, type = float,
                            help = '(default value 1.0)')
        parser.add_argument('-t', action = 'store', dest = 'pd_type', default = 'all',
                            help = 'Packing Defect Type (all/shallow/deep) default all')
        parser.add_argument('-p', action = 'store', dest = 'paramFile', default = 'param.txt',
                            help = 'File for lipid parameters')
        parser.add_argument('-n', action = 'store', dest = 'indexFile',
                            help = 'File for index (Gromacs ndx style). Only Lower/Upper group accepted')
        parser.add_argument('-b', action = 'store', dest = 'start', type=int,
                            help = 'frame to start the analysis')
        parser.add_argument('-e', action = 'store', dest = 'end', type=int,
                            help = 'frame to end the analysis')
        parser.add_argument('-pdb', dest = 'pdbout', action = 'store_true',
                            help = 'Get .pdb outputs of the packing defects')
        parser.add_argument('-prot', dest = 'protein', action = 'store_true',
                            help = 'Analyse the packing defects close/far of the protein')
                            
        args = parser.parse_args()
        
        if len(sys.argv) == 1:
            parser.print_usage()
            sys.exit()
        
        # DICT_3L = {'DOP_N': 'DOP', 'DOE_N': 'DOE', 'DPP_N': 'DPP', etc}
        # RESNAME_GLYC = {'DOP': 'C2', 'DOE': 'C2', 'DPP': 'C2', etc}
        # LIPID = ['DOP', 'DOE', 'DPP', etc]
        DICT_3L, RESNAME_GLYC, LIPID = p.set_params(args.paramFile)

        # Determine packing defedct type (all/shallow/deep)
        if args.pd_type == 'all':
            FlagPDtype = 0
        elif args.pd_type == 'deep' :
            FlagPDtype = 1
        elif args.pd_type == 'shallow':
            FlagPDtype = 2
        else: 
            print('ERROR : Packing defect type not known. Please choose : all/deep/shallow')
            sys.exit()
            
        if args.dist_suppl_Z < 0.0 :
            print('ERROR : The distance for the -d option must be > 0')
            sys.exit()
        
        # ['REMARK    GENERATED BY TRJCONV\n',
        #  'TITLE     Title t= 100.00000 step= 50000\n',
        #  'REMARK    THIS IS A SIMULATION BOX\n',
        #  'CRYST1   56.638   56.638   74.530  90.00  90.00  90.00 P 1           1\n',
        #  'MODEL        2\n', 
        # 'ATOM      1  N   ARG A  43      28.420  41.300  49.850  1.00  0.00           N\n', etc]
        #pdblines = bfrg.read_file(args.filename) # remplace with MDAnalysis
        # {'ALA N': 1.85, 'ALA HN': 0.22, 'ALA HT1': 0.22, etc}
        # for the amino acids then the lipids
        radius = bfrg.read_radius(args.filesrad)
        # {'ALA N': 'n', 'ALA HN': 'n', 'ALA HT1': 'n', etc }
        # for the amino acids then the lipids
        aliphatic = bfrg.read_aliphatic(args.filesrad)
        
    except:
        print('Command line: PackMem.py -i file.pdb -r fileRadius.txt \
               -p param.txt -o output -d distGlyc -t deep/all/shallow [-n index file] [-v] or -h for help')
        sys.exit()

    ######### LOAD UNIVERSE  #########
    u = mda.Universe(args.topo, args.traj)
    # select the lipids in system
    system = u.select_atoms(f"resname {args.lipid}")
    
    ######### Extract UPPER/LOWER leaflet ##########
    atom_mb = RESNAME_GLYC[args.lipid]
    # If no index file seperating the leaflets
    if args.indexFile == None:
        # Create lists of the residue number  for upper and lower leaflets
        L = LeafletFinder(system, f'name {atom_mb}')
        upper_leaflet = L.groups(0).resids
        lower_leaflet = L.groups(1).resids
    else:
        lower_leaflet, upper_leaflet = p.read_ndx(args.indexFile)

    for ts in u.trajectory[args.start:args.end+1]:
        print(f"Working on frame {ts.frame:3d}")
        # select all atoms in systems
        #system = u.atoms
        system = u.select_atoms(f"resname {args.lipid} or protein")

        # Get all the residues in the membrane
        res_ids = system.resids
        dicoMb = list(set(res_ids))

        # If there are a lot of lipids (> 9999)
        if len(dicoMb) > 9999:
            print("The number of lipids in your membrane is very high (>9999)!")
            print("Please check that you only give to PackMem your membrane (+ protein)")
            sys.exit()

        # Get all x, y and z
        x_atoms = system.positions[:,0].round(2)
        y_atoms = system.positions[:,1].round(2)
        z_atoms = system.positions[:,2].round(2)
        # Get membrane dimension
        xmin, xmax, xmean = l.min_max(x_atoms)
        ymin, ymax, ymean = l.min_max(y_atoms)
        zmin, zmax, zmean = l.min_max(z_atoms)
        
        upper_listZ = create_listZ(system.residues, upper_leaflet, atom_mb, args.dist_suppl_Z, zmax)
        lower_listZ = create_listZ(system.residues, lower_leaflet, atom_mb, args.dist_suppl_Z, zmin, up=False)
                    
        # Build a lists from xmin-1 to xmax+1 every 1.0
        listX = l.create_array(int(xmin-1), int(xmax+2), m.SIZE)
        # Build a lists from ymin-1 to ymax+1 every 1.0
        listY = l.create_array(int(ymin-1), int(ymax+2), m.SIZE)

        # Initialize 2 matrixes for Upper and Lower leaflet of length listX,listY
        MatrixUp = m.initialize_matrix2D(len(listX), len(listY), np.nan)
        MatrixLo = m.initialize_matrix2D(len(listX), len(listY), np.nan)
                    
        ####################  Compute Matrix    #################
        # Search v cells around coord
        v = 5
        # For each atoms of lipids
        for i in range(len(res_ids)):
            atom_name = system.names[i]
            res_name = system.resnames[i]
            res_id = system.resids[i]
            radius_res=m.get_radius(radius, res_name, atom_name)
            # Get if the residue is aliphatic or not
            aliph_atom=m.get_aliphatic(aliphatic, res_name, atom_name)
            # Get the coordinates X Y Z
            coordtmp = [x_atoms[i], y_atoms[i], z_atoms[i]]
            # Upper leaflet ########################
            if res_id in upper_leaflet :
                # dfZ = z_C2_res - z_atom
                dfZ = round(m.find_Z(coordtmp[2], upper_listZ[res_id]),2)
                # If dfZ < 5
                # To limit the search around an atom to 5 cells
                if dfZ < (1. * v):
                    # Fill the matrix with value 0 < a < 1 for aliphatic
                    # Or with > 1 if polar OR deep
                    # Defects = 0
                    MatrixUp= m.fill_matrix(MatrixUp, coordtmp, listX, listY,
                                            upper_listZ[res_id], radius_res,
                                            FlagPDtype, aliph_atom)
            # Lower leaflet ########################
            if res_id in lower_leaflet :
                # dfZ = z_C2_res - z_atom
                dfZ = round(m.find_Z(coordtmp[2], lower_listZ[res_id]),2)
                # If dfZ > -5
                # To limit the search around an atom to 5 cells
                if dfZ > (-1. * v):
                    # Fill the matrix with value 0 < a < 1 for aliphatic
                    # Or with > 1 if polar OR deep
                    # Defects = 0
                    MatrixLo = m.fill_matrix(MatrixLo, coordtmp, listX, listY,
                                            lower_listZ[res_id], radius_res,
                                            FlagPDtype, aliph_atom)

        # Preliminary process for shallow defect and problem of the edges 
        # To eliminate shallow defects on edges first: binarize on all defects and storage edges coord
        # If shallow defects (2)
        if FlagPDtype == 2:
            # Initalise matrices to 0.0
            Matrix_UpbinM = m.initialize_matrix2D(len(listX),len(listY),0.)
            Matrix_LobinM = m.initialize_matrix2D(len(listX),len(listY),0.)
            # Binarise these matrices, with 0 for aliphatic atoms + packing defects and 1 otherwise
            Matrix_UpbinM = m.binarize_matrix_without0(MatrixUp, Matrix_UpbinM , -0.01, 0.99)
            Matrix_LobinM = m.binarize_matrix_without0(MatrixLo, Matrix_LobinM, -0.01, 0.99)
        
            # Get temporary packing defects
            # Connect the packing defects + label them + count the area
            Matrix_labels_UpM, root_labels_UpM, area_clusters_UpM, coor_clusters_UpM = \
                cc.get_connected_components(Matrix_UpbinM)
            Matrix_labels_LoM, root_labels_LoM, area_clusters_LoM, coor_clusters_LoM = \
                cc.get_connected_components(Matrix_LobinM)
            # Get the cluster on the edge
            clust_edge_UpM=cc.get_clusters_on_the_edge(Matrix_labels_UpM)
            clust_edge_LoM=cc.get_clusters_on_the_edge(Matrix_labels_LoM)

        # Binarisation matrices ###################################################
        # Initalise matrices to 0.0
        Matrix_Upbin=m.initialize_matrix2D(len(listX),len(listY),0.)
        Matrix_Lobin=m.initialize_matrix2D(len(listX),len(listY),0.)

        # If shallow (2)
        if FlagPDtype == 2:
            # Binarise these matrices, with 0 for aliphatic atoms and 1 otherwise
            Matrix_Upbin = m.binarize_matrix_without0(MatrixUp, Matrix_Upbin , 0, 0.99)
            Matrix_Lobin = m.binarize_matrix_without0(MatrixLo, Matrix_Lobin, 0, 0.99)
        # If all (0)
        elif FlagPDtype == 0:
            # Binarise these matrices, with 0 for aliphatic atoms + packing defects and 1 otherwise
            Matrix_Upbin = m.binarize_matrix_without0(MatrixUp, Matrix_Upbin , -0.01, 0.99)
            Matrix_Lobin = m.binarize_matrix_without0(MatrixLo, Matrix_Lobin, -0.01, 0.99)
        # If deep (1)
        else:
            # Binarise these matrices, with 0 if deep defect, 1 otherwise
            Matrix_Upbin = m.binarize_matrix(MatrixUp, Matrix_Upbin, 0.)
            Matrix_Lobin = m.binarize_matrix(MatrixLo, Matrix_Lobin, 0.)

        # If shallow defects (2)
        if FlagPDtype == 2:
            # Modify the binary matrix to take account edges (determined by all packing defects)
            # The clusters that had their labels on the edge are put to 0.0
            Matrix_Upbin = m.modify_matrix(Matrix_labels_UpM, Matrix_Upbin, clust_edge_UpM)
            Matrix_Lobin = m.modify_matrix(Matrix_labels_LoM, Matrix_Lobin, clust_edge_LoM)


        # Packing defects determination  ##########################################
        # Connect the packing defects + label them + count the area
        Matrix_labels_Up, root_labels_Up, area_clusters_Up, coor_clusters_Up = \
            cc.get_connected_components(Matrix_Upbin)
        Matrix_labels_Lo, root_labels_Lo, area_clusters_Lo, coor_clusters_Lo = \
            cc.get_connected_components(Matrix_Lobin)

        # Get cluster on the edge
        clust_edge_Up = cc.get_clusters_on_the_edge(Matrix_labels_Up)
        clust_edge_Lo = cc.get_clusters_on_the_edge(Matrix_labels_Lo)

        # Count area of the edge
        total_edge_Up = 0
        for key in clust_edge_Up:
            total_edge_Up += area_clusters_Up[key]
        total_edge_Lo = 0
        for key in clust_edge_Lo:
            total_edge_Lo += area_clusters_Lo[key]
        
        # Clean dico defects (without edge)
        area_clusters_Up = d.del_key_dico(area_clusters_Up, clust_edge_Up)
        coor_clusters_Up = d.del_key_dico(coor_clusters_Up, clust_edge_Up)
        area_clusters_Lo = d.del_key_dico(area_clusters_Lo, clust_edge_Lo)
        coor_clusters_Lo = d.del_key_dico(coor_clusters_Lo, clust_edge_Lo)

        # If shallow (2)
        if FlagPDtype == 2:
            # Eliminate nan inside (deep not shallow defect)
            Matrix_labels_Up, total_edge_Up, clustPb_Up = \
                m.clean_NA_inside(Matrix_labels_Up, clust_edge_Up, MatrixUp, total_edge_Up)
            root_labels_Up, area_clusters_Up, coor_clusters_Up = \
                cc.delete_NApoints_inside(clustPb_Up, Matrix_labels_Up,
                                        root_labels_Up, area_clusters_Up)

            Matrix_labels_Lo, total_edge_Lo, clustPb_Lo = \
                m.clean_NA_inside(Matrix_labels_Lo, clust_edge_Lo, MatrixLo, total_edge_Lo)
            root_labels_Lo, area_clusters_Lo, coor_clusters_Lo = \
                cc.delete_NApoints_inside(clustPb_Lo, Matrix_labels_Lo,
                                        root_labels_Lo, area_clusters_Lo)
            
            # Reclean dico defects (without edge)
            area_clusters_Up = d.del_key_dico(area_clusters_Up, clust_edge_Up)
            coor_clusters_Up = d.del_key_dico(coor_clusters_Up, clust_edge_Up)
            area_clusters_Lo = d.del_key_dico(area_clusters_Lo, clust_edge_Lo)
            coor_clusters_Lo = d.del_key_dico(coor_clusters_Lo, clust_edge_Lo)
        
        # Output text file  #######################################################
        total_size = len(listX) * len(listY)
        pdb.outputTXT_defects(f"{args.outputname}{ts.frame}", FlagPDtype, "Up", area_clusters_Up, 
                            coor_clusters_Up, total_size, total_edge_Up, listX, listY)
        pdb.outputTXT_defects(f"{args.outputname}{ts.frame}", FlagPDtype, "Lo", area_clusters_Lo, 
                            coor_clusters_Lo, total_size, total_edge_Lo, listX, listY)

        # Output PDB files ########################################################
        # final matrix values PD (X,Y) with Z cooresponding to the max(Upper/lowerZlevel)
        if args.pdbout :
            valzmax=float(d.max_value_dico(upper_listZ))
            valzmin=float(d.min_value_dico(lower_listZ))
            pdb.outputPDB_leaflet(system.resids, system.resnames, system.names,
                                  system.positions, upper_leaflet, 
                                  f"{args.outputname}{ts.frame}" + "_Upper.pdb", 
                                  ts.frame)
            pdb.outputPDB_leaflet(system.resids, system.names, system.resnames,
                                  system.positions, lower_leaflet, 
                                  f"{args.outputname}{ts.frame}" + "_Lower.pdb", 
                                  ts.frame)
                            
            pdb.outputPDB_Total_matrix(f"{args.outputname}{ts.frame}", FlagPDtype, "Up", ts.frame,
                                    listX, listY, valzmax, MatrixUp)
            pdb.outputPDB_Total_matrix(f"{args.outputname}{ts.frame}", FlagPDtype, "Lo", ts.frame,
                                    listX, listY, valzmin, MatrixLo)
                                    
            pdb.outputPDB_defects(f"{args.outputname}{ts.frame}", FlagPDtype, "Up", ts.frame,
                                listX, listY, valzmax, Matrix_labels_Up, clust_edge_Up)
            pdb.outputPDB_defects(f"{args.outputname}{ts.frame}", FlagPDtype, "Lo", ts.frame,
                                listX, listY, valzmin, Matrix_labels_Lo, clust_edge_Lo)
        
        # Distance from the protein ###############################################
        if args.protein :
            # Label 1 corresponds to the void around the simulation, so we remove it.
            Matrix_labels_Up = np.where(Matrix_labels_Up == 1, 0, Matrix_labels_Up)
            Matrix_labels_Lo = np.where(Matrix_labels_Lo == 1, 0, Matrix_labels_Lo)
            
            # Retrieve the protein    
            list_code3L_AA = ['ALA', 'ASP', 'ARG', 'ASN', 'CYS',
                                'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                                'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                                'SER', 'THR', 'TRP', 'TYR', 'VAL']

            # Create empty arrays
            array2d_prot_Up = np.zeros_like(Matrix_labels_Up)
            array2d_prot_Lo = np.zeros_like(Matrix_labels_Lo)

            # Find where the protein is in the simulation box
            for i in range(len(res_ids)) :
                res_name = system.resnames[i]
                # If it is an amino acid atom
                if res_name in list_code3L_AA:
                    atom_name = system.names[i]
                    # Get the coordinates X Y Z
                    coordtmp = [x_atoms[i], y_atoms[i], z_atoms[i]]
                    iX,iY = m.find_X_Y(coordtmp, listX, listY)
                    # Upper leaflet
                    if coordtmp[2] > zmean:
                        array2d_prot_Up[iX, iY] = 1
                    # Lower leaflet
                    if coordtmp[2] < zmean:
                        array2d_prot_Lo[iX, iY] = 1
            
            # Classification of Packing Defects by distance group 
            # Get the coordinates of the matrix where the edges of the packing defects are located.
            # dictionnary label coords {lab1 = [(x1, y1), (x2, y2), ...],
            #                           lab2 = [(x1, y1), (x2, y2), ...], 
            #                           ...                              }
            dict_labels_coor_Up = pdist.find_pd_border(Matrix_labels_Up)
            #dict_labels_coor_Lo = pdist.find_pd_border(Matrix_labels_Lo)

            # Get the coordinates of the matrix where the edges of the protein are located.
            # list of tuples [(x1, y1), (x2, y2), ...]
            list_edge_coor_prot_Up = pdist.find_prot_border(array2d_prot_Up)
            #list_edge_coor_prot_Lo = pdist.find_prot_border(array2d_prot_Lo)

            # Assign distance group for each packing defect, "far" or "close". Default threshold = 10 A.
            # dict {lab1 : 'group', lab2 : 'group', ... }
            pd_labels_group_Up = pdist.assign_dist_group(list_edge_coor_prot_Up, dict_labels_coor_Up, 10)
            #pd_labels_group_Lo = pdist.assign_dist_group(list_edge_coor_prot_Lo, dict_labels_coor_Lo, 10)

            # Results
            # Ouput text file
            # header : label,dist_group,area
            pdist.outputTXT_defects_prot(f"{args.outputname}{ts.frame}", FlagPDtype, "Up", pd_labels_group_Up, area_clusters_Up)
            #pdist.outputTXT_defects_prot(f"{args.outputname}{ts.frame}", FlagPDtype, "Lo", pd_labels_group_Lo, area_clusters_Lo)
