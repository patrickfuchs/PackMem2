#!/usr/bin/python3
#-*- coding: utf-8 -*-
# Pg compute packing defects on membranes simulations
# R. Gautier A. Bacle april 2016
# R. Gerard may 2024
# M. Zygadlo august 2024

import sys
import argparse
import warnings
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
    ####### PARAMETRES et INPUT #####
    try:
        outputname = "output"
        pdbout = 0
        parser = argparse.ArgumentParser(description = 'Arguments for the app')
        parser.add_argument('-f', action='store', dest='traj',
                            help = 'Trajectory file (.xtc)')
        parser.add_argument('-s', action='store', dest='topo',
                            help = 'Topology file (.gro)')
        parser.add_argument('-l', action='store', dest='lipid',
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
        
        if len(sys.argv) == 1:
            parser.print_usage()
            sys.exit()
            
        if args.dist_suppl_Z < 0.0 :
            print('ERROR : The distance for the -d option must be > 0')
            sys.exit()
        
    except:
        print('Command line: PackMem_prot.py -f file.xtc -s file.gro -l lipid_name\
              -b start_frame -e end_frame -r radius.txt -p param.txt -o output\
              -d distGlyc [-n index file] [-pdb] [-prot] or -h for help')
        sys.exit()
    

    ######### READ PARAM FILES #########
    # DICT_3L = {'DOP_N': 'DOP', 'DOE_N': 'DOE', 'DPP_N': 'DPP', etc}
    # RESNAME_GLYC = {'DOP': 'C2', 'DOE': 'C2', 'DPP': 'C2', etc}
    # LIPID = ['DOP', 'DOE', 'DPP', etc]
    DICT_3L, RESNAME_GLYC, LIPID = p.set_params(args.paramFile)
    # {'ALA N': 1.85, 'ALA HN': 0.22, 'ALA HT1': 0.22, etc}
    # for the amino acids then the lipids
    radius = bfrg.read_radius(args.filesrad)
    # {'ALA N': 'n', 'ALA HN': 'n', 'ALA HT1': 'n', etc }
    # for the amino acids then the lipids
    aliphatic = bfrg.read_aliphatic(args.filesrad)


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


        ####################  Compute Matrix    #################
        # Initialize 2 matrixes for Upper and Lower leaflet of length listX,listY
        MatrixUp_Deep = m.initialize_matrix2D(len(listX), len(listY), np.nan)
        MatrixLo_Deep = m.initialize_matrix2D(len(listX), len(listY), np.nan)

        MatrixUp_Shallow = m.initialize_matrix2D(len(listX), len(listY), np.nan)
        MatrixLo_Shallow = m.initialize_matrix2D(len(listX), len(listY), np.nan)

        MatrixUp_All = m.initialize_matrix2D(len(listX), len(listY), np.nan)
        MatrixLo_All = m.initialize_matrix2D(len(listX), len(listY), np.nan)
                    
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
                dfZ = round(m.diff_Z(upper_listZ[res_id], coordtmp[2]),2)
                # If dfZ < 5
                # To limit the search around an atom to 5 cells
                if dfZ < (1. * v):
                    # Fill the matrix with value 0 < a < 1 for aliphatic
                    # Or with > 1 if polar OR deep
                    # Defects = 0
                    MatrixUp_Deep = m.fill_matrix(MatrixUp_Deep, coordtmp, listX, listY,
                                            upper_listZ[res_id], radius_res,
                                            "deep", aliph_atom)
                    MatrixUp_Shallow = m.fill_matrix(MatrixUp_Shallow, coordtmp, listX, listY,
                                            upper_listZ[res_id], radius_res,
                                            "shallow", aliph_atom)
                    MatrixUp_All = m.fill_matrix(MatrixUp_All, coordtmp, listX, listY,
                                            upper_listZ[res_id], radius_res,
                                            "all", aliph_atom)
            # Lower leaflet ########################
            if res_id in lower_leaflet :
                # dfZ = z_C2_res - z_atom
                dfZ = round(m.diff_Z(lower_listZ[res_id], coordtmp[2]),2)
                # If dfZ > -5
                # To limit the search around an atom to 5 cells
                if dfZ > (-1. * v):
                    # Fill the matrix with value 0 < a < 1 for aliphatic
                    # Or with > 1 if polar OR deep
                    # Defects = 0
                    MatrixLo_Deep = m.fill_matrix(MatrixLo_Deep, coordtmp, listX, listY,
                                            lower_listZ[res_id], radius_res,
                                            "deep", aliph_atom)
                    MatrixLo_Shallow = m.fill_matrix(MatrixLo_Shallow, coordtmp, listX, listY,
                                            lower_listZ[res_id], radius_res,
                                            "shallow", aliph_atom)
                    MatrixLo_All = m.fill_matrix(MatrixLo_All, coordtmp, listX, listY,
                                            lower_listZ[res_id], radius_res,
                                            "all", aliph_atom)


        ####################  Binarise the matrix    #################
        # Preliminary process for shallow defect and problem of the edges 
        # To eliminate shallow defects on edges first: binarize on all defects and storage edges coord
        # Initalise matrices to 0.0
        MatrixUp_Shallowbin = m.initialize_matrix2D(len(listX),len(listY),0.)
        MatrixLo_Shallowbin = m.initialize_matrix2D(len(listX),len(listY),0.)
        # Binarise these matrices, with 0 for aliphatic atoms + packing defects and 1 otherwise
        MatrixUp_binM = m.binarize_matrix_without0(MatrixUp_Shallow, MatrixUp_Shallowbin , -0.01, 0.99)
        MatrixLo_binM = m.binarize_matrix_without0(MatrixLo_Shallow, MatrixLo_Shallowbin, -0.01, 0.99)
    
        # Get temporary packing defects
        # Connect the packing defects + label them + count the area
        MatrixUp_labels_Shallow, rootUp_labels_Shallow, areaUp_clusters_Shallow, coorUp_clusters_Shallow = \
            cc.get_connected_components(MatrixUp_Shallowbin)
        MatrixLo_labels_Shallow, rootLo_labels_Shallow, areaLo_clusters_Shallow, coorLo_clusters_Shallow = \
            cc.get_connected_components(MatrixLo_Shallowbin)
        # Get the cluster on the edge
        labelsUp_edge_Shallow=cc.get_clusters_on_the_edge(MatrixUp_labels_Shallow)
        labelsLo_edge_Shallow=cc.get_clusters_on_the_edge(MatrixLo_labels_Shallow)


        # Initalise matrices to 0.0
        MatrixUp_Deepbin=m.initialize_matrix2D(len(listX),len(listY),0.)
        MatrixLo_Deepbin=m.initialize_matrix2D(len(listX),len(listY),0.)
        MatrixUp_Shallowbin=m.initialize_matrix2D(len(listX),len(listY),0.)
        MatrixLo_Shallowbin=m.initialize_matrix2D(len(listX),len(listY),0.)
        MatrixUp_Allbin=m.initialize_matrix2D(len(listX),len(listY),0.)
        MatrixLo_Allbin=m.initialize_matrix2D(len(listX),len(listY),0.)

        #### Shallow ####
        # Binarise these matrices, with 0 for aliphatic atoms and 1 otherwise
        MatrixUp_Shallowbin = m.binarize_matrix_without0(MatrixUp_Shallow, MatrixUp_Shallowbin , 0, 0.99)
        MatrixLo_Shallowbin = m.binarize_matrix_without0(MatrixLo_Shallow, MatrixLo_Shallowbin, 0, 0.99)
        # Modify the binary matrix to take account edges (determined by all packing defects)
        # The clusters that had their labels on the edge are put to 0.0
        MatrixUp_Shallowbin = m.modify_matrix(MatrixUp_labels_Shallow, MatrixUp_Shallowbin, labelsUp_edge_Shallow)
        MatrixLo_Shallowbin = m.modify_matrix(MatrixLo_labels_Shallow, MatrixLo_Shallowbin, labelsLo_edge_Shallow)
        # Packing defects determination
        # Connect the packing defects + label them + count the area
        MatrixUp_labels_Shallow, rootUp_labels_Shallow, areaUp_clusters_Shallow, coorUp_clusters_Shallow = \
            cc.get_connected_components(MatrixUp_Shallowbin)
        MatrixLo_labels_Shallow, rootLo_labels_Shallow, areaLo_clusters_Shallow, coorLo_clusters_Shallow = \
            cc.get_connected_components(MatrixLo_Shallowbin)
        # Get cluster on the edge
        labelsUp_edge_Shallow = cc.get_clusters_on_the_edge(MatrixUp_labels_Shallow)
        labelsLo_edge_Shallow = cc.get_clusters_on_the_edge(MatrixLo_labels_Shallow)
        # Count area of the edge
        totalUp_edge_Shallow = m.count_edge_area(areaUp_clusters_Shallow, labelsUp_edge_Shallow)
        totalLo_edge_Shallow = m.count_edge_area(areaLo_clusters_Shallow, labelsLo_edge_Shallow)
        # Clean dico defects (without edge)
        areaUp_clusters_Shallow = d.del_key_dico(areaUp_clusters_Shallow, labelsUp_edge_Shallow)
        coorUp_clusters_Shallow = d.del_key_dico(coorUp_clusters_Shallow, labelsUp_edge_Shallow)
        areaLo_clusters_Shallow = d.del_key_dico(areaLo_clusters_Shallow, labelsLo_edge_Shallow)
        coorLo_clusters_Shallow = d.del_key_dico(coorLo_clusters_Shallow, labelsLo_edge_Shallow)
        # Eliminate nan inside (deep not shallow defect)
        MatrixUp_labels_Shallow, totalUp_edge_Shallow, clustPb_Up_Shallow = \
            m.clean_NA_inside(MatrixUp_labels_Shallow, labelsUp_edge_Shallow, MatrixUp_Shallow, totalUp_edge_Shallow)
        rootUp_labels_Shallow, areaUp_clusters_Shallow, coorUp_clusters_Shallow = \
            cc.delete_NApoints_inside(clustPb_Up_Shallow, MatrixUp_labels_Shallow,
                                    rootUp_labels_Shallow, areaUp_clusters_Shallow)

        MatrixLo_labels_Shallow, totalLo_edge_Shallow, clustPb_Lo_Shallow = \
            m.clean_NA_inside(MatrixLo_labels_Shallow, labelsLo_edge_Shallow, MatrixLo_Shallow, totalLo_edge_Shallow)
        rootLo_labels_Shallow, areaLo_clusters_Shallow, coorLo_clusters_Shallow = \
            cc.delete_NApoints_inside(clustPb_Lo_Shallow, MatrixLo_labels_Shallow,
                                    rootLo_labels_Shallow, areaLo_clusters_Shallow)
        # Reclean dico defects (without edge)
        areaUp_clusters_Shallow = d.del_key_dico(areaUp_clusters_Shallow, labelsUp_edge_Shallow)
        coorUp_clusters_Shallow = d.del_key_dico(coorUp_clusters_Shallow, labelsUp_edge_Shallow)
        areaLo_clusters_Shallow = d.del_key_dico(areaLo_clusters_Shallow, labelsLo_edge_Shallow)
        coorLo_clusters_Shallow = d.del_key_dico(coorLo_clusters_Shallow, labelsLo_edge_Shallow)


        #### All ####
        # Binarise these matrices, with 0 for aliphatic atoms + packing defects and 1 otherwise
        MatrixUp_Allbin = m.binarize_matrix_without0(MatrixUp_All, MatrixUp_Allbin , -0.01, 0.99)
        MatrixLo_Allbin = m.binarize_matrix_without0(MatrixLo_All, MatrixLo_Allbin, -0.01, 0.99)
        # Packing defects determination
        # Connect the packing defects + label them + count the area
        MatrixUp_labels_All, rootUp_labels_All, areaUp_clusters_All, coorUp_clusters_All = \
            cc.get_connected_components(MatrixUp_Allbin)
        MatrixLo_labels_All, rootLo_labels_All, areaLo_clusters_All, coorLo_clusters_All = \
            cc.get_connected_components(MatrixLo_Allbin)
        # Get cluster on the edge
        labelsUp_edge_All = cc.get_clusters_on_the_edge(MatrixUp_labels_All)
        labelsLo_edge_All = cc.get_clusters_on_the_edge(MatrixLo_labels_All)
        # Count area of the edge
        totalUp_edge_All = m.count_edge_area(areaUp_clusters_All, labelsUp_edge_All)
        totalLo_edge_All = m.count_edge_area(areaLo_clusters_All, labelsLo_edge_All)
        # Clean dico defects (without edge)
        areaUp_clusters_All = d.del_key_dico(areaUp_clusters_All, labelsUp_edge_All)
        coorUp_clusters_All = d.del_key_dico(coorUp_clusters_All, labelsUp_edge_All)
        areaLo_clusters_All = d.del_key_dico(areaLo_clusters_All, labelsLo_edge_All)
        coorLo_clusters_All = d.del_key_dico(coorLo_clusters_All, labelsLo_edge_All)


        #### Deep ####
        # Binarise these matrices, with 0 if deep defect, 1 otherwise
        MatrixUp_Deepbin = m.binarize_matrix_without0(MatrixUp_Deep, MatrixUp_Deepbin, -0.01, 0.001)
        MatrixLo_Deepbin = m.binarize_matrix_without0(MatrixLo_Deep, MatrixLo_Deepbin, -0.01, 0.001)
        # Packing defects determination
        # Connect the packing defects + label them + count the area
        MatrixUp_labels_Deep, rootUp_labels_Deep, areaUp_clusters_Deep, coorUp_clusters_Deep = \
            cc.get_connected_components(MatrixUp_Deepbin)
        MatrixLo_labels_Deep, rootLo_labels_Deep, areaLo_clusters_Deep, coorLo_clusters_Deep = \
            cc.get_connected_components(MatrixLo_Deepbin)
        # Get cluster on the edge
        labelsUp_edge_Deep = cc.get_clusters_on_the_edge(MatrixUp_labels_Deep)
        labelsLo_edge_Deep = cc.get_clusters_on_the_edge(MatrixLo_labels_Deep)
        # Count area of the edge
        totalUp_edge_Deep = m.count_edge_area(areaUp_clusters_Deep, labelsUp_edge_Deep)
        totalLo_edge_Deep = m.count_edge_area(areaLo_clusters_Deep, labelsLo_edge_Deep)
        # Clean dico defects (without edge)
        areaUp_clusters_Deep = d.del_key_dico(areaUp_clusters_Deep, labelsUp_edge_Deep)
        coorUp_clusters_Deep = d.del_key_dico(coorUp_clusters_Deep, labelsUp_edge_Deep)
        areaLo_clusters_Deep = d.del_key_dico(areaLo_clusters_Deep, labelsLo_edge_Deep)
        coorLo_clusters_Deep = d.del_key_dico(coorLo_clusters_Deep, labelsLo_edge_Deep)

        
        
        ####################  Output text file  #################
        # Compute the total area of the matrix
        total_size = len(listX) * len(listY)
        pdb.outputTXT_defects(f"{args.outputname}{ts.frame}", "deep", "Up", areaUp_clusters_Deep, 
                            coorUp_clusters_Deep, total_size, totalUp_edge_Deep, listX, listY)
        pdb.outputTXT_defects(f"{args.outputname}{ts.frame}", "deep", "Lo", areaLo_clusters_Deep, 
                            coorLo_clusters_Deep, total_size, totalLo_edge_Deep, listX, listY)
        
        pdb.outputTXT_defects(f"{args.outputname}{ts.frame}", "shallow", "Up", areaUp_clusters_Shallow, 
                            coorUp_clusters_Shallow, total_size, totalUp_edge_Shallow, listX, listY)
        pdb.outputTXT_defects(f"{args.outputname}{ts.frame}", "shallow", "Lo", areaLo_clusters_Shallow, 
                            coorLo_clusters_Shallow, total_size, totalLo_edge_Shallow, listX, listY)
        
        pdb.outputTXT_defects(f"{args.outputname}{ts.frame}", "all", "Up", areaUp_clusters_All, 
                            coorUp_clusters_All, total_size, totalUp_edge_All, listX, listY)
        pdb.outputTXT_defects(f"{args.outputname}{ts.frame}", "all", "Lo", areaLo_clusters_All, 
                            coorLo_clusters_All, total_size, totalLo_edge_All, listX, listY)


        ####################  Output PDB files  #################
        # final matrix values PD (X,Y) with Z cooresponding to the max(Upper/lowerZlevel)
        if args.pdbout :
            # Get the max/min of the z_coord+1
            valzmax=float(d.max_value_dico(upper_listZ))
            valzmin=float(d.min_value_dico(lower_listZ))

            # Write the leaflets into a PDB
            # To ignore the warnings when writing a PDB
            warnings.filterwarnings("ignore", category=UserWarning)
            # Write the leaflets into a PDB
            u.select_atoms(f"resid {upper_leaflet[0]}:{upper_leaflet[-1]}").write(f"{args.outputname}{ts.frame}_Upper_leaflet.pdb")
            u.select_atoms(f"resid {lower_leaflet[0]}:{lower_leaflet[-1]}").write(f"{args.outputname}{ts.frame}_Lower_leaflet.pdb")
            
            # Write the Matrix (each cell) into a PDB 
            pdb.outputPDB_Total_matrix(f"{args.outputname}{ts.frame}", "deep", "Up", ts.frame,
                                    listX, listY, valzmax, MatrixUp_Deep)
            pdb.outputPDB_Total_matrix(f"{args.outputname}{ts.frame}", "deep", "Lo", ts.frame,
                                    listX, listY, valzmin, MatrixLo_Deep)
            pdb.outputPDB_Total_matrix(f"{args.outputname}{ts.frame}", "shallow", "Up", ts.frame,
                                    listX, listY, valzmax, MatrixUp_Shallow)
            pdb.outputPDB_Total_matrix(f"{args.outputname}{ts.frame}", "shallow", "Lo", ts.frame,
                                    listX, listY, valzmin, MatrixLo_Shallow)
            pdb.outputPDB_Total_matrix(f"{args.outputname}{ts.frame}", "all", "Up", ts.frame,
                                    listX, listY, valzmax, MatrixUp_All)
            pdb.outputPDB_Total_matrix(f"{args.outputname}{ts.frame}", "all", "Lo", ts.frame,
                                    listX, listY, valzmin, MatrixLo_All)
            
            # Write the defects into a PDB
            pdb.outputPDB_defects(f"{args.outputname}{ts.frame}", "deep", "Up", ts.frame,
                                listX, listY, valzmax, MatrixUp_labels_Deep, labelsUp_edge_Deep)
            pdb.outputPDB_defects(f"{args.outputname}{ts.frame}", "deep", "Lo", ts.frame,
                                listX, listY, valzmin, MatrixLo_labels_Deep, labelsLo_edge_Deep)
            pdb.outputPDB_defects(f"{args.outputname}{ts.frame}", "shallow", "Up", ts.frame,
                                listX, listY, valzmax, MatrixUp_labels_Shallow, labelsUp_edge_Shallow)
            pdb.outputPDB_defects(f"{args.outputname}{ts.frame}", "shallow", "Lo", ts.frame,
                                listX, listY, valzmin, MatrixLo_labels_Shallow, labelsLo_edge_Shallow)
            pdb.outputPDB_defects(f"{args.outputname}{ts.frame}", "all", "Up", ts.frame,
                                listX, listY, valzmax, MatrixUp_labels_All, labelsUp_edge_All)
            pdb.outputPDB_defects(f"{args.outputname}{ts.frame}", "all", "Lo", ts.frame,
                                listX, listY, valzmin, MatrixLo_labels_All, labelsLo_edge_All)
        

        ####################  Distance from the protein  #################
        if args.protein :
            # Select only the protein
            protein = u.select_atoms("protein")
            pos_prot = protein.positions

            # Label 1 corresponds to the void around the simulation, so we remove it.
            MatrixUp_labels_Deep = np.where(MatrixUp_labels_Deep == 1, 0, MatrixUp_labels_Deep)
            MatrixLo_labels_Deep = np.where(MatrixLo_labels_Deep == 1, 0, MatrixLo_labels_Deep)
            MatrixUp_labels_Shallow = np.where(MatrixUp_labels_Shallow == 1, 0, MatrixUp_labels_Shallow)
            MatrixLo_labels_Shallow = np.where(MatrixLo_labels_Shallow == 1, 0, MatrixLo_labels_Shallow)
            MatrixUp_labels_All = np.where(MatrixUp_labels_All == 1, 0, MatrixUp_labels_All)
            MatrixLo_labels_All = np.where(MatrixLo_labels_All == 1, 0, MatrixLo_labels_All)

            # Create empty arrays
            MatrixUp_prot_Deep = m.initialize_matrix2D(len(listX), len(listY), 0)
            MatrixLo_prot_Deep = m.initialize_matrix2D(len(listX), len(listY), 0)
            MatrixUp_prot_Shallow = m.initialize_matrix2D(len(listX), len(listY), 0)
            MatrixLo_prot_Shallow = m.initialize_matrix2D(len(listX), len(listY), 0)
            MatrixUp_prot_All = m.initialize_matrix2D(len(listX), len(listY), 0)
            MatrixLo_prot_All = m.initialize_matrix2D(len(listX), len(listY), 0)

            # Find where the protein is in the matrix
            MatrixUp_prot_Deep = pdist.find_protein(MatrixUp_prot_Deep, 'up', protein, listX, listY, zmean)
            MatrixUp_prot_Shallow = pdist.find_protein(MatrixUp_prot_Shallow, 'up', protein, listX, listY, zmean)
            MatrixUp_prot_All = pdist.find_protein(MatrixUp_prot_All, 'up', protein, listX, listY, zmean)
            MatrixLo_prot_Deep = pdist.find_protein(MatrixLo_prot_Deep, 'lo', protein, listX, listY, zmean)
            MatrixLo_prot_Shallow = pdist.find_protein(MatrixLo_prot_Shallow, 'lo', protein, listX, listY, zmean)
            MatrixLo_prot_All = pdist.find_protein(MatrixLo_prot_All, 'lo', protein, listX, listY, zmean)
            

            # If there are proteins on the upper leaflet
            if len(np.argwhere(MatrixUp_prot_All > 0)) > 0:
                # Classification of Packing Defects by distance group 
                # Get the coordinates of the matrix where the edges of the packing defects are located.
                # dictionnary label coords {lab1 = [(x1, y1), (x2, y2), ...],
                #                           lab2 = [(x1, y1), (x2, y2), ...], 
                #                           ...                              }
                labelsUp_coor_Deep = pdist.find_pd_border(MatrixUp_labels_Deep)
                labelsUp_coor_Shallow = pdist.find_pd_border(MatrixUp_labels_Shallow)
                labelsUp_coor_All = pdist.find_pd_border(MatrixUp_labels_All)

                # Get the coordinates of the matrix where the edges of the protein are located.
                # list of tuples [(x1, y1), (x2, y2), ...]
                edgeUp_coor_prot_Deep = pdist.find_prot_border(MatrixUp_prot_Deep)
                edgeUp_coor_prot_Shallow = pdist.find_prot_border(MatrixUp_prot_Shallow)
                edgeUp_coor_prot_All = pdist.find_prot_border(MatrixUp_prot_All)

                # Assign distance group for each packing defect, "far" or "close". Default threshold = 10 A.
                # dict {lab1 : 'group', lab2 : 'group', ... }
                DefectsUp_labels_group_Deep = pdist.assign_dist_group(edgeUp_coor_prot_Deep, labelsUp_coor_Deep, 10)
                DefectsUp_labels_group_Shallow = pdist.assign_dist_group(edgeUp_coor_prot_Shallow, labelsUp_coor_Shallow, 10)
                DefectsUp_labels_group_All = pdist.assign_dist_group(edgeUp_coor_prot_All, labelsUp_coor_All, 10)

                # Write the result in a text file
                # format : label,dist_group,area
                pdist.outputTXT_defects_prot(f"Prot_{args.outputname}{ts.frame}", "deep", "Up", DefectsUp_labels_group_Deep, areaUp_clusters_Deep)
                pdist.outputTXT_defects_prot(f"Prot_{args.outputname}{ts.frame}", "shallow", "Up", DefectsUp_labels_group_Shallow, areaUp_clusters_Shallow)
                pdist.outputTXT_defects_prot(f"Prot_{args.outputname}{ts.frame}", "all", "Up", DefectsUp_labels_group_All, areaUp_clusters_All)

            # If there are proteins on the lower leaflet
            elif len(np.argwhere(MatrixLo_prot_All > 0)) > 0:
                labelsLo_coor_Deep = pdist.find_pd_border(MatrixLo_labels_Deep)
                labelsLo_coor_Shallow = pdist.find_pd_border(MatrixLo_labels_Shallow)
                labelsLo_coor_All = pdist.find_pd_border(MatrixLo_labels_All)

                edgeLo_coor_prot_Deep = pdist.find_prot_border(MatrixLo_prot_Deep)
                edgeLo_coor_prot_Shallow = pdist.find_prot_border(MatrixLo_prot_Shallow)
                edgeLo_coor_prot_All = pdist.find_prot_border(MatrixLo_prot_All)

                DefectsLo_labels_group_Deep = pdist.assign_dist_group(edgeLo_coor_prot_Deep, labelsLo_coor_Deep, 10)
                DefectsLo_labels_group_Shallow = pdist.assign_dist_group(edgeLo_coor_prot_Shallow, labelsLo_coor_Shallow, 10)
                DefectsLo_labels_group_All = pdist.assign_dist_group(edgeLo_coor_prot_All, labelsLo_coor_All, 10)

                pdist.outputTXT_defects_prot(f"Prot_{args.outputname}{ts.frame}", "deep", "Lo", DefectsLo_labels_group_Deep, areaLo_clusters_Deep)
                pdist.outputTXT_defects_prot(f"Prot_{args.outputname}{ts.frame}", "shallow", "Lo", DefectsLo_labels_group_Shallow, areaLo_clusters_Shallow)
                pdist.outputTXT_defects_prot(f"Prot_{args.outputname}{ts.frame}", "all", "Lo", DefectsLo_labels_group_All, areaLo_clusters_All)
