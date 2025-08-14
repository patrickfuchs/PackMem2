#!/usr/bin/python3
#-*- coding: utf-8 -*-
# Pg compute packing defects on membranes simulations
# R. Gautier A. Bacle april 2016
# R. Gerard may 2024
# M. Zygadlo august 2024

import sys
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

##########################################################################################
##### main
if __name__ == '__main__':
    """
    Python
    Script to compute Packing defect in flat bilayers
    Lipids parameters adapted to either Berger lipid (with corrections), CHARMM36 (with Klauda corrections) or Martini FF. Be careful to the atoms name if you use other lipids
    fileRadius.txt example:
    DOPC  C02 1.875 a (aliphatic)
    DOPC  O8  1.48 n (non aliphatic)
    DOPC  C25 1.98 n (non aliphatic)

    Output files:
        outputname_Up/Lo_Shallow/Deep/All_result.txt
        ## MatrixSize  9646  9801         # Membrane matrix size, Total matrix size
        ## Total   51   582 11.41 6.034   # number of packing defects, total area of packing defects, average size, pourcent of membrane (Membrane matrix size)
        1    1    77.00     2.00          # for each packing defect num size x_position (first pixel) y_mean (first pixel)
        2    1    21.00     3.00 
        3    8    59.12     7.00 

        if -pdb option:
            outputname_TotalUp/Lo_Shallow/Deep/All.pdb 
            Matrix x,y with, for each cell, the value of "packing defect" (in B_factor column, the last column) 
            if 0 = Deep packing defects The value increases with the number of atoms in the cell.
            if > 0 and < 1 = Shallow defect
            if -1 = Edge
            ATOM      3   H1 EDG     3       4.000  54.000  56.940  1.00 -1.00
            ATOM      3   H1 MAT     3       4.000  54.000  56.940  1.00  0.00
            ATOM      3   H1 MAT     3       4.000  54.000  56.940  1.00  0.00

            outputname+"_DefectUp/Lo_Shallow/Deep/All.pdb
            Packing defects in pdb format:
            the residue number corresponds to the different packing defects.
            ATOM      5   H1 DEF     2       6.000  22.000  56.940  1.00  2.00
            ATOM      5   H1 DEF     3       6.000  85.000  56.940  1.00  3.00
            ATOM      5   H1 DEF     3       6.000  86.000  56.940  1.00  3.00
            ATOM      6   H1 DEF     3       7.000  86.000  56.940  1.00  3.00

        if -prot option:
            Prot_outputnameFrame_Up/Lo_Deep/Shallow/All_prot.txt
    """
    ####### PARAMETERS and INPUT #####
    try:
        args = p.get_args()

    except (Exception, FileNotFoundError) as error:
        # Don't print Traceback 
        sys.exit(error)
    

    ######### READ PARAM FILES #########
    # RESNAME_GLYC = {'DOP': 'C2', 'DOE': 'C2', 'DPP': 'C2', etc}
    RESNAME_GLYC = p.set_params(args.paramFile)
    # {'ALA N': 1.85, 'ALA HN': 0.22, 'ALA HT1': 0.22, etc}
    # for the amino acids then the lipids
    radius = bfrg.read_radius(args.filesrad)
    # {'ALA N': 'n', 'ALA HN': 'n', 'ALA HT1': 'n', etc }
    # for the amino acids then the lipids
    aliphatic = bfrg.read_aliphatic(args.filesrad)


    ######### LOAD UNIVERSE  #########
    u = mda.Universe(args.topo, args.traj)
    # If multiple lipids
    lipid_names = args.lipid.replace('_', ' ')
    # select the lipids in system
    lipids = u.select_atoms(f"resname {lipid_names}")
    

    ######### Extract UPPER/LOWER leaflet ##########
    # Get the glycerol atom name(s)
    glyc_mb = l.get_glyc_lipids(args.lipid.split('_'), RESNAME_GLYC)
    # If no index file seperating the leaflets
    if args.indexFile == None:
        # Create lists of the residue number  for upper and lower leaflets
        L = LeafletFinder(lipids, f'name {glyc_mb}')
        upper_leaflet = L.groups(0).resids
        lower_leaflet = L.groups(1).resids
    else:
        lower_leaflet, upper_leaflet = p.read_ndx(args.indexFile)


    ############################## Main loop ##################################
    for ts in u.trajectory[args.start:args.end+1]:
        # select all atoms in systems
        system = u.select_atoms(f"resname {lipid_names} or protein")

        # Get all the residues in the membrane
        res_ids = system.resids
        md_uniq_ids = list(set(res_ids))

        # If there are a lot of lipids (> 9999)
        if len(md_uniq_ids) > 9999:
            print("The number of lipids in your membrane is very high (>9999)!")
            print("Please check that you only give to PackMem your membrane (+ protein)")
            sys.exit()

        # Get all x, y and z
        x_atoms = system.positions[:,0].round(2)
        y_atoms = system.positions[:,1].round(2)
        z_atoms = system.positions[:,2].round(2)
        coords = np.stack((x_atoms, y_atoms, z_atoms), axis=1)
        # Get membrane dimension
        xmin, xmax, xmean = l.min_max(x_atoms)
        ymin, ymax, ymean = l.min_max(y_atoms)
        zmin, zmax, zmean = l.min_max(z_atoms)
        
        upper_arrayZ = l.create_arrayZ(system.residues, upper_leaflet, glyc_mb, args.dist_suppl_Z, zmax)
        lower_arrayZ = l.create_arrayZ(system.residues, lower_leaflet, glyc_mb, args.dist_suppl_Z, zmin, up=False)

        # Build an array from xmin-1 to xmax+1 every 1.0
        arrayX = l.create_array(int(xmin-1), int(xmax+2), m.SIZE)
        # Build an array from ymin-1 to ymax+1 every 1.0
        arrayY = l.create_array(int(ymin-1), int(ymax+2), m.SIZE)


        ####################  Compute Matrix    #################
        # Initialize 2 matrixes for Upper and Lower leaflet of length arrayX,arrayY
        MatrixUp_Deep = m.initialize_matrix2D(len(arrayX), len(arrayY), np.nan)
        MatrixLo_Deep = m.initialize_matrix2D(len(arrayX), len(arrayY), np.nan)

        MatrixUp_Shallow = m.initialize_matrix2D(len(arrayX), len(arrayY), np.nan)
        MatrixLo_Shallow = m.initialize_matrix2D(len(arrayX), len(arrayY), np.nan)

        MatrixUp_All = m.initialize_matrix2D(len(arrayX), len(arrayY), np.nan)
        MatrixLo_All = m.initialize_matrix2D(len(arrayX), len(arrayY), np.nan)
        
        # Search v cells around coord
        v = 5.0
        # For each atoms of lipids
        for i, (res_id, atom_name, res_name) in enumerate(zip(res_ids, system.names, system.resnames)):
            radius_atm = m.get_radius(radius, res_name, atom_name)
            aliph_atom = m.get_aliphatic(aliphatic, res_name, atom_name)
            coordtmp = coords[i]
            #### Upper leaflet ####
            if res_id in upper_leaflet :
                # dZ = z_C2_res - z_atom
                dZ = round(m.diff_Z(upper_arrayZ[res_id], coordtmp[2]), 2)
                # If dZ < 5.0
                # To limit the search around an atom to 5 cells
                if dZ < v:
                    # Fill the matrix with value 0 < a < 1 for aliphatic
                    # Or with > 1 if polar OR deep
                    # Defects = 0
                    MatrixUp_Deep = m.fill_matrix(MatrixUp_Deep, coordtmp, arrayX, arrayY,
                                            upper_arrayZ[res_id], radius_atm,
                                            "deep", aliph_atom)
                    MatrixUp_Shallow = m.fill_matrix(MatrixUp_Shallow, coordtmp, arrayX, arrayY,
                                            upper_arrayZ[res_id], radius_atm,
                                            "shallow", aliph_atom)
                    MatrixUp_All = m.fill_matrix(MatrixUp_All, coordtmp, arrayX, arrayY,
                                            upper_arrayZ[res_id], radius_atm,
                                            "all", aliph_atom)
            #### Lower leaflet ####
            if res_id in lower_leaflet :
                # dZ = z_C2_res - z_atom
                dZ = round(m.diff_Z(lower_arrayZ[res_id], coordtmp[2]), 2)
                # If dfZ > -5.0
                # To limit the search around an atom to 5 cells
                if dZ > -v:
                    # Fill the matrix with value 0 < a < 1 for aliphatic
                    # Or with > 1 if polar OR deep
                    # Defects = 0
                    MatrixLo_Deep = m.fill_matrix(MatrixLo_Deep, coordtmp, arrayX, arrayY,
                                            lower_arrayZ[res_id], radius_atm,
                                            "deep", aliph_atom)
                    MatrixLo_Shallow = m.fill_matrix(MatrixLo_Shallow, coordtmp, arrayX, arrayY,
                                            lower_arrayZ[res_id], radius_atm,
                                            "shallow", aliph_atom)
                    MatrixLo_All = m.fill_matrix(MatrixLo_All, coordtmp, arrayX, arrayY,
                                            lower_arrayZ[res_id], radius_atm,
                                            "all", aliph_atom)


        ####################  Binarise the matrix    #################
        # Preliminary process for shallow defect and problem of the edges 
        # To eliminate shallow defects on edges first: binarize on all defects and storage edges coord
        # Initalise matrices to 0.0
        MatrixUp_Shallowbin = m.initialize_matrix2D(len(arrayX), len(arrayY), 0.)
        MatrixLo_Shallowbin = m.initialize_matrix2D(len(arrayX), len(arrayY), 0.)
        # Binarise these matrices, with 0 for aliphatic atoms + packing defects and 1 otherwise
        MatrixUp_Shallowbin = m.binarize_matrix_without0(MatrixUp_Shallow, MatrixUp_Shallowbin, -0.01, 0.99)
        MatrixLo_Shallowbin = m.binarize_matrix_without0(MatrixLo_Shallow, MatrixLo_Shallowbin, -0.01, 0.99)
    
        # Get temporary packing defects
        # Connect the packing defects + label them + count the area
        MatrixUp_labels_Shallow, set_labelsUp_Shallow, area_labelsUp_Shallow, firstCoorUp_labels_Shallow = \
            cc.get_connected_components(MatrixUp_Shallowbin)
        MatrixLo_labels_Shallow, set_labelsLo_Shallow, area_labelsLo_Shallow, firstCoorLo_labels_Shallow = \
            cc.get_connected_components(MatrixLo_Shallowbin)
        # Get the labels on the edge
        edge_labelsUp_Shallow=cc.get_clusters_on_the_edge(MatrixUp_labels_Shallow)
        edge_labelsLo_Shallow=cc.get_clusters_on_the_edge(MatrixLo_labels_Shallow)


        # Initalise matrices to 0.0
        MatrixUp_Deepbin = m.initialize_matrix2D(len(arrayX), len(arrayY), 0.)
        MatrixLo_Deepbin = m.initialize_matrix2D(len(arrayX), len(arrayY), 0.)
        MatrixUp_Shallowbin = m.initialize_matrix2D(len(arrayX), len(arrayY), 0.)
        MatrixLo_Shallowbin = m.initialize_matrix2D(len(arrayX), len(arrayY), 0.)
        MatrixUp_Allbin = m.initialize_matrix2D(len(arrayX), len(arrayY), 0.)
        MatrixLo_Allbin = m.initialize_matrix2D(len(arrayX), len(arrayY), 0.)

        #### Deep ####
        # Binarise these matrices, with 0 if deep defect, 1 otherwise
        MatrixUp_Deepbin = m.binarize_matrix_without0(MatrixUp_Deep, MatrixUp_Deepbin, -0.01, 0.001)
        MatrixLo_Deepbin = m.binarize_matrix_without0(MatrixLo_Deep, MatrixLo_Deepbin, -0.01, 0.001)
        # Packing defects determination
        # Connect the packing defects + label them + count the area
        MatrixUp_labels_Deep, set_labelsUp_Deep, area_labelsUp_Deep, firstCoorUp_labels_Deep = \
            cc.get_connected_components(MatrixUp_Deepbin)
        MatrixLo_labels_Deep, set_labelsLo_Deep, area_labelsLo_Deep, firstCoorLo_labels_Deep = \
            cc.get_connected_components(MatrixLo_Deepbin)
        # Get cluster on the edge
        edge_labelsUp_Deep = cc.get_clusters_on_the_edge(MatrixUp_labels_Deep)
        edge_labelsLo_Deep = cc.get_clusters_on_the_edge(MatrixLo_labels_Deep)
        # Count area of the edge
        area_edgeUp_Deep = m.count_edge_area(area_labelsUp_Deep, edge_labelsUp_Deep)
        area_edgeLo_Deep = m.count_edge_area(area_labelsLo_Deep, edge_labelsLo_Deep)
        # Clean dico defects (without edge)
        area_labelsUp_Deep = d.del_key_dico(area_labelsUp_Deep, edge_labelsUp_Deep)
        firstCoorUp_labels_Deep = d.del_key_dico(firstCoorUp_labels_Deep, edge_labelsUp_Deep)
        area_labelsLo_Deep = d.del_key_dico(area_labelsLo_Deep, edge_labelsLo_Deep)
        firstCoorLo_labels_Deep = d.del_key_dico(firstCoorLo_labels_Deep, edge_labelsLo_Deep)
        

        #### Shallow ####
        # Binarise these matrices, with 0 for aliphatic atoms and 1 otherwise
        MatrixUp_Shallowbin = m.binarize_matrix_without0(MatrixUp_Shallow, MatrixUp_Shallowbin, 0, 0.99)
        MatrixLo_Shallowbin = m.binarize_matrix_without0(MatrixLo_Shallow, MatrixLo_Shallowbin, 0, 0.99)
        # Modify the binary matrix to take account edges (determined by all packing defects)
        # The clusters that had their labels on the edge are put to 0.0
        MatrixUp_Shallowbin = m.modify_matrix(MatrixUp_labels_Shallow, MatrixUp_Shallowbin, edge_labelsUp_Shallow)
        MatrixLo_Shallowbin = m.modify_matrix(MatrixLo_labels_Shallow, MatrixLo_Shallowbin, edge_labelsLo_Shallow)
        # Packing defects determination
        # Connect the packing defects + label them + count the area
        MatrixUp_labels_Shallow, set_labelsUp_Shallow, area_labelsUp_Shallow, firstCoorUp_labels_Shallow = \
            cc.get_connected_components(MatrixUp_Shallowbin)
        MatrixLo_labels_Shallow, set_labelsLo_Shallow, area_labelsLo_Shallow, firstCoorLo_labels_Shallow = \
            cc.get_connected_components(MatrixLo_Shallowbin)
        # Get cluster on the edge
        edge_labelsUp_Shallow = cc.get_clusters_on_the_edge(MatrixUp_labels_Shallow)
        edge_labelsLo_Shallow = cc.get_clusters_on_the_edge(MatrixLo_labels_Shallow)
        # Count area of the edge
        area_edgeUp_Shallow = m.count_edge_area(area_labelsUp_Shallow, edge_labelsUp_Shallow)
        area_edgeLo_Shallow = m.count_edge_area(area_labelsLo_Shallow, edge_labelsLo_Shallow)
        # Clean dico defects (without edge)
        area_labelsUp_Shallow = d.del_key_dico(area_labelsUp_Shallow, edge_labelsUp_Shallow)
        firstCoorUp_labels_Shallow = d.del_key_dico(firstCoorUp_labels_Shallow, edge_labelsUp_Shallow)
        area_labelsLo_Shallow = d.del_key_dico(area_labelsLo_Shallow, edge_labelsLo_Shallow)
        firstCoorLo_labels_Shallow = d.del_key_dico(firstCoorLo_labels_Shallow, edge_labelsLo_Shallow)
        # Eliminate nan inside (deep not shallow defect)
        MatrixUp_labels_Shallow, area_edgeUp_Shallow, labelsPb_Up_Shallow = \
            m.clean_NA_inside(MatrixUp_labels_Shallow, edge_labelsUp_Shallow,
                              MatrixUp_Shallow, area_edgeUp_Shallow)
        set_labelsUp_Shallow, area_labelsUp_Shallow, firstCoorUp_labels_Shallow = \
            cc.delete_NApoints_inside(labelsPb_Up_Shallow, MatrixUp_labels_Shallow,
                                    set_labelsUp_Shallow, area_labelsUp_Shallow)

        MatrixLo_labels_Shallow, area_edgeLo_Shallow, labelsPb_Lo_Shallow = \
            m.clean_NA_inside(MatrixLo_labels_Shallow, edge_labelsLo_Shallow,
                              MatrixLo_Shallow, area_edgeLo_Shallow)
        set_labelsLo_Shallow, area_labelsLo_Shallow, firstCoorLo_labels_Shallow = \
            cc.delete_NApoints_inside(labelsPb_Lo_Shallow, MatrixLo_labels_Shallow,
                                    set_labelsLo_Shallow, area_labelsLo_Shallow)
        # Reclean dico defects (without edge)
        area_labelsUp_Shallow = d.del_key_dico(area_labelsUp_Shallow, edge_labelsUp_Shallow)
        firstCoorUp_labels_Shallow = d.del_key_dico(firstCoorUp_labels_Shallow, edge_labelsUp_Shallow)
        area_labelsLo_Shallow = d.del_key_dico(area_labelsLo_Shallow, edge_labelsLo_Shallow)
        firstCoorLo_labels_Shallow = d.del_key_dico(firstCoorLo_labels_Shallow, edge_labelsLo_Shallow)


        #### All ####
        # Binarise these matrices, with 0 for aliphatic atoms + packing defects and 1 otherwise
        MatrixUp_Allbin = m.binarize_matrix_without0(MatrixUp_All, MatrixUp_Allbin, -0.01, 0.99)
        MatrixLo_Allbin = m.binarize_matrix_without0(MatrixLo_All, MatrixLo_Allbin, -0.01, 0.99)
        # Packing defects determination
        # Connect the packing defects + label them + count the area
        MatrixUp_labels_All, set_labelsUp_All, area_labelsUp_All, firstCoorUp_labels_All = \
            cc.get_connected_components(MatrixUp_Allbin)
        MatrixLo_labels_All, set_labelsLo_All, area_labelsLo_All, firstCoorLo_labels_All = \
            cc.get_connected_components(MatrixLo_Allbin)
        # Get cluster on the edge
        edge_labelsUp_All = cc.get_clusters_on_the_edge(MatrixUp_labels_All)
        edge_labelsLo_All = cc.get_clusters_on_the_edge(MatrixLo_labels_All)
        # Count area of the edge
        area_edgeUp_All = m.count_edge_area(area_labelsUp_All, edge_labelsUp_All)
        area_edgeLo_All = m.count_edge_area(area_labelsLo_All, edge_labelsLo_All)
        # Clean dico defects (without edge)
        area_labelsUp_All = d.del_key_dico(area_labelsUp_All, edge_labelsUp_All)
        firstCoorUp_labels_All = d.del_key_dico(firstCoorUp_labels_All, edge_labelsUp_All)
        area_labelsLo_All = d.del_key_dico(area_labelsLo_All, edge_labelsLo_All)
        firstCoorLo_labels_All = d.del_key_dico(firstCoorLo_labels_All, edge_labelsLo_All)


        ####################  Output text file  #################
        # Compute the total area of the matrix
        total_area = len(arrayX) * len(arrayY)
        pdb.outputTXT_defects(f"{args.outputname}{ts.frame}", "deep", "Up", area_labelsUp_Deep, 
                            firstCoorUp_labels_Deep, total_area, area_edgeUp_Deep, arrayX, arrayY)
        pdb.outputTXT_defects(f"{args.outputname}{ts.frame}", "deep", "Lo", area_labelsLo_Deep, 
                            firstCoorLo_labels_Deep, total_area, area_edgeLo_Deep, arrayX, arrayY)
        
        pdb.outputTXT_defects(f"{args.outputname}{ts.frame}", "shallow", "Up", area_labelsUp_Shallow, 
                            firstCoorUp_labels_Shallow, total_area, area_edgeUp_Shallow, arrayX, arrayY)
        pdb.outputTXT_defects(f"{args.outputname}{ts.frame}", "shallow", "Lo", area_labelsLo_Shallow, 
                            firstCoorLo_labels_Shallow, total_area, area_edgeLo_Shallow, arrayX, arrayY)
        
        pdb.outputTXT_defects(f"{args.outputname}{ts.frame}", "all", "Up", area_labelsUp_All, 
                            firstCoorUp_labels_All, total_area, area_edgeUp_All, arrayX, arrayY)
        pdb.outputTXT_defects(f"{args.outputname}{ts.frame}", "all", "Lo", area_labelsLo_All, 
                            firstCoorLo_labels_All, total_area, area_edgeLo_All, arrayX, arrayY)


        ####################  Output PDB files  #################
        # final matrix values PD (X,Y) with Z cooresponding to the max(Upper/lowerZlevel)
        if args.pdbout :
            # Get the max/min of the z_coord+1
            valZmax=float(d.max_value_dico(upper_arrayZ))
            valZmin=float(d.min_value_dico(lower_arrayZ))

            # Write the leaflets into a PDB
            # To ignore the warnings when writing a PDB
            warnings.filterwarnings("ignore", category=UserWarning)
            # Write the leaflets into a PDB
            u.select_atoms(f"resid {upper_leaflet[0]}:{upper_leaflet[-1]}").write(f"{args.outputname}{ts.frame}_Upper_leaflet.pdb")
            u.select_atoms(f"resid {lower_leaflet[0]}:{lower_leaflet[-1]}").write(f"{args.outputname}{ts.frame}_Lower_leaflet.pdb")
            
            # Write the Matrix (each cell) into a PDB
            pdb.outputPDB_Total_matrix(f"{args.outputname}{ts.frame}", "deep", "Up", ts.frame,
                                    arrayX, arrayY, valZmax, MatrixUp_Deep)
            pdb.outputPDB_Total_matrix(f"{args.outputname}{ts.frame}", "deep", "Lo", ts.frame,
                                    arrayX, arrayY, valZmin, MatrixLo_Deep)
            pdb.outputPDB_Total_matrix(f"{args.outputname}{ts.frame}", "shallow", "Up", ts.frame,
                                    arrayX, arrayY, valZmax, MatrixUp_Shallow)
            pdb.outputPDB_Total_matrix(f"{args.outputname}{ts.frame}", "shallow", "Lo", ts.frame,
                                    arrayX, arrayY, valZmin, MatrixLo_Shallow)
            pdb.outputPDB_Total_matrix(f"{args.outputname}{ts.frame}", "all", "Up", ts.frame,
                                    arrayX, arrayY, valZmax, MatrixUp_All)
            pdb.outputPDB_Total_matrix(f"{args.outputname}{ts.frame}", "all", "Lo", ts.frame,
                                    arrayX, arrayY, valZmin, MatrixLo_All)
            
            # Write the defects into a PDB
            pdb.outputPDB_defects(f"{args.outputname}{ts.frame}", "deep", "Up", ts.frame,
                                arrayX, arrayY, valZmax, MatrixUp_labels_Deep, edge_labelsUp_Deep)
            pdb.outputPDB_defects(f"{args.outputname}{ts.frame}", "deep", "Lo", ts.frame,
                                arrayX, arrayY, valZmin, MatrixLo_labels_Deep, edge_labelsLo_Deep)
            pdb.outputPDB_defects(f"{args.outputname}{ts.frame}", "shallow", "Up", ts.frame,
                                arrayX, arrayY, valZmax, MatrixUp_labels_Shallow, edge_labelsUp_Shallow)
            pdb.outputPDB_defects(f"{args.outputname}{ts.frame}", "shallow", "Lo", ts.frame,
                                arrayX, arrayY, valZmin, MatrixLo_labels_Shallow, edge_labelsLo_Shallow)
            pdb.outputPDB_defects(f"{args.outputname}{ts.frame}", "all", "Up", ts.frame,
                                arrayX, arrayY, valZmax, MatrixUp_labels_All, edge_labelsUp_All)
            pdb.outputPDB_defects(f"{args.outputname}{ts.frame}", "all", "Lo", ts.frame,
                                arrayX, arrayY, valZmin, MatrixLo_labels_All, edge_labelsLo_All)
        

        ####################  Distance from the protein  #################
        if args.protein :
            # Select only the protein
            protein = u.select_atoms("protein")
            pos_prot = protein.positions

            # Label 1 corresponds to the void around the simulation
            # So we replace it by 0
            MatrixUp_labels_Deep = np.where(MatrixUp_labels_Deep == 1, 0, MatrixUp_labels_Deep)
            MatrixLo_labels_Deep = np.where(MatrixLo_labels_Deep == 1, 0, MatrixLo_labels_Deep)
            MatrixUp_labels_Shallow = np.where(MatrixUp_labels_Shallow == 1, 0, MatrixUp_labels_Shallow)
            MatrixLo_labels_Shallow = np.where(MatrixLo_labels_Shallow == 1, 0, MatrixLo_labels_Shallow)
            MatrixUp_labels_All = np.where(MatrixUp_labels_All == 1, 0, MatrixUp_labels_All)
            MatrixLo_labels_All = np.where(MatrixLo_labels_All == 1, 0, MatrixLo_labels_All)

            # Create empty arrays
            MatrixUp_prot_Deep = m.initialize_matrix2D(len(arrayX), len(arrayY), 0)
            MatrixLo_prot_Deep = m.initialize_matrix2D(len(arrayX), len(arrayY), 0)
            MatrixUp_prot_Shallow = m.initialize_matrix2D(len(arrayX), len(arrayY), 0)
            MatrixLo_prot_Shallow = m.initialize_matrix2D(len(arrayX), len(arrayY), 0)
            MatrixUp_prot_All = m.initialize_matrix2D(len(arrayX), len(arrayY), 0)
            MatrixLo_prot_All = m.initialize_matrix2D(len(arrayX), len(arrayY), 0)

            # Find where the protein is in the matrix
            MatrixUp_prot_Deep = pdist.find_protein(MatrixUp_prot_Deep, 'up', protein, arrayX, arrayY, zmean)
            MatrixUp_prot_Shallow = pdist.find_protein(MatrixUp_prot_Shallow, 'up', protein, arrayX, arrayY, zmean)
            MatrixUp_prot_All = pdist.find_protein(MatrixUp_prot_All, 'up', protein, arrayX, arrayY, zmean)
            MatrixLo_prot_Deep = pdist.find_protein(MatrixLo_prot_Deep, 'lo', protein, arrayX, arrayY, zmean)
            MatrixLo_prot_Shallow = pdist.find_protein(MatrixLo_prot_Shallow, 'lo', protein, arrayX, arrayY, zmean)
            MatrixLo_prot_All = pdist.find_protein(MatrixLo_prot_All, 'lo', protein, arrayX, arrayY, zmean)
            

            # If there are proteins on the upper leaflet
            if len(np.argwhere(MatrixUp_prot_All > 0)) > 0:
                # Classification of Packing Defects by distance group 
                # Get the coordinates of the matrix where the edges of the packing defects are located.
                # dictionnary label coords {lab1 = [(x1, y1), (x2, y2), ...],
                #                           lab2 = [(x1, y1), (x2, y2), ...], 
                #                           ...                              }
                coor_label_bordersUp_Deep = pdist.find_pd_border(MatrixUp_labels_Deep)
                coor_label_bordersUp_Shallow = pdist.find_pd_border(MatrixUp_labels_Shallow)
                coor_label_bordersUp_All = pdist.find_pd_border(MatrixUp_labels_All)

                # Get the coordinates of the matrix where the edges of the protein are located.
                # list of tuples [(x1, y1), (x2, y2), ...]
                coor_prot_edgeUp_Deep = pdist.find_prot_border(MatrixUp_prot_Deep)
                coor_prot_edgeUp_Shallow = pdist.find_prot_border(MatrixUp_prot_Shallow)
                coor_prot_edgeUp_All = pdist.find_prot_border(MatrixUp_prot_All)

                # Assign distance group for each packing defect, "far" or "close". Default threshold = 10 A.
                # dict {lab1 : 'group', lab2 : 'group', ... }
                DefectsUp_labels_group_Deep = pdist.assign_dist_group(coor_prot_edgeUp_Deep, coor_label_bordersUp_Deep, 10)
                DefectsUp_labels_group_Shallow = pdist.assign_dist_group(coor_prot_edgeUp_Shallow, coor_label_bordersUp_Shallow, 10)
                DefectsUp_labels_group_All = pdist.assign_dist_group(coor_prot_edgeUp_All, coor_label_bordersUp_All, 10)

                # Write the result in a text file
                # format : label,dist_group,area
                pdist.outputTXT_defects_prot(f"Prot_{args.outputname}{ts.frame}", "deep", "Up", DefectsUp_labels_group_Deep, area_labelsUp_Deep)
                pdist.outputTXT_defects_prot(f"Prot_{args.outputname}{ts.frame}", "shallow", "Up", DefectsUp_labels_group_Shallow, area_labelsUp_Shallow)
                pdist.outputTXT_defects_prot(f"Prot_{args.outputname}{ts.frame}", "all", "Up", DefectsUp_labels_group_All, area_labelsUp_All)

            # If there are proteins on the lower leaflet
            elif len(np.argwhere(MatrixLo_prot_All > 0)) > 0:
                coor_label_bordersLo_Deep = pdist.find_pd_border(MatrixLo_labels_Deep)
                coor_label_bordersLo_Shallow = pdist.find_pd_border(MatrixLo_labels_Shallow)
                coor_label_bordersLo_All = pdist.find_pd_border(MatrixLo_labels_All)

                coor_prot_edgeLo_Deep = pdist.find_prot_border(MatrixLo_prot_Deep)
                coor_prot_edgeLo_Shallow = pdist.find_prot_border(MatrixLo_prot_Shallow)
                coor_prot_edgeLo_All = pdist.find_prot_border(MatrixLo_prot_All)

                DefectsLo_labels_group_Deep = pdist.assign_dist_group(coor_prot_edgeLo_Deep, coor_label_bordersLo_Deep, 10)
                DefectsLo_labels_group_Shallow = pdist.assign_dist_group(coor_prot_edgeLo_Shallow, coor_label_bordersLo_Shallow, 10)
                DefectsLo_labels_group_All = pdist.assign_dist_group(coor_prot_edgeLo_All, coor_label_bordersLo_All, 10)

                pdist.outputTXT_defects_prot(f"Prot_{args.outputname}{ts.frame}", "deep", "Lo", DefectsLo_labels_group_Deep, area_labelsLo_Deep)
                pdist.outputTXT_defects_prot(f"Prot_{args.outputname}{ts.frame}", "shallow", "Lo", DefectsLo_labels_group_Shallow, area_labelsLo_Shallow)
                pdist.outputTXT_defects_prot(f"Prot_{args.outputname}{ts.frame}", "all", "Lo", DefectsLo_labels_group_All, area_labelsLo_All)
