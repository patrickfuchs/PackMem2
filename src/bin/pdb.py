#-*- coding: utf-8 -*-
# Functions about PDB data and output files
# R. Gautier A. Bacle 2015
# M. Zygadlo 2025

import numpy as np


def write_a_pdb_line(nb_atm, atm_name, nb_res, coords, nb_defect):
    """
    Write a PDB line in a given file.

    --------------------
    INPUT
    file: file
        An open file to write into
    nb_atm: int
        The atom number
    atm_name: string
        The atom name
    nb_res: int
        The residue number
    coords: list
        Contains the coordinates x,y,z of the atom
    nb_defect: int
        The label of the defect
    """
    return f"{'ATOM  ':6s}{nb_atm:5d} {'  H1':4s} {atm_name:3s}  {nb_res:4d}    {float(coords[0]):8.3f}{float(coords[1]):8.3f}{float(coords[2]):8.3f}{1.0:6.2f}{nb_defect:6.2f}\n"

def outputPDB_Total_matrix(out_name, num_frame, arrayX, arrayY, z_extr, Matrix_final):
    """
    Create output file with Total Matrix in pdb format.

    --------------------
    INPUT
    out_name: string
        Name of the output PDB file
    num_frame: int
        The  frame number
    arrayX: numpy array
        Array from xmin-1 to xmax+1 by step of 1.0
    arrayY: numpy array
        Array from ymin-1 to ymax+1 by step of 1.0
    z_extr: float
        The maximum or  minimum z value
    Matrix_final: numpy array
        Contains the label of the defects or np.nan for the edges
    """
    with open(f'{out_name}.pdb',"w") as f_tot:
        f_tot.write(f"MODEL      {num_frame:3d}\n")
        nb=0
        for i, ind_matX in enumerate(arrayX):
            nb+=1
            for j, ind_matY in enumerate(arrayY):
                coordtmp = [ind_matX, ind_matY, z_extr]
                if np.isnan(Matrix_final[i][j]):
                    f_tot.write(write_a_pdb_line(nb, "EDG", nb, coordtmp, -1))
                else:
                    f_tot.write(write_a_pdb_line(nb, "MAT", nb, coordtmp, Matrix_final[i][j]))
        f_tot.write("ENDMDL\n")

def outputPDB_defects(out_name, num_frame, arrayX, arrayY, z_extr, 
                        Matrix_final, cluster_edge):
    """
    Create output file with defects only in pdb format.

    --------------------
    INPUT
    out_name: string
        Name of the output PDB file
    num_frame: int
        The  frame number
    arrayX: numpy array
        Array from xmin-1 to xmax+1 by step of 1.0
    arrayY: numpy array
        Array from ymin-1 to ymax+1 by step of 1.0
    z_extr: float
        The maximum or  minimum z value
    Matrix_final: numpy array
        Contains the label of the defects or np.nan for the edges
    edge_labels: list
        Contains the labels on the edge of the matrix
    """
    with open(f"{out_name}.pdb","w") as f:
        f.write(f"MODEL      {num_frame:3d}\n")
        nb=0
        dict_Def={}
        for i, ind_matX in enumerate(arrayX):
            for j, ind_matY in enumerate(arrayY):
                label_def = Matrix_final[i][j]
                if label_def !=0. and label_def not in cluster_edge:
                    if int(label_def) not in dict_Def:
                        dict_Def[int(label_def)]=[]
                    nb+=1
                    coordtmp=[ind_matX, ind_matY, z_extr]
                    dict_Def[int(label_def)].append(write_a_pdb_line(nb, 'DEF', int(label_def), coordtmp, label_def))
        keys_def= list(dict_Def.keys())
        keys_def.sort()
        nb=0
        for key in keys_def:
            for line in dict_Def[key]:
                nb+=1
                f.write(f"{line[:6]}{nb:5d}{line[11:]}")
        f.write("ENDMDL\n")

def outputTXT_defects(out_name, dict_def_area, dict_def_coor, 
                        total_size, total_edge, arrayX, arrayY):
    """
    Create output file for packing defects in TXT format.

    ## MatrixSize  9646  9801           # Membrane_matrix_size  Total_matrix_size
    ## Total   51   582 11.41 6.034     # Total_packing_defect_number, total_packing_defect_area, average_size, pourcent_of_membrane
    1    7     7.00    55.00            # Defect_label  defect_size  X_position  Y_position
    2    1     7.00    58.00

    --------------------
    INPUT
    out_name: string
        Name of the output result txt file
    dict_def_area: dictionary
        Contains the area of each label
    dict_def_coor: dictionary
        Contains the first appearance of the label in the matrix
    total_size: int
        The total area of the matrix
    total_edge: int
        The sum of the defect areas on the edge
    arrayX: numpy array
        Array from xmin-1 to xmax+1 by step of 1.0
    arrayY: numpy array
        Array from ymin-1 to ymax+1 by step of 1.0
    """
    with open(f"{out_name}.txt","w") as f_result:
        f_result.write(f"## MatrixSize {total_size-total_edge:5d} {total_size:5d} \n")
        if len(dict_def_coor) != 0:
            f_result.write(f"## Total {len(dict_def_area):4d} {sum(dict_def_area.values()):5d} {sum(dict_def_area.values())/float(len(dict_def_area)):4.2f} {(sum(dict_def_area.values())*100.)/(total_size-total_edge):5.3f}\n")
        else:
            #exception if no defect 
            f_result.write("## Total 0 0 0.0 0.0\n")
        i=0
        for key in dict_def_coor:
            f_result.write(f"{i+1:3d} {dict_def_area[key]:4d}   {arrayX[dict_def_coor[key][0]]:6.2f}   {arrayY[dict_def_coor[key][1]]:6.2f} \n")
            i+=1

