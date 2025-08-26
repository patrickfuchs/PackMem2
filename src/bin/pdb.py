#-*- coding: utf-8 -*-
# Functions about PDB data and output files
# R. Gautier A. Bacle 2015
# M. Zygadlo 2025

import numpy as np


def write_a_pdb_line(file, nb_atm, atm_name, nb_res, coords, nb_defect):
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
    file.write(f"{'ATOM  ':6s}{nb_atm:5d} {'  H1':4s} {atm_name:3s}  {nb_res:4d}    {float(coords[0]):8.3f}{float(coords[1]):8.3f}{float(coords[2]):8.3f}{1.0:6.2f}{nb_defect:6.2f}\n")

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
    with open(f'{out_name}.pdb',"w") as f:
        f.write(f"MODEL      {num_frame:3d}\n")
        nb=0
        for j, ind_matX in enumerate(arrayX):
            nb+=1
            for k, ind_matY in enumerate(arrayY):
                coordtmp = [ind_matX, ind_matY, z_extr]
                if np.isnan(Matrix_final[j][k]):
                    write_a_pdb_line(f, nb, "EDG", nb, coordtmp, -1)
                else:
                    write_a_pdb_line(f, nb, "MAT", nb, coordtmp, Matrix_final[j][k])
        f.write("ENDMDL\n")


# create output file with defects only in pdb format
def outputPDB_defects(outputname, FlagPDtype, leaflet, num_frame, listX, listY, Rpos, 
                        Matrix_fin, cluster_edge):
    outputname = outputname + "_Defect" + leaflet+ "_"+ FlagPDtype +".pdb"
    with open(outputname,"w") as f:
        f.write("MODEL      %3d\n"%(num_frame))
        nb=0
        dicoDef={}
        for j, sliceX in enumerate(listX):
            for k, sliceY in enumerate(listY):
                num_def = Matrix_fin[j][k]
                if num_def !=0. and num_def not in cluster_edge:
                    if int(num_def) not in dicoDef:
                        dicoDef[int(num_def)]=[]
                    nb+=1
                    coordtmp=[sliceX,sliceY, Rpos]
                    dicoDef[int(num_def)].append(
                                "%6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%
                                ("ATOM  ", nb,"  H1", "DEF", int(num_def),
                                float(coordtmp[0]),float(coordtmp[1]),float(coordtmp[2]),
                                1.0,num_def))
        listdef = list(dicoDef.keys())
        listdef.sort()
        nb=0
        for l in listdef:
            for line in dicoDef[l]:
                nb+=1
                f.write("%s%5d%s"%(line[:6],nb,line[11:]))
        f.write("ENDMDL\n")


#create output file for packing defects in TXT format
## MatrixSize  9646  9801         # Membrane matrix size, Total matrix size
## Total   51   582 11.41 6.034   # number of packing defects, total area of packing defects, average size, pourcent of membrane (Membrane matrix size)
# 1    7     7.00    55.00        #defect size Xposition   Yposition
# 2    1     7.00    58.00
def outputTXT_defects(outputname, FlagPDtype, leaflet, dico_def_area, dico_def_coor, 
                        total_size, total_edge, listX, listY):
    outputname = outputname + "_" + leaflet+ "_"+FlagPDtype+"_result.txt"
    with open(outputname,"w") as f:
        f.write("## MatrixSize %5d %5d \n"%(total_size-total_edge,total_size))
        if len(dico_def_coor) != 0:
            f.write("## Total %4d %5d %4.2f %5.3f\n"%(len(dico_def_area),sum(dico_def_area.values()), sum(dico_def_area.values())/float(len(dico_def_area)), (sum(dico_def_area.values())*100.)/(total_size-total_edge)))
        else:
            #exception if no defect 
            f.write("## Total 0 0 0.0 0.0\n")
        i=0
        for key in dico_def_coor:
            # version with matrice position
            #f.write("%3d %4d   %6.2f   %6.2f \n"%(i+1, dico_def_area[key], dico_def_coor[key][0], dico_def_coor[key][1]))
            # version with x, y coordinates
            f.write("%3d %4d   %6.2f   %6.2f \n"%(i+1, dico_def_area[key], listX[dico_def_coor[key][0]], listY[dico_def_coor[key][1]]))
            i+=1

