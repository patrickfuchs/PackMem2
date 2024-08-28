#!/bin/bash

# Script to launch PackMem analysis on multiple PDB files.
# P. Fuchs, R. Gautier, 2016

# set PackMem absolute path
# change this path according to your system 
PackMemPATH=/home/maya/Documents/tools/New_PackMem/

# Get the files needed
traj=$1
topo=$2

# set the starting and ending frames
frame_start=$3
frame_stop=$4

###################
# No change below #
###################
# print counter to screen
#echo "$(date): PackMem running on pdb/${prefix}${pdbnum}.pdb"
# launch PackMem for the 3 types of packing defects
python ${PackMemPATH}/src/PackMem_prot.py -f ${traj} -s ${topo} \
                            -r ${PackMemPATH}/data/vdw_radii_Charmm.txt \
                            -p ${PackMemPATH}/data/param_Charmm.txt \
                            -o packmemout -d 1.0 -t deep -b ${frame_start} -e ${frame_stop}
python ${PackMemPATH}/src/PackMem_prot.py -f ${traj} -s ${topo} \
                            -r ${PackMemPATH}/data/vdw_radii_Charmm.txt \
                            -p ${PackMemPATH}/data/param_Charmm.txt \
                            -o packmemout -d 1.0 -t shallow -b ${frame_start} -e ${frame_stop}
python ${PackMemPATH}/src/PackMem_prot.py -f ${traj} -s ${topo} \
                            -r ${PackMemPATH}/data/vdw_radii_Charmm.txt \
                            -p ${PackMemPATH}/data/param_Charmm.txt \
                            -o packmemout -d 1.0 -t all -b ${frame_start} -e ${frame_stop}
