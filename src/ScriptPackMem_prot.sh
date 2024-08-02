#!/bin/python3

# Script to launch PackMem analysis on multiple PDB files.
# P. Fuchs, R. Gautier, 2016

# set PackMem absolute path
# change this path according to your system 
PackMemPATH=/home/mayazygadlo/These/tools/packmem_prot/New_PackMem/src

# set prefix of the pdb files
prefix=$1

# set the starting and ending frames
frame_start=$2
frame_stop=$3

###################
# No change below #
###################

# loop over the frames
# for example, if we have 400ns with one frame every 1000 ps,
# we loop from 0 to 400 frames (0-400 ns)
for pdbnum in $(seq ${frame_start} ${frame_stop})
do
    # print counter to screen
    echo "$(date): PackMem running on pdb/${prefix}${pdbnum}.pdb"
    # launch PackMem for the 3 types of packing defects
    ${PackMemPATH}/PackMem_prot.py -i ${PackMemPATH}/${prefix}/pdb/${prefix}${pdbnum}.pdb \
                              -r ${PackMemPATH}/vdw_radii_Charmm.txt \
                              -p ${PackMemPATH}/param_Charmm.txt \
                              -o ${prefix}${pdbnum} -d 1.0 -t deep
    ${PackMemPATH}/PackMem_prot.py -i ${PackMemPATH}/${prefix}/pdb/${prefix}${pdbnum}.pdb \
                              -r ${PackMemPATH}/vdw_radii_Charmm.txt \
                              -p ${PackMemPATH}/param_Charmm.txt \
                              -o ${prefix}${pdbnum} -d 1.0 -t shallow
    ${PackMemPATH}/PackMem_prot.py -i ${PackMemPATH}/${prefix}/pdb/${prefix}${pdbnum}.pdb \
                              -r ${PackMemPATH}/vdw_radii_Charmm.txt \
                              -p ${PackMemPATH}/param_Charmm.txt \
                              -o ${prefix}${pdbnum} -d 1.0 -t all
done
