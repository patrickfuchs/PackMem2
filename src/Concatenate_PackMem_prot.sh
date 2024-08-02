#!/bin/bash

# Script to concatenate the PackMem analysis on multiple PDB files.

# set prefix of the pdb files
prefix=$1

# set the first frame
first_frame=$2

# set the last frame
last_frame=$3


# The following Total_*.txt files will contain the statistics for each type of
# packing defects accumulated over all frames (one file per membrane leaflet)
# Delete those files if they already exist before launching the main loop
# (the -f flag avoids an error if the file doesn't exist)
#rm -f Total_Up_${prefix}_deep.txt
#rm -f Total_Lo_${prefix}_deep.txt
#rm -f Total_Up_${prefix}_shallow.txt
#rm -f Total_Lo_${prefix}_shallow.txt
#rm -f Total_Up_${prefix}_all.txt
#rm -f Total_Lo_${prefix}_all.txt



for pdbnum in $(seq ${first_frame} ${last_frame})
do
    # Accumulate packing defects of the current frame in Total_*.txt files per leaflet
     cat  ${prefix}${pdbnum}_Up_All_prot.txt >> Total_Up_${prefix}_all_prot.txt
     cat  ${prefix}${pdbnum}_Up_Deep_prot.txt >> Total_Up_${prefix}_deep_prot.txt
     cat  ${prefix}${pdbnum}_Up_Shallow_prot.txt >> Total_Up_${prefix}_shallow_prot.txt

    # we no longer need the defects of the current frame
     rm -f ${prefix}${pdbnum}_Up_All_prot.txt
     rm -f ${prefix}${pdbnum}_Up_Deep_prot.txt
     rm -f ${prefix}${pdbnum}_Up_Shallow_prot.txt
done