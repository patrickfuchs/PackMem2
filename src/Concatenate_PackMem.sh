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
     cat  ${prefix}${pdbnum}_Up_Deep_result.txt >> Total_Up_${prefix}_deep.txt
     cat  ${prefix}${pdbnum}_Lo_Deep_result.txt >> Total_Lo_${prefix}_deep.txt
     cat  ${prefix}${pdbnum}_Up_Shallow_result.txt >> Total_Up_${prefix}_shallow.txt
     cat  ${prefix}${pdbnum}_Lo_Shallow_result.txt >> Total_Lo_${prefix}_shallow.txt
     cat  ${prefix}${pdbnum}_Up_All_result.txt >> Total_Up_${prefix}_all.txt
     cat  ${prefix}${pdbnum}_Lo_All_result.txt >> Total_Lo_${prefix}_all.txt
    # we no longer need the defects of the current frame
     rm -f ${prefix}${pdbnum}_Up_Deep_result.txt
     rm -f ${prefix}${pdbnum}_Lo_Deep_result.txt
     rm -f ${prefix}${pdbnum}_Up_Shallow_result.txt
     rm -f ${prefix}${pdbnum}_Lo_Shallow_result.txt
     rm -f ${prefix}${pdbnum}_Up_All_result.txt
     rm -f ${prefix}${pdbnum}_Lo_All_result.txt
done

# Concatenate the packing defects statistics of the lower and upper leaflets
cat Total_Up_${prefix}_deep.txt Total_Lo_${prefix}_deep.txt > Total_${prefix}_deep.txt
cat Total_Up_${prefix}_shallow.txt Total_Lo_${prefix}_shallow.txt > Total_${prefix}_shallow.txt
cat Total_Up_${prefix}_all.txt Total_Lo_${prefix}_all.txt > Total_${prefix}_all.txt

# Clean txt file for python script
egrep -v "#" Total_${prefix}_deep.txt > Total_${prefix}_deep_clean.txt
egrep -v "#" Total_${prefix}_shallow.txt > Total_${prefix}_shallow_clean.txt
egrep -v "#" Total_${prefix}_all.txt > Total_${prefix}_all_clean.txt

# Clean leaflet txt file for python script
egrep -v "#" Total_Up_${prefix}_deep.txt > Total_Up_${prefix}_deep_clean.txt
egrep -v "#" Total_Up_${prefix}_shallow.txt > Total_Up_${prefix}_shallow_clean.txt
egrep -v "#" Total_Up_${prefix}_all.txt > Total_Up_${prefix}_all_clean.txt
egrep -v "#" Total_Lo_${prefix}_deep.txt > Total_Lo_${prefix}_deep_clean.txt
egrep -v "#" Total_Lo_${prefix}_shallow.txt > Total_Lo_${prefix}_shallow_clean.txt
egrep -v "#" Total_Lo_${prefix}_all.txt > Total_Lo_${prefix}_all_clean.txt

# Finally remove the the non clean .txt files
rm -f Total_*${prefix}_deep.txt
rm -f Total_*${prefix}_shallow.txt
rm -f Total_*${prefix}_all.txt

