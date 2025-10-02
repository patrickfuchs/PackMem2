#!/bin/bash

# Script to concatenate the PackMem analysis on multiple PDB files.

usage(){
cat <<EOF

usage : $0 -l lipid_name -b first_frame -e last_frame options

This script concatenates the PackMem analysis on multiple frame.

Options to specify input files:
    -l lipid_name of the files
    -b The first frame of the PackMem analysis
    -e The last frame of the PackMem analysis

OPTIONS:
    -p If there is a protein (Default : false)

EOF
}

prot=false
first_frame=0

while getopts "hl:b:e:p" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        l)
            lipid_name=$OPTARG
            ;;
        b)
            first_frame=$OPTARG
            ;;
        e)
            last_frame=$OPTARG
            ;;
        p)
            prot=true
            ;;
        ?)
            usage
            exit
            ;;
    esac
done


# The following Total_*.txt files will contain the statistics for each type of
# packing defects accumulated over all frames (one file per membrane leaflet)
# Delete those files if they already exist before launching the main loop
# (the -f flag avoids an error if the file doesn't exist)
rm -f Total_Up_${lipid_name}_Deep.txt
rm -f Total_Lo_${lipid_name}_Deep.txt
rm -f Total_Up_${lipid_name}_Shallow.txt
rm -f Total_Lo_${lipid_name}_Shallow.txt
rm -f Total_Up_${lipid_name}_All.txt
rm -f Total_Lo_${lipid_name}_All.txt



for pdbnum in $(seq ${first_frame} ${last_frame})
do
    # Accumulate packing defects of the current frame in Total_*.txt files per leaflet
    cat  ${lipid_name}${pdbnum}_Up_Deep_result.txt >> Total_Up_${lipid_name}_Deep.txt
    cat  ${lipid_name}${pdbnum}_Lo_Deep_result.txt >> Total_Lo_${lipid_name}_Deep.txt
    cat  ${lipid_name}${pdbnum}_Up_Shallow_result.txt >> Total_Up_${lipid_name}_Shallow.txt
    cat  ${lipid_name}${pdbnum}_Lo_Shallow_result.txt >> Total_Lo_${lipid_name}_Shallow.txt
    cat  ${lipid_name}${pdbnum}_Up_All_result.txt >> Total_Up_${lipid_name}_All.txt
    cat  ${lipid_name}${pdbnum}_Lo_All_result.txt >> Total_Lo_${lipid_name}_All.txt

    # we no longer need the defects of the current frame
    rm -f ${lipid_name}${pdbnum}_Up_Deep_result.txt
    rm -f ${lipid_name}${pdbnum}_Lo_Deep_result.txt
    rm -f ${lipid_name}${pdbnum}_Up_Shallow_result.txt
    rm -f ${lipid_name}${pdbnum}_Lo_Shallow_result.txt
    rm -f ${lipid_name}${pdbnum}_Up_All_result.txt
    rm -f ${lipid_name}${pdbnum}_Lo_All_result.txt

    if [ ${prot} = true ]
    then
        # Accumulate packing defects of the current frame in Total_*.txt files per leaflet
        cat  Prot_${lipid_name}${pdbnum}_Up_All.txt >> Total_Up_${lipid_name}_All_prot.txt
        cat  Prot_${lipid_name}${pdbnum}_Up_Deep.txt >> Total_Up_${lipid_name}_Deep_prot.txt
        cat  Prot_${lipid_name}${pdbnum}_Up_Shallow.txt >> Total_Up_${lipid_name}_Shallow_prot.txt

        # we no longer need the defects of the current frame
        rm -f Prot_${lipid_name}${pdbnum}_Up_All_prot.txt
        rm -f Prot_${lipid_name}${pdbnum}_Up_Deep_prot.txt
        rm -f Prot_${lipid_name}${pdbnum}_Up_Shallow_prot.txt
    fi
done

# Concatenate the packing defects statistics of the lower and upper leaflets
cat Total_Up_${lipid_name}_Deep.txt Total_Lo_${lipid_name}_Deep.txt > Total_${lipid_name}_Deep.txt
cat Total_Up_${lipid_name}_Shallow.txt Total_Lo_${lipid_name}_Shallow.txt > Total_${lipid_name}_Shallow.txt
cat Total_Up_${lipid_name}_All.txt Total_Lo_${lipid_name}_All.txt > Total_${lipid_name}_All.txt

# Clean txt file for python script
egrep -v "#" Total_${lipid_name}_Deep.txt > Total_${lipid_name}_Deep_clean.txt
egrep -v "#" Total_${lipid_name}_Shallow.txt > Total_${lipid_name}_Shallow_clean.txt
egrep -v "#" Total_${lipid_name}_All.txt > Total_${lipid_name}_All_clean.txt

# Clean leaflet txt file for python script
egrep -v "#" Total_Up_${lipid_name}_Deep.txt > Total_Up_${lipid_name}_Deep_clean.txt
egrep -v "#" Total_Up_${lipid_name}_Shallow.txt > Total_Up_${lipid_name}_Shallow_clean.txt
egrep -v "#" Total_Up_${lipid_name}_All.txt > Total_Up_${lipid_name}_All_clean.txt
egrep -v "#" Total_Lo_${lipid_name}_Deep.txt > Total_Lo_${lipid_name}_Deep_clean.txt
egrep -v "#" Total_Lo_${lipid_name}_Shallow.txt > Total_Lo_${lipid_name}_Shallow_clean.txt
egrep -v "#" Total_Lo_${lipid_name}_All.txt > Total_Lo_${lipid_name}_All_clean.txt

# Finally remove the the non clean .txt files
rm -f Total_*${lipid_name}_Deep.txt
rm -f Total_*${lipid_name}_Shallow.txt
rm -f Total_*${lipid_name}_All.txt

