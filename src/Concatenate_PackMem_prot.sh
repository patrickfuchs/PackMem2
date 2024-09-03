#!/bin/bash

usage(){
cat <<EOF

usage : $0 -n prefix -b first_frame -e last_frame options

This script concatenates the PackMem analysis on multiple frame.

Options to specify input files:
    -n prefix of the files
    -b The first frame of the PackMem analysis
    -e The last frame of the PackMem analysis

OPTIONS:
    -p If there is a protein (Dafult : false)

EOF
}

prot=false

while getopts "hn:b:e:p" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        n)
            prefix=$OPTARG
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
    # rm -f ${prefix}${pdbnum}_Up_Deep_result.txt
    # rm -f ${prefix}${pdbnum}_Lo_Deep_result.txt
    # rm -f ${prefix}${pdbnum}_Up_Shallow_result.txt
    # rm -f ${prefix}${pdbnum}_Lo_Shallow_result.txt
    # rm -f ${prefix}${pdbnum}_Up_All_result.txt
    # rm -f ${prefix}${pdbnum}_Lo_All_result.txt
done

if [ ${prot} = true ]
then
    # Accumulate packing defects of the current frame in Total_*.txt files per leaflet
     cat  ${prefix}${pdbnum}_Up_All_prot.txt >> Total_Up_${prefix}_all_prot.txt
     cat  ${prefix}${pdbnum}_Up_Deep_prot.txt >> Total_Up_${prefix}_deep_prot.txt
     cat  ${prefix}${pdbnum}_Up_Shallow_prot.txt >> Total_Up_${prefix}_shallow_prot.txt

    # we no longer need the defects of the current frame
     rm -f ${prefix}${pdbnum}_Up_All_prot.txt
     rm -f ${prefix}${pdbnum}_Up_Deep_prot.txt
     rm -f ${prefix}${pdbnum}_Up_Shallow_prot.txt
fi

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
