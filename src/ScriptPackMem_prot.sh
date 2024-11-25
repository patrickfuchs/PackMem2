#!/bin/bash

# Script to launch PackMem analysis on multiple PDB files.
# P. Fuchs, R. Gautier, 2016

# set PackMem absolute path
# change this path according to your system 
PackMemPATH=/home/maya/Documents/tools/New_PackMem/

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

while getopts "hf:s:n:b:e:p" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        f)
            traj=$OPTARG
            ;;
        s)
            topo=$OPTARG
            ;;
        n)
            prefix=$OPTARG
            ;;
        b)
            frame_start=$OPTARG
            ;;
        e)
            frame_stop=$OPTARG
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


###################
# No change below #
###################
BEFORE=`date +%s`
# print counter to screen
#echo "$(date): PackMem running on pdb/${prefix}${pdbnum}.pdb"
# launch PackMem for the 3 types of packing defects
python ${PackMemPATH}/src/PackMem_prot.py -f ${traj} -s ${topo} \
                            -r ${PackMemPATH}/data/vdw_radii_Charmm.txt \
                            -p ${PackMemPATH}/data/param_Charmm.txt \
                            -o ${prefix} -d 1.0 -t deep -b ${frame_start} -e ${frame_stop}
AFTER=`date +%s`
ELAPSED=$(((($AFTER-$BEFORE))/60))
echo "Deep analysis ran for $ELAPSED minutes."
BEFORE=`date +%s`
python ${PackMemPATH}/src/PackMem_prot.py -f ${traj} -s ${topo} \
                            -r ${PackMemPATH}/data/vdw_radii_Charmm.txt \
                            -p ${PackMemPATH}/data/param_Charmm.txt \
                            -o ${prefix} -d 1.0 -t shallow -b ${frame_start} -e ${frame_stop}
AFTER=`date +%s`
ELAPSED=$(((($AFTER-$BEFORE))/60))
echo "Shallow analysis ran for $ELAPSED minutes."
BEFORE=`date +%s`
python ${PackMemPATH}/src/PackMem_prot.py -f ${traj} -s ${topo} \
                            -r ${PackMemPATH}/data/vdw_radii_Charmm.txt \
                            -p ${PackMemPATH}/data/param_Charmm.txt \
                            -o ${prefix} -d 1.0 -t all -b ${frame_start} -e ${frame_stop}
AFTER=`date +%s`
ELAPSED=$(((($AFTER-$BEFORE))/60))
echo "All analysis ran for $ELAPSED minutes."
