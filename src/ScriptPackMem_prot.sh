#!/bin/bash

# Script to launch PackMem analysis on multiple PDB files.
# P. Fuchs, R. Gautier, 2016

# set PackMem absolute path
# change this path according to your system 
PackMemPATH=/home/maya/Documents/tools/New_PackMem

usage(){
cat <<EOF

usage : $0 -f traj.xtc -s topo.gro -n prefix -b first_frame -e last_frame options

This script concatenates the PackMem analysis on multiple frame.

Options to specify input files:
    -f Trajectory file (.xtc)
    -s Topology file (.gro)
    -n prefix of the files
    -b The first frame of the PackMem analysis
    -e The last frame of the PackMem analysis

OPTIONS:
    -p Flag. If put, compute the default around the protein (Default : False)

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
if [ ${prot} == true ]
then
    # launch PackMem for the 3 types of packing defects
    python ${PackMemPATH}/src/PackMem_prot.py -f ${traj} -s ${topo} -l ${prefix} \
                            -r ${PackMemPATH}/data/vdw_radii_Charmm.txt \
                            -p ${PackMemPATH}/data/param_Charmm.txt \
                            -o ${prefix} -d 1.0 -b ${frame_start} -e ${frame_stop} -prot
else
    # launch PackMem for the 3 types of packing defects
    python ${PackMemPATH}/src/PackMem_prot.py -f ${traj} -s ${topo} -l ${prefix} \
                            -r ${PackMemPATH}/data/vdw_radii_Charmm.txt \
                            -p ${PackMemPATH}/data/param_Charmm.txt \
                            -o ${prefix} -d 1.0 -b ${frame_start} -e ${frame_stop}
fi
AFTER=`date +%s`
ELAPSED=$(((($AFTER-$BEFORE))/60))
echo "All the analysis ran for $ELAPSED minutes."
