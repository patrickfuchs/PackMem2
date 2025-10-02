# New_PackMem
This repository contains the script for the new version of PackMem.\
This version presents a code written in python 3 and some new mechanisms such as the use of MDAnalysis package to work directly with the trajectory file or the possibility to see the defects around a protein inserted into the membrane.

-------------------------------------------------------------------------------

## Packages to install
These are the packages required to get this version of PackMem running :\
* argparse
* numpy
* MDAnalysis\

All of them are already put in the yaml file and the environment can be set up with:\
`conda env create -f PackMem_env.yml`

-------------------------------------------------------------------------------

## How to run the code
The code is simply launched with the command:\
`python [path]/New_PackMem/src/PackMem_prot.py -f [file.xtc] -s [file.gro] -l [lipid name] -r [path]/New_PackMem/data/vdw_radii_Charmm.txt -p [path]/New_PackMem/data/param_Charmm.txt -b [frame start] -e [frame end] -o [lipid name]`\
Note that you can also use the `-prot` option if you also want to analysis the packing defects near/far from the protein.\

If your computer has multiple cores, we advise you to parallelise the code using the `gen_packmem_launch_prot.py`:\
`python src/gen_packmem_launch_prot.py -f [file.xtc] -s [file.gro] -c [number of cores] -fm [number of frames] -l [lipid name]`\
You can also use the `-prot` option here.\
This script then gives you the command line(s) to copy/paste in the terminal to launch the actual analysis of the packing defects.\
The command line(s) given should look like this :\
`nice -19 python [path]src/PackMem_prot.py -f [file.xtc] -s [file.gro] -l [lipid name] -r [path]data/vdw_radii_Charmm.txt -p [path]data/param_Charmm.txt -b [frame start] -e [frame end] -o [lipid name] >& OUT_packmem0 &`\
With each command line corresponding to a core and the number of frames divided into these number of cores.\
Again, if you chose the `-prot` option ealier, then you will have the `-prot` option written in the command line(s).\

This script `PackMem_prot.py` will compute the packing defects, divided into 3 categories: deep, shallow and all packing defects.\

If you parallelised this script by using multiple cores, you then need to do one more action, that is to run the `Concatenate_PackMem_prot.sh` script.\
`bash src/Concatenate_PackMem_prot.sh -l [lipid name] -b [frame start] -e [frame stop]`\
Here again, you can also use the `-p` option if you ran the analysis of the packing defects near/far from the protein.\
This script will gather all the data file generated and put them into one file per type of defect (deep, shallow, all).\

## Analysis
In the end, you should have at least 9 files :\
* 3 for the Total defects - we mean by that the whole membrane defects - which are named `Total_[lipid name]_[type of defects]_clean.txt`
* 3 for the Upper leaflet defects - which are named `Total_Up_[lipid name]_[type of defects]_clean.txt`
* 3 for the Lower lealfet defects - which are named `Total_Lo_[lipid name]_[type of defects]_clean.txt`\

If you had the `-prot` option, you also have 3 or 6 files:
* 3 or 6 for the defects that are in near and far from the protein - which are named `Total_[leaflet]_[lipid name]_[type of defects]_prot.txt`

Once you have all the data files, you can run the plot-making script `Fit_plot.py`.\
`python analysis/Fit_plot.py -l [lipid name] -o [output name]`\
Here again, you can also use the `-prot` option if you ran the analysis of the packing defects near/far from the protein.\

-------------------------------------------------------------------------------

## Possible problems
If your lipid isn't present in the parameter files (data/param_Charmm.txt or vdw_radii_Charmm.txt), then you can add it on your own using the same convention.\
In param_Charmm.txt :\
* LIPID NAME_first atom three letter code (two frist letter and last one of the lipid name). Example : DMPC_N DMC
* LIPID NAME glycerol. Example : DMPC C2

In the vdw_radii_Charmm.txt:\
* LIPID NAME atom name  vdw radii n/a
You do this for every atoms in the lipid. If you don't know the van der Waals radii of each atoms, the script `Calculate_radii.py`is made to compute it for you for each atom.\
`python Calculate_radii.py -mol [file.itp of the lipid/molecule] -ff [forcefield.itp]`\

-------------------------------------------------------------------------------

## Contact
If there are any other problems, please contact me at :\
`maya.zygadlo@sorbonne-universite.fr`
