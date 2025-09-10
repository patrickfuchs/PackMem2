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
In order to have the script running smoothly, you first need to launch the script `gen_packmem_launch_prot.py` :\
`python src/gen_packmem_launch_prot.py -f [file.xtc] -s [file.gro] -c [number of cores] -fm [number of frames] -n [lipid name]`\
Note that you can also use the `-prot` option if you also want to analysis the packing defects near/far from the protein.\

This script then gives you the command line(s) to copy/paste in the terminal to launch the actual analysis of the packing defects.\
The command line(s) given should look like this :\
`nice -19 bash src/ScriptPackMem_prot.sh -f [file.xtc] -s [file.gro] -n [lipid name] -b [frame start] -e [frame end] >& OUT_packmem0 &`\
With each command line corresponding to a core and the number of frames divided into these number of cores.\
Again, if you chose the `-prot` option ealier, then you will have the `-p` option written in the command line(s).\
Just be careful to change the path of the `.sh` script if you want to launch it from elsewhere.\

This script `ScriptPackMem_prot.sh` launches the actual PackMem script `PackMem_prot.py` that will compute the packing defects, first for the deep defects, then the shallow ons and finaly all of them sequentially.\

If you parallelised this script by using multiple cores, you then need to do one more action, that is to run the `Concatenate_PackMem_prot.sh` script.\
`bash src/Concatenate_PackMem_prot.sh -n [lipid name] -b [frame start] -e [frame stop]`\
Here again, you can also use the `-p` option if you ran the analysis of the packing defects near / far from the protein.\
This script will gather all the data file generated and put them into one file per type of defect (deep, shallow, all).\
In the end, you should have 9 files :\
* 3 for the Total defects - we mean by that the whole membrane defects - which are named `Total_[lipid name]_[type of defects]_clean.txt`
* 3 for the Upper leaflet defects - which are named `Total_Up_[lipid name]_[type of defects]_clean.txt`
* 3 for the Lower lealfet defects - which are named `Total_Lo_[lipid name]_[type of defects]_clean.txt`

Once you have all the data files, you can run the plot-making script `Fit_plot.py` or its R equivalent `Script_fit_and_plot.R`.\
`python analysis/Fit_plot.py -f Total_[lipid name]`\
OR\
`R --vanilla < Script_fit_and_plot.R`\

If you studied also the packing defacts near/far from the protein, you should run `Script_fit_and_plot_prot.R`\
`R --vanilla < Script_fit_and_plot_prot.R`\
For the R scripts, be careful to change the variables inside the script.\

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
