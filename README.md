# PackMem2

[PackMem](https://packmem.ipmc.cnrs.fr/) ([link to the github repo](https://github.com/rogautier/packmem)) is a computational tool that computes the distributions of packing defects in lipid bilayers from molecular dynamics trajectories. The present repository contains the new version of PackMem called **PackMem2**.

PackMem2 has been completely refactorized compared to the original PackMem (notably with the use of Numpy, MDAnalysis, and the multiprocessing modules), which heavily improved performances. PackMem2 also implements some new features, notably computing packing defects around a protein inserted into the membrane.

-------------------------------------------------------------------------------

## Packages to install
These are the packages required to get this version of PackMem running :\
* argparse
* numpy
* pandas
* matplotlib
* MDAnalysis\

All of them are already put in the .yml file and the environment can be set up with:\
`conda env create -f PackMem_env.yml`

-------------------------------------------------------------------------------

## How to run the code
The code is simply launched with the command:\
`python [path]/PackMem2/src/packmem2/launch_packmem2.py -f [file.xtc] -s [file.gro] -c [number of cores] -fm [number of frames] -l [lipid name]`\
You can also choose the parameter files (depending on your simulation) between martini and Charmm (knowing that Charmm is the default one) with the `-p` ans `-r` options.\
Note that you can also use the `-prot` option if you also want to analysis the packing defects near/far from the protein.\

This script launches the actual PackMem script, then the concatenation script and finaly, the analysis script.

Of course, if you encounter problems along the way, each scripts are individually launchable with their respective commands:
`python [path]/PackMem2/src/packmem2/packmem2.py -f [file.xtc] -s [file.gro] -l [lipid name] -r [path]/PackMem2/data/vdw_radii_Charmm.txt -p [path]/PackMem2/data/param_Charmm.txt -b [frame start] -e [frame end] -o [lipid name]`\
Here again, the `-prot` option is usable.\
This script `packmem2.py` will compute the packing defects, divided into 3 categories: deep, shallow and all packing defects.\

`python [path]/PackMem2/src/packmem2/concatenate.py -l [lipid name] -b [frame start] -e [frame stop]`\
Here again, you can also use the `-prot` option if you ran the analysis of the packing defects near/far from the protein.\
This script will gather all the data file generated and put them into one file per type of defect (deep, shallow, all).\

`python [path]/PackMem2/src/packmem2/analysis.py -o [output name]`\
Here again, you can also use the `-prot` option if you ran the analysis of the packing defects near/far from the protein.\

## Results
In the end, you should have at least 9 files that were used for the analysis :\
* 3 for the Total defects - we mean by that the whole membrane defects - which are named `Total_[type of defects].csv`
* 3 for the Upper leaflet defects - which are named `Total_Up_[type of defects].csv`
* 3 for the Lower lealfet defects - which are named `Total_Lo_[type of defects]csv`\

If you had the `-prot` option, you also have 3 or 6 files used for the analysis:
* 3 or 6 for the defects that are in near and far from the protein - which are named `Total_[leaflet]_[type of defects]_prot.csv`\

And the result file:\
* `Res_membrane.pdf` that contains all the plots of the defects analysis



-------------------------------------------------------------------------------

## Possible problems

If your lipid isn't present in the parameter files (in the data folder), then you can add it on your own using the same convention.\
In param_Charmm.txt or param_Martini.txt:\
* LIPID NAME glycerol / central atom. Example : DMPC C2

In the vdw_radii_Charmm.txt or vdw_radii_Martini*.txt:\
* LIPID NAME atom name  vdw radii n/a
You do this for every atoms in the lipid. If you don't know the van der Waals radii of each atoms, the script `calculate_radii.py`is made to compute it for you for each atom.\
`python [path]/PackMem2/src/packmem2/calculate_radii.py -mol [file.itp of the lipid/molecule] -ff [forcefield.itp]`\
Note that you can use the `-martini` or `-martini3` option to precise how to compute the radii.\

-------------------------------------------------------------------------------

## Contact

If there are any other problems, please contact: `maya.zygadlo@sorbonne-universite.fr`


## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md)
