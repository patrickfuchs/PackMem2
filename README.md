# PackMem2

[PackMem](https://packmem.ipmc.cnrs.fr/) ([link to the github repo](https://github.com/rogautier/packmem)) is a computational tool that computes the distributions of packing defects in lipid bilayers from molecular dynamics trajectories. The present repository contains the new version of PackMem called **PackMem2**.

PackMem2 has been completely refactorized compared to the original PackMem (notably with the use of Numpy, MDAnalysis, and the multiprocessing modules), which heavily improved performances. PackMem2 also implements some new features, notably computing packing defects around a protein inserted into the membrane.

-------------------------------------------------------------------------------

## Setup environment

We use [uv](https://docs.astral.sh/uv/) to manage all dependencies and the project environment.

To get PackMem2 up and running, follow the next steps.\
Clone this repository:
```
git clone https://github.com/patrickfuchs/PackMem2.git
```

Then move to PackMem2 directory:
```
cd PackMem2
```

Create the environment:
```
uv sync
```

-------------------------------------------------------------------------------

## How to run the code

The package is launched with the minimal command:
```
uv run launch_packmem2 -f [trajectory] -s [structure] -l [lipid name] -e [ending frame]
```

You can see the full explanation of the options in [launch_packmem2](docs/launch_packmem2.md).

This script runs the actual PackMem2, the concatenation of the output files and the analysis.

If any problems is encountered along the way, each scripts can be individually run with their respective commands:
* [packmem2](docs/packmem2.md)
* [concatenate](docs/concatenate.md)
* [analysis](docs/analysis.md)

## Results

In the end, you should have the following 9 basic files:
* 3 for the Total defects - we mean by that the whole membrane defects - which are named `Total_[type of defects].csv`
* 3 for the Upper leaflet defects - which are named `Total_Up_[type of defects].csv`
* 3 for the Lower lealfet defects - which are named `Total_Lo_[type of defects]csv`

If you used the `-prot` flag, 3 or 6 more files will be in your directory:
* 3 or 6 for the defects that are near and far from the protein - which are named `Total_[leaflet]_[type of defects]_prot.csv`

And the result .pdf file:
* `Res_membrane.pdf` that contains all the plots of the defects analysis.

-------------------------------------------------------------------------------

## Possible problems

If your lipid isn't present in the parameter files (param and vdw_radii files in the `data` folder), then you can add it on your own using the same convention.

In the [CHARMM parameter file](data/param_Charmm.txt) or the [MARTINI parameter file](data//param_Martini.txt):
* LIPID NAME glycerol / central atom. Example : DMPC C2

And in the [CHARMM radii files](data/vdw_radii_Charmm.txt) or [MARTINI radii files](data/vdw_radii_Martini_old.txt) - Note that there are multiple for the MARTINI force field:
* LIPID NAME atom name  vdw radii n/a
This is done for every lipid atoms.

If you don't know the van der Waals radii of each atoms, run the following script:
```
uv run calculate_radii -mol [lipid.itp] -ff [forcefield.itp]
```
Here again see its dedicated documentation.\
[calculate_radii](docs/calculate_radii.md)

-------------------------------------------------------------------------------

## Contact

If there are any other problems, please contact: `maya.zygadlo@sorbonne-universite.fr`


## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md)
