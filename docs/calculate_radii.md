# calculate_radii Documentation

This script computes the van der Waals radii from .itp files.

The output is written in an output file:
* `vdw_out.txt` that contains the van der Waals radii for the molecule file given.

## Command line

The complete options for analysis are:
```
calculate_radii -m [molecule file] -f [forcefield file] -ff [forcefield] -o [output name] [-h]
```

## Options

Input files options:
- `-m` molecule parameter file [.itp]
- `-f` force field parameter file [.itp]
- `-ff` force field of the given files [charmm, martini, martini3]

Ouput options:
- `-o` Output name (Default: vdw_out)

## Trouble shooting

This script is here for helping compute the radius of each atom. It adds the information of the aliphathic nature of the atom. This should always be double checked as the order the atom appears it the molecule.itp file may vary !