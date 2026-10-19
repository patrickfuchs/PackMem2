# PackMem2

**PackMem2** is a computational tool that computes lipid packing defects and their distribution from of packing defects in lipid bilayers from a molecular dynamics trajectory. **PackMem2** is an upgrade of the former [PackMem](https://packmem.ipmc.cnrs.fr/). It has been completely refactorized compared to the original PackMem (notably with the use of Numpy, MDAnalysis, and the multiprocessing modules), which heavily improved performances. **PackMem2** also implements some new features, notably computing packing defects around a protein inserted into the membrane.

## Installation

**PackMem2** is intended to run under Linux systems. We highly recommand the use of [uv](https://docs.astral.sh/uv/) to manage all dependencies and the project environment. First, [install uv](https://docs.astral.sh/uv/getting-started/installation/). Once done, you can simply use:

```bash
uv tool install packmem2
```

This command installs all the scripts and dependencies in your `~/.local` directory. Since `~/.local/bin` is generally on your PATH environment variable under Linux systems, all PackMem2 scripts are directly accessible on the command line. More on installation [here]().

## Running PackMem2

The main script is called `launch_packmem2` which will do the full analysis. A minimal command on a pure lipid bilayer (say made of POPC lipids) simulated with GROMACS would look like this:

```bash
launch_packmem2 -f trajectory.xtc -s topol.tpr -l POPC
```

This command will loop over the trajectory, calculate all types of defects (deep, shallow, all) on each frame, concatenate defects from all frames, perform the analysis and write results into csv files, and generate a pdf with multiple plots. You can see the full explanation of this commend and its many options in the [doc]().

`launch_packmem2` is a wrapper which uses some other scripts that can be launched individually: 

- [packmem2](): computation of packing defects on a single frame;
- [concatenate](): concatenation of packing defects from all frames;
- [analysis](): analysis of results and generation of plots.

`launch_packmem2` produces the following files:

- `Total_[type_of_defects].csv`: all defects of the trajectory whatever the leaflet.
- `Total_Up_[type_of_defects].csv`: all defects of the trajectory of the upper leaflet.
- `Total_Lo_[type_of_defects].csv`: all defects of the trajectory of the lower leaflet.

where `type_of_defects` is within `(deep|shallow|all)`.

Each line of these csv files contains information on a single defect: defect number, defect area in units of $\AA^2$, an approximate `x` and `y` coordinate where the defect stands.

- `Res_membrane.pdf`: contains all plots from which the packing defect constants are computed.

PackMem2 allows the analysis of simulation containing a membrane protein or peptide with the `-prot` flag. More on this feature [here]().

## Documentation

The full documentation is available at [Read the Docs](https://readthedocs.org/projects/packmem2).

## Contributors

- Maya Zygadlo
- Pierre Poulain
- Patrick Fuchs

## Licence

buildH is licensed under the [BSD License](https://github.com/patrickfuchs/buildH/blob/master/LICENSE.txt).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md)

## Citing PackMem2

If you use PackMem2 for your research, please cite:

```txt
Ref to come
```
