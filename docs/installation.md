# Installation

## Requirements and compatibility

PackMem2 requires at least Python 3.10 and needs the following dependencies :

- argparse
- matplotlib
- mdanalysis
- networkx
- numpy
- pandas

## Installation

We highly recommand the use of [uv](https://docs.astral.sh/uv/) to manage all dependencies and the project environment. First, [install uv](https://docs.astral.sh/uv/getting-started/installation/). Once done, you can simply use:

```bash
uv tool install packmem2
```

This uv command will install PackMem2 and all its dependencies in your `~/.local/bin` and `~/.local/share/uv/tools/packmem2`. Since `~/.local/bin` should be in your `PATH` environment variable, all PackMem2 scripts will be available right away. To test if it works, just invoke:

```bash
$ launch_packmem2
usage: launch_packmem2 [-h] -f TRAJ -s TOPO -l LIPID [-b START] [-e END] [-r RADIIFILE]
                       [-p PARAMFILE] [-o OUTPUTNAME] [-od OUTPUT_DIR] [-d DIST_SUPPL_Z]
                       [-n INDEXFILE] [-pdb] [-prot] [-c CORES] [-nb NB_BLOCK] [-prec PRECISION]
                       [-lx LIMX] [-ly LIMY]
```

If you get the error `launch_packmem: command not found`, you can add `~/.local/bin` to your PATH:

```bash
export PATH="$PATH:~/.local/bin"
```

## Installation with pip

TODO

## Installation of a development version

Clone the repository:

```bash
git clone https://github.com/patrickfuchs/PackMem2.git
```

Then move to PackMem2 directory:

```bash
cd PackMem2
```

We use [uv](https://docs.astral.sh/uv/) to manage all dependencies. Sync the environment:

```
uv sync
```

Then all the PackMem2 scripts (`launch_packmem2`, `packmem2`, etc.) will be available with:

```bash
uv run launch_packmem2 [options]
```
