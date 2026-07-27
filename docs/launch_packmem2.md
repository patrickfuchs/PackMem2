# Launch_packmem2 Documentation

This script runs the PackMem2 script, then the concatenation one and finaly, the analysis script.

The complete options for launch_packmem2 are:
```
launch_packmem2 -f [trajectory] -s [structure] -l [lipid name] -b [starting frame] -e [ending frame] -r [radii file] -p [parameter file] -c [core number] -d [distance in z] -n [index file] -prec [precision wanted] -lx [limit in x] -ly [limit in y] -o [output name] -od [output directory] [-prot] [-pdb] [-h]
```

## Options

Input files options:
- `-f` Trajectory file [.xtc/.trr/...]
- `-s` topology file [.gro/...]
- `-r` Van der Waals radii file (Default: `vdw_radii_Charmm.txt`)\
Other options are: `vdw_radii_Martini_old.txt`, `vdw_radii_Martini2.txt`, `vdw_radii_Martini2P.txt`, `vdw_radii_Martini3.txt`
- `-p` Parameter file (Default: `param_Charmm.txt`)\
Other option is: `param_Martini.txt`
- `-n` Index file (Optional)\
Index to the upper and lower leaflet

Input options:
- `-l` Lipid name(s)\
If multiple lipids, seperate them with `_`
- `-b` Starting Frame number (Default: 0)
- `-e` Ending Frame number 
- `-c` Number of cores given for the analysis (Default: maximum cores on the computer - 2) 
- `-d` Distance to differenciate Deep from Shallow defects in Angstrom (Default: 1.0)
- `-prec` Precision for defect constants in the output (Default: 2)
- `-lx` Lowest defect area considered for the fit (Default: 15)
- `-ly` Lowest probability considered for the fit (Default: 1e-4)
- `-prot` Protein in the simulation (Default: False)
- `-pdb` Create .pdb output (Default: False)

Ouput options:
- `-o` Output name (Default: output)
- `-od` Output directory (Default: ./)
