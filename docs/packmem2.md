# packmem2 Documentation

This script runs PackMem2's main loop. It reads the trajectory and computes if there are defects and their type between three:
* Shallow: superfluous defects 
* Deep: big defects
* All: the combination of the previous two defects

The output comprises of multiple .txt file:
* [outputname][Frame]\_[Up/Lo]\_[Shallow/Deep/All]\_result.txt\
        ## MatrixSize  9646  9801   # Membrane matrix size, Total matrix size\
        ## Total   51   582 11.41 6.034   # number of packing defects, total area of packing defects, average size, pourcent of membrane (Membrane matrix size)\
        1    1    77.00     2.00          # for each packing defect number, size, x_position (first cell), y_mean (first cell)\
        2    1    21.00     3.00\
        3    8    59.12     7.00
* Prot_[outputname][Frame]\_[Up/Lo]\_[Deep/Shallow/All].txt (if `-prot` flag)\
        1,far,2 # defects number, distance category, size\
        2,far,3\
        3,close,6
* [outputname]\_Total[Up/Lo]\_[Shallow/Deep/All].pdb\
            Matrix x,y with, for each cell, the value of "packing defect" (in B_factor column, the last column)\
            if 0 = Deep packing defects The value increases with the number of atoms in the cell.\
            if > 0 and < 1 = Shallow defect\
            if -1 = Edge\
            ATOM      3   H1 EDG     3       4.000  54.000  56.940  1.00 -1.00\
            ATOM      3   H1 MAT     3       4.000  54.000  56.940  1.00  0.00\
            ATOM      3   H1 MAT     3       4.000  54.000  56.940  1.00  0.00
* [outputname]\_Defect[Up/Lo]\_[Shallow/Deep/All].pdb
            Packing defects in pdb format:\
            the residue number corresponds to the different packing defects.\
            ATOM      5   H1 DEF     2       6.000  22.000  56.940  1.00  2.00\
            ATOM      5   H1 DEF     3       6.000  85.000  56.940  1.00  3.00\
            ATOM      5   H1 DEF     3       6.000  86.000  56.940  1.00  3.00\
            ATOM      6   H1 DEF     3       7.000  86.000  56.940  1.00  3.00

## Command line

The complete options for packmem2 are:
```
packmem2 -f [trajectory] -s [structure] -l [lipid name] -b [starting frame] -e [ending frame] -r [radii file] -p [parameter file] -d [distance in z] -n [index file] -o [output name] -od [output directory] [-prot] [-pdb] [-h]
```

## Options

Input files options:
- `-f` Trajectory file [.xtc/.trr/...]
- `-s` topology file [.gro/...]
- `-r` Van der Waals radii file (Default: `data/vdw_radii_Charmm.txt`)\
Other options are: `data/vdw_radii_Martini_old.txt`, `data/vdw_radii_Martini2.txt`, `data/vdw_radii_Martini2P.txt`, `data/vdw_radii_Martini3.txt`
- `-p` Parameter file (Default: `data/param_Charmm.txt`)\
Other option is: `data/param_Martini.txt`
- `-n` Index file (Optional)\
Index to the upper and lower leaflet

Input options:
- `-l` Lipid name(s)\
If multiple lipids, seperate them with `_`
- `-b` Starting Frame number (Default: 0)
- `-e` Ending Frame number 
- `-d` Distance to differenciate Deep from Shallow defects in Angstrom (Default: 1.0)
- `-prot` Protein in the simulation (Default: False)
- `-pdb` Create .pdb output (Default: False)

Ouput options:
- `-o` Output name (Default: output)
- `-od` Output directory (Default: ./)
