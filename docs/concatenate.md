# concatenate Documentation

This script concatenates packmem2 results. It reads the output .txt files and concatenates the date into 3 files file by defect type (Shallow / Deep / All).

The output comprises of multiple .csv file:
* Total\_[Up/Lo/-]\_[type of defects].csv\
        1,1,77.00,2.00 # packing defect number, size, x_position (first cell), y_mean (first cell)\
        2,1,21.00,3.00\
        3,8,59.12,7.00
* Total\_[Up/Lo/-]\_[type of defects]\_prot.csv\ (if `-prot` flag)\
        1,far,2 # defects number, distance category, size\
        2,far,3\
        3,close,6

## Command line

```
concatenate -l [lipid name] -b [starting frame] -e [ending frame] -od [output directory] [-prot] [-h]
```

## Options

Input options:
- `-l` Prefix (usually the lipid name(s))\
If multiple lipids, seperate them with `_`
- `-b` Starting Frame number (Default: 0)
- `-e` Ending Frame number 
- `-prot` Protein in the simulation (Default: False)

Ouput options:
- `-od` Output directory (Default: ./)
