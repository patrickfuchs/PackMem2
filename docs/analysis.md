# analysis Documentation

This script runs the analysis of the packmem2 results.

The output is a single .pdf file:
* `Res_membrane.pdf` that contains all the plots of the defects analysis.\
        In order: the fitted decay plots, the block averaging bar plots and the final bar plots\
        For each defect type

## Command line

The complete options for analysis are:
```
analysis -p [precision wanted] -lx [limit in x] -ly [limit in y] -o [output name] -od [output directory] [-prot] [-h]
```

## Options

Input options:
- `-p` Precision for defect constants in the output (Default: 2)
- `-lx` Lowest defect area considered for the fit (Default: 15)
- `-ly` Lowest probability considered for the fit (Default: 1e-4)
- `-prot` Protein in the simulation (Default: False)

Ouput options:
- `-o` Output name (Default: Res_membrane)
- `-od` Output directory (Default: ./)
