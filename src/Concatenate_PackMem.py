"""Script to concatenate the PackMem analysis on multiple PDB files."""

import argparse
import os
import pandas as pd

def get_args():
    """
    Get the arguments for the script and check that the inputfiles are valid.

    --------------------
    OUTPUT
    parser.parse_args
        Contains all the arguments for the script
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', action='store', dest='prefix', 
                        help='The prefix of PackMem output file names.')
    parser.add_argument('-b', action = 'store', dest = 'start', type=int,
                        help = 'The first frame.')
    parser.add_argument('-e', action = 'store', dest = 'end',
                        help = 'The number of frames.',
                        type = int)
    parser.add_argument('-prot', action='store_true', dest='prot', 
                        help='Put if you want to also concatenante Protein files.')
    args = parser.parse_args()

    return args

def concat_files(prefix, suffix):
    """
    Concatenate .txt files from PackMem_prot.py

    --------------------
    INPUT
    prefix: str
        The name of the lipid usually but could just be the common name of the outputfiles
    suffix: str
        The end name of the .txt files
    
    --------------------
    OUTPUT
    pandas DataFrame
        Contains the informations in the .txt files
    """
    file_concat = pd.read_csv(f"{prefix}0{suffix}", sep=r'\s+', header=None, skiprows=[0,1])
    for pdbnum in range(args.start+1, args.end+1):
        file = pd.read_csv(f"{prefix}{pdbnum}{suffix}", sep=r'\s+', header=None, skiprows=[0,1])
        file_concat = pd.concat([file_concat, file], axis=0, ignore_index=True)
    
    return file_concat

def concat_files_prot(prefix, suffix):
    """
    Concatenate .txt files from PackMem_prot.py for the protein

    --------------------
    INPUT
    prefix: str
        The name of the lipid usually but could just be the common name of the outputfiles
    suffix: str
        The end name of the .txt files
    
    --------------------
    OUTPUT
    pandas DataFrame
        Contains the informations in the .txt files
    """
    file_concat = pd.read_csv(f"{prefix}0{suffix}", header=None)
    for pdbnum in range(args.start+1, args.end+1):
        file = pd.read_csv(f"{prefix}{pdbnum}{suffix}", header=None)
        file_concat = pd.concat([file_concat, file], axis=0, ignore_index=True)
    
    return file_concat

if __name__=="__main__":
    args = get_args()

    Total_Up_Deep = concat_files(args.prefix, "_Up_Deep_result.txt")
    Total_Lo_Deep = concat_files(args.prefix, "_Lo_Deep_result.txt")
    Total_Up_All = concat_files(args.prefix, "_Up_All_result.txt")
    Total_Lo_All = concat_files(args.prefix, "_Lo_All_result.txt")
    Total_Up_Shallow = concat_files(args.prefix, "_Up_Shallow_result.txt")
    Total_Lo_Shallow = concat_files(args.prefix, "_Lo_Shallow_result.txt")

    if args.prot:
        if os.path.isfile(f"Prot_{args.prefix}0_Up_All.txt"):
            Total_Up_Deep_prot = concat_files_prot(f"Prot_{args.prefix}", "_Up_Deep.txt")
            Total_Up_All_prot = concat_files_prot(f"Prot_{args.prefix}", "_Up_All.txt")
            Total_Up_Shallow_prot = concat_files_prot(f"Prot_{args.prefix}", "_Up_Shallow.txt")
            # Save files
            Total_Up_Deep_prot.to_csv("Total_Up_Deep_prot.csv", header=False, index=False)
            Total_Up_All_prot.to_csv("Total_Up_All_prot.csv", header=False, index=False)
            Total_Up_Shallow_prot.to_csv("Total_Up_Shallow_prot.csv", header=False, index=False)
            # Remove the files
            for pdbnum in range(args.start, args.end+1):
                os.remove(f"Prot_{args.prefix}{pdbnum}_Up_Deep.txt")
                os.remove(f"Prot_{args.prefix}{pdbnum}_Up_All.txt")
                os.remove(f"Prot_{args.prefix}{pdbnum}_Up_Shallow.txt")
        elif os.path.isfile(f"Prot_{args.prefix}0_Lo_All.txt"):
            Total_Lo_Deep_prot = concat_files_prot(f"Prot_{args.prefix}", "_Lo_Deep.txt")
            Total_Lo_All_prot = concat_files_prot(f"Prot_{args.prefix}", "_Lo_All.txt")
            Total_Lo_Shallow_prot = concat_files_prot(f"Prot_{args.prefix}", "_Lo_Shallow.txt")
            # Save files
            Total_Lo_Deep_prot.to_csv("Total_Lo_Deep_prot.csv", header=False, index=False)
            Total_Lo_All_prot.to_csv("Total_Up_All_prot.csv", header=False, index=False)
            Total_Lo_Shallow_prot.to_csv("Total_Up_Shallow_prot.csv", header=False, index=False)
            # Remove the files
            for pdbnum in range(args.start, args.end+1):
                os.remove(f"Prot_{args.prefix}{pdbnum}_Lo_Deep.txt")
                os.remove(f"Prot_{args.prefix}{pdbnum}_Lo_All.txt")
                os.remove(f"Prot_{args.prefix}{pdbnum}_Lo_Shallow.txt")
        
    # Concatenate the lealfets' results
    Total_Deep =  pd.concat([Total_Up_Deep, Total_Lo_Deep], axis=0, ignore_index=True)
    Total_All =  pd.concat([Total_Up_All, Total_Lo_All], axis=0, ignore_index=True)
    Total_Shallow =  pd.concat([Total_Up_Shallow, Total_Lo_Shallow], axis=0, ignore_index=True)

    # Save all files to .csv
    Total_Up_Deep.to_csv("Total_Up_Deep.csv", header=False, index=False)
    Total_Lo_Deep.to_csv("Total_Lo_Deep.csv", header=False, index=False)
    Total_Up_All.to_csv("Total_Up_All.csv", header=False, index=False)
    Total_Lo_All.to_csv("Total_Lo_All.csv", header=False, index=False)
    Total_Up_Shallow.to_csv("Total_Up_Shallow.csv", header=False, index=False)
    Total_Lo_Shallow.to_csv("Total_Lo_Shallow.csv", header=False, index=False)

    Total_Deep.to_csv("Total_Deep.csv", header=False, index=False)
    Total_All.to_csv("Total_All.csv", header=False, index=False)
    Total_Shallow.to_csv("Total_Shallow.csv", header=False, index=False)

    # Remove the files
    for pdbnum in range(args.start, args.end+1):
        os.remove(f"{args.prefix}{pdbnum}_Up_Deep_result.txt")
        os.remove(f"{args.prefix}{pdbnum}_Lo_Deep_result.txt")
        os.remove(f"{args.prefix}{pdbnum}_Up_All_result.txt")
        os.remove(f"{args.prefix}{pdbnum}_Lo_All_result.txt")
        os.remove(f"{args.prefix}{pdbnum}_Up_Shallow_result.txt")
        os.remove(f"{args.prefix}{pdbnum}_Lo_Shallow_result.txt")
