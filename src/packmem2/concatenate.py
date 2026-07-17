"""Script to concatenate the PackMem analysis on multiple PDB files."""

import argparse
import os
import pandas as pd
from pathlib import Path

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
                        required=True,
                        help='The prefix of PackMem output file names.')
    parser.add_argument('-od', action='store', dest='output_dir',
                        default = './',
                        help = 'Name for output directory (default: ./)')
    parser.add_argument('-b', action = 'store', dest = 'start', type=int,
                        required=True,
                        help = 'The first frame.')
    parser.add_argument('-e', action = 'store', dest = 'end', type = int,
                        required=True,
                        help = 'The number of frames.')
    parser.add_argument('-prot', action='store_true', dest='protein', 
                        help='Put if you want to also concatenante Protein files.')
    args = parser.parse_args()

    return args

def check_file(filename):
    """
    Check if there is content in the file.

    --------------------
    INPUT
    filename: str
        The name the input file

    --------------------
    OUTPUT
    Boolean
        If the files has content
    """
    with open(filename, 'r') as file_in:
        count = 0
        for line in file_in:
            count+=1
            if count == 3:
                return True
        return False

def start_file(prefix, suffix, start, end):
    """
    Find the first file to read

    --------------------
    INPUT
    prefix: string
        The begining of the file to be read
    suffix: string
        The ending of the file to be read
    start: int
        The minimum frame number in the filenames
    end: int
        The maximum frame number in the filenames
    
    --------------------
    OUTPUT
    string
        The first filename to be read
    int
        The frame number of the first filename to be read
    """
    for i in range(start, end+1):
        if check_file(f"{prefix}{i}{suffix}"):
            filename = f"{prefix}{i}{suffix}"
            return filename, i

def concat_files(prefix, suffix, start, end):
    """
    Concatenate .txt files from PackMem_prot.py

    --------------------
    INPUT
    prefix: str
        The name of the lipid usually but could just be the common name of the outputfiles
    suffix: str
        The end name of the .txt files
    start: int
        The minimum frame number in the filenames
    end: int
        The maximum frame number in the filenames
    
    --------------------
    OUTPUT
    pandas DataFrame
        Contains the informations in the .txt files
    """
    filename, start = start_file(prefix, suffix, start, end)
    file_concat = pd.read_csv(filename, sep=r'\s+', header=None, skiprows=[0,1])

    for pdbnum in range(start+1, end+1):
        if check_file(f"{prefix}{pdbnum}{suffix}"):
            file = pd.read_csv(f"{prefix}{pdbnum}{suffix}", sep=r'\s+', header=None, skiprows=[0,1])
            file_concat = pd.concat([file_concat, file], axis=0, ignore_index=True)

    return file_concat

def concat_files_prot(prefix, suffix, start, end):
    """
    Concatenate .txt files from PackMem_prot.py for the protein

    --------------------
    INPUT
    prefix: str
        The name of the lipid usually but could just be the common name of the outputfiles
    suffix: str
        The end name of the .txt files
    start: int
        The minimum frame number in the filenames
    end: int
        The maximum frame number in the filenames

    --------------------
    OUTPUT
    pandas DataFrame
        Contains the informations in the .txt files
    """
    if suffix[:4] == "_Up_":
        other_suffix = f"_Lo_{suffix[4:]}"
    else:
        other_suffix = f"_Up_{suffix[4:]}"
    
    filename, start = start_file(prefix, suffix, start, end)
    file_concat = pd.read_csv(filename, header=None)
    
    for pdbnum in range(start+1, end+1):
        if Path(f"{prefix}{pdbnum}{suffix}").is_file():
            if check_file(f"{prefix}{pdbnum}{suffix}"):
                file = pd.read_csv(f"{prefix}{pdbnum}{suffix}", header=None)
                file_concat = pd.concat([file_concat, file], axis=0, ignore_index=True)
        elif Path(f"{prefix}{pdbnum}{other_suffix}").is_file():
            if check_file(f"{prefix}{pdbnum}{other_suffix}"):
                file = pd.read_csv(f"{prefix}{pdbnum}{other_suffix}", header=None)
                file_concat = pd.concat([file_concat, file], axis=0, ignore_index=True)

    return file_concat

def launch(output_dir, prefix, start, end, prot):
    """
    Launch the concatenation of the PackMem2 produced files
    """
    Total_Up_Deep = concat_files(f"{output_dir}/{prefix}", "_Up_Deep_result.txt", start, end)
    Total_Lo_Deep = concat_files(f"{output_dir}/{prefix}", "_Lo_Deep_result.txt", start, end)
    Total_Up_All = concat_files(f"{output_dir}/{prefix}", "_Up_All_result.txt", start, end)
    Total_Lo_All = concat_files(f"{output_dir}/{prefix}", "_Lo_All_result.txt", start, end)
    Total_Up_Shallow = concat_files(f"{output_dir}/{prefix}", "_Up_Shallow_result.txt", start, end)
    Total_Lo_Shallow = concat_files(f"{output_dir}/{prefix}", "_Lo_Shallow_result.txt", start, end)

    if prot:
        if Path(f"{output_dir}/Prot_{prefix}0_Up_All.txt").is_file():
            Total_Up_Deep_prot = concat_files_prot(f"{output_dir}/Prot_{prefix}", "_Up_Deep.txt", start, end)
            Total_Up_All_prot = concat_files_prot(f"{output_dir}/Prot_{prefix}", "_Up_All.txt", start, end)
            Total_Up_Shallow_prot = concat_files_prot(f"{output_dir}/Prot_{prefix}", "_Up_Shallow.txt", start, end)
            # Save files
            Total_Up_Deep_prot.to_csv(f"{output_dir}/Total_Up_Deep_prot.csv", header=False, index=False)
            Total_Up_All_prot.to_csv(f"{output_dir}/Total_Up_All_prot.csv", header=False, index=False)
            Total_Up_Shallow_prot.to_csv(f"{output_dir}/Total_Up_Shallow_prot.csv", header=False, index=False)
            # Remove the files
            for pdbnum in range(start, end+1):
                os.remove(f"{output_dir}/Prot_{prefix}{pdbnum}_Up_Deep.txt")
                os.remove(f"{output_dir}/Prot_{prefix}{pdbnum}_Up_All.txt")
                os.remove(f"{output_dir}/Prot_{prefix}{pdbnum}_Up_Shallow.txt")
        elif Path(f"{output_dir}/Prot_{prefix}0_Lo_All.txt").is_file():
            Total_Lo_Deep_prot = concat_files_prot(f"{output_dir}/Prot_{prefix}", "_Lo_Deep.txt", start, end)
            Total_Lo_All_prot = concat_files_prot(f"{output_dir}/Prot_{prefix}", "_Lo_All.txt", start, end)
            Total_Lo_Shallow_prot = concat_files_prot(f"{output_dir}/Prot_{prefix}", "_Lo_Shallow.txt", start, end)
            # Save files
            Total_Lo_Deep_prot.to_csv(f"{output_dir}/Total_Lo_Deep_prot.csv", header=False, index=False)
            Total_Lo_All_prot.to_csv(f"{output_dir}/Total_Lo_All_prot.csv", header=False, index=False)
            Total_Lo_Shallow_prot.to_csv(f"{output_dir}/Total_Lo_Shallow_prot.csv", header=False, index=False)
            # Remove the files
            for pdbnum in range(start, end+1):
                os.remove(f"{output_dir}/Prot_{prefix}{pdbnum}_Lo_Deep.txt")
                os.remove(f"{output_dir}/Prot_{prefix}{pdbnum}_Lo_All.txt")
                os.remove(f"{output_dir}/Prot_{prefix}{pdbnum}_Lo_Shallow.txt")
        
    # Concatenate the lealfets' results
    Total_Deep =  pd.concat([Total_Up_Deep, Total_Lo_Deep], axis=0, ignore_index=True)
    Total_All =  pd.concat([Total_Up_All, Total_Lo_All], axis=0, ignore_index=True)
    Total_Shallow =  pd.concat([Total_Up_Shallow, Total_Lo_Shallow], axis=0, ignore_index=True)

    # Save all files to .csv
    Total_Up_Deep.to_csv(f"{output_dir}/Total_Up_Deep.csv", header=False, index=False)
    Total_Lo_Deep.to_csv(f"{output_dir}/Total_Lo_Deep.csv", header=False, index=False)
    Total_Up_All.to_csv(f"{output_dir}/Total_Up_All.csv", header=False, index=False)
    Total_Lo_All.to_csv(f"{output_dir}/Total_Lo_All.csv", header=False, index=False)
    Total_Up_Shallow.to_csv(f"{output_dir}/Total_Up_Shallow.csv", header=False, index=False)
    Total_Lo_Shallow.to_csv(f"{output_dir}/Total_Lo_Shallow.csv", header=False, index=False)

    Total_Deep.to_csv(f"{output_dir}/Total_Deep.csv", header=False, index=False)
    Total_All.to_csv(f"{output_dir}/Total_All.csv", header=False, index=False)
    Total_Shallow.to_csv(f"{output_dir}/Total_Shallow.csv", header=False, index=False)

    # Remove the files
    for pdbnum in range(start, end+1):
        os.remove(f"{output_dir}/{prefix}{pdbnum}_Up_Deep_result.txt")
        os.remove(f"{output_dir}/{prefix}{pdbnum}_Lo_Deep_result.txt")
        os.remove(f"{output_dir}/{prefix}{pdbnum}_Up_All_result.txt")
        os.remove(f"{output_dir}/{prefix}{pdbnum}_Lo_All_result.txt")
        os.remove(f"{output_dir}/{prefix}{pdbnum}_Up_Shallow_result.txt")
        os.remove(f"{output_dir}/{prefix}{pdbnum}_Lo_Shallow_result.txt")

def main():
    args = get_args()

    launch(args.output_dir, args.prefix, args.start, args.end, args.protein)


if __name__=="__main__":
    main()
