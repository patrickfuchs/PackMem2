import math
import argparse

def get_args():
    """
    Get the arguments for the script and check that the inputfiles are valid.

    --------------------
    OUTPUT
    parser.parse_args
        Contains all the arguments for the script
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-mol', action = 'store', dest = 'fname1',
                    help = 'the .ipt file for the molecule')
    parser.add_argument('-ff', action = 'store', dest = 'fname2',
                    help = 'forcefield.itp')
    args = parser.parse_args()
    return args

def aliphatic(atom_name):
    """
    Checks if the atom if aliphatic or not

    -------------------
    INPUT
    atom_name: str
        The name of the atom
    
    -------------------
    OUTPUT
    str
        The nature of the atom (a / n)
    """
    if len(atom_name) < 3:
        return 'n'
    else:
        if atom_name[0] != 'H' and len(atom_name) > 3:
            if atom_name[1] == "1" and (atom_name[3] == 'A' or atom_name[3] == 'B' or atom_name[3] == 'C'):
                return 'n'
            else:
                return 'a'
        elif atom_name[0] == 'H' and len(atom_name) == 3:
            if atom_name[1] == 'N' or atom_name[1] == 'O':
                return 'n'
            else:
                return 'a'
        elif atom_name[0] == 'H' and len(atom_name) > 3:
            if atom_name[1] == "1" and (atom_name[3] == 'A' or atom_name[3] == 'B' or atom_name[3] == 'C'):
                return 'n'
        elif atom_name[0] != 'H' and len(atom_name) == 3:
            if atom_name[1] == "1" and int(atom_name[2]) > 0:
                return 'n'
            elif atom_name[1] == "C":
                return 'n'
            elif int(atom_name[1]) > 1 and atom_name[2] == "1":
                return 'n'
            elif atom_name == 'O22' or atom_name == 'O32':
                return 'n'
        return 'a'



if __name__=="__main__":
    args = get_args()
    
    sigma = {}
    flag = False

    # Get the value of sigma from the forcefield .itp
    with open(args.fname2, "r") as file_ff:
        for line in file_ff:
            if "[ nonbond_params ]" in line:
                flag = False 
            if flag and len(line.strip()) != 0:
                sigma[line.split()[0]] = float(line.split()[5])       
            if "; name	at.num" in line:
                flag = True

    # Calculate the radii and prints the radii of each atom type
    with open(args.fname1, "r") as file_mol:
        for line in file_mol:
            if "[ bonds ]" in line:
                flag = False 
            if flag and len(line.strip()) != 0:
                radii = ((math.pow(2,(1/6))*sigma[line.split()[1]])/2)*10
                aliph = aliphatic(line.split()[4])
                print(f"{line.split()[3]}  {line.split()[4]:4s} {radii:.2f} {aliph}")
            if "; nr	type" in line:
                flag = True

