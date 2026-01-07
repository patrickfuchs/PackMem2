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
    parser.add_argument('-mol', action = 'store', dest = 'mol',
                    help = 'the .ipt file for the molecule')
    parser.add_argument('-ff', action = 'store', dest = 'ff',
                    help = 'forcefield.itp')
    parser.add_argument('-martini', action = 'store_true', dest = 'martini',
                    help = 'If you want to compute the radii with the MARTINI forcefield')
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
    with open(args.ff, "r") as file_ff:
        for line in file_ff:
            if args.martini:
                if "; cross terms" in line:
                    flag = False 
                if flag and len(line.strip()) != 0:
                    # sigma = (C12/C6)^(1/6)
                    C12 = float(line.split()[4])
                    C6 = float(line.split()[3])
                    if C6 >= 0.00001 or C6 <= -0.00001:
                        sigma[line.split()[0]] = (C12 / C6)**(1/6)
                    else:
                        print(C6)
                        sigma[line.split()[0]] = 0.0
                if "; self terms" in line:
                    flag = True
            else:
                if "[ nonbond_params ]" in line:
                    flag = False 
                if flag and len(line.strip()) != 0:
                    sigma[line.split()[0]] = float(line.split()[5])       
                if "; name	at.num" in line:
                    flag = True

    # Calculate the radii and prints the radii of each atom type
    with open(args.mol, "r") as file_mol:
        for line in file_mol:
            if "[ bonds ]" in line:
                flag = False 
            if flag and len(line.strip()) != 0:
                radii = ((math.pow(2,(1/6))*sigma[line.split()[1]])/2)*10
                aliph = aliphatic(line.split()[4])
                print(f"{line.split()[3]}  {line.split()[4]:4s} {radii:.2f} {aliph}")
            if (args.martini and "[ atoms ]" in line) or (not args.martini and "; nr	type" in line):
                flag = True

