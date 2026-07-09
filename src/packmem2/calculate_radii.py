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
                        required = True,
                        help = 'the .ipt file for the molecule')
    parser.add_argument('-ff', action = 'store', dest = 'ff',
                        required = True,
                        help = 'forcefield.itp')
    parser.add_argument('-o', action = 'store', dest = 'output',
                        default = 'vdw_out',
                        help = 'The name of the output .txt file')
    parser.add_argument('-martini', action = 'store_true', dest = 'martini',
                        help = 'If you want to compute the radii with the MARTINI forcefield')
    parser.add_argument('-martini3', action = 'store_true', dest = 'martini3',
                        help = 'If you specificaly have the MARTINI3 force field')
    args = parser.parse_args()
    return args

def get_aliphatic(res_name, res_name_prev, atom_name, flag):
    """
    Check if the atom or bead is aliphatic or not

    --------------------
    INPUT
    res_name: str
        The name of the residue
    res_name_prev: str
        The name of the residue on the line before
    atom_name: str
        The name of the atom or bead
    flag: bool
        If the central atom of glycerol or equivalent already passed
    
    --------------------
    OUTPUT
    str, bool
        The aliphatic corresponding letter ['n' for neutral or 'a' for apolar]
        The new value of the flag
    """
    # Check if we changed residue
    if res_name != res_name_prev:
        flag = False

    # Check if we passed the central atom of glycerol or equivalent already passed
    if ((args.martini or args.martini3) and (atom_name == 'GL2' or atom_name == 'AM1' or atom_name == 'R1')) or ((not args.martini and not args.martini3) and atom_name == 'C2'):
        flag = True
        return 'n', flag
    elif flag:
        return 'a', flag
    else:
        return 'n', flag

def main():
    """
    Calcultate the radii of given molecules
    """
    args = get_args()
    
    sigma = {}
    flag_ff = False
    flag_mol = False
    flag_aliph = False
    res_name_prev = ""

    # Get the value of sigma from the forcefield .itp
    with open(args.ff, "r") as file_ff:
        for line in file_ff:
            if line.startswith(';'):
                continue
            if line.startswith('['):
                flag_ff = False 
            if flag_ff and len(line.strip()) != 0:
                if (args.martini or args.martini3) and line.split()[0] == line.split()[1]:
                    if args.martini3:
                        sigma[line.split()[0]] = float(line.split()[3])
                    else:
                        # sigma = (C12/C6)^(1/6)
                        C12 = float(line.split()[4])
                        C6 = float(line.split()[3])
                        if C6 >= 0.00001 or C6 <= -0.00001:
                            sigma[line.split()[0]] = (C12 / C6)**(1/6)
                        else:
                            sigma[line.split()[0]] = 0.0
                elif not args.martini and not args.martini3:
                    sigma[line.split()[0]] = float(line.split()[5])  
            if ((args.martini or args.martini3)  and "[ nonbond_params ]" in line) or ((not args.martini and not args.martini3) and "[ atomtypes ]" in line):
                flag_ff = True

    # Calculate the radii and prints the radii of each atom type
    with open(args.mol, "r") as file_mol, open(f"{args.output}.txt", "w") as file_out:
        for line in file_mol:
            if line.startswith(';'):
                    continue
            if line.startswith('['):
                flag_mol = False 
            if flag_mol and len(line.strip()) != 0:
                radii = ((math.pow(2,(1/6))*sigma[line.split()[1]])/2)*10
                res_name = line.split()[3]
                atom_name = line.split()[4]
                aliph, flag_aliph = get_aliphatic(res_name, res_name_prev, atom_name, flag_aliph)
                res_name_prev = res_name
                file_out.write(f"{res_name:<6s}{atom_name:4s} {radii:.3f} {aliph}\n")
            if "[ atoms ]" in line or "[atoms]" in line:
                flag_mol = True

if __name__=="__main__":
    main()