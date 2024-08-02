import math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-mol', action = 'store', dest = 'fname1',
                help = 'the .ipt file for the molecule')
parser.add_argument('-ff', action = 'store', dest = 'fname2',
                help = 'forcefield.itp')
args = parser.parse_args()

if __name__=="__main__":
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
                print(f"{line.split()[3]}  {line.split()[4]:4s} {radii:.2f} n")
            if "; nr	type" in line:
                flag = True

