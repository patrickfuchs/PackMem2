"""Calculate the radii for PackMem2 van der Waals data"""

import math
import argparse
import importlib.metadata


def get_args() -> argparse.Namespace:
    """
    Get the arguments for the script and check that the inputfiles are valid.

    --------------------
    OUTPUT
    parser.parse_args
        Contains all the arguments for the script
    """
    packmem2_version = importlib.metadata.version('packmem2')
    # Getting the arguments
    parser = argparse.ArgumentParser(description=f"Arguments for calculate_radii of PackMem2 version {packmem2_version}")
    parser.add_argument(
        "-m",
        action="store",
        dest="mol_file",
        required=True,
        help="the .ipt file for the molecule",
    )
    parser.add_argument(
        "-f", action="store", dest="ff_file", required=True, help="forcefield.itp"
    )
    parser.add_argument(
        "-ff",
        action="store",
        dest="forcefield",
        required=True,
        help="The force field of the files given. Options: charmm, martini, martini3",
    )
    parser.add_argument(
        "-o",
        action="store",
        dest="output",
        default="vdw_out",
        help="The name of the output .txt file",
    )
    parser.add_argument(
        "--version",
        action='version',
        version=packmem2_version
    )
    
    args = parser.parse_args()
    return args


def get_aliphatic(
    res_name: str,
    res_name_prev: str,
    atom_name: str,
    flag: bool,
    forcefield: str
) -> tuple[str, bool]:
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
    forcefield :str
        The forcefield of the .itp files given

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
    if (
        ("martini" in forcefield)
        and (atom_name == "GL2" or atom_name == "AM1" or atom_name == "R1")
    ) or ((forcefield == "charmm") and atom_name == "C2"):
        flag = True
        return "p", flag
    elif flag:
        return "a", flag
    else:
        return "p", flag


def launch(ff_file: str, mol_file: str, forcefield: str, output: str) -> None:
    """
    Calcultate the radii of given molecules

    --------------------
    INPUT
    ff_file: str
        The name of the forcef field .itp file
    mol_file: str
        The name of the molecule .itp file
    forcefield: str
        The forcefield of the .itp files
        Can be: charmm, martini or martini3
    output: str
        The name of the output file
    """
    sigma = {}
    flag_ff = False
    flag_mol = False
    flag_aliph = False
    res_name_prev = ""

    # Get the value of sigma from the forcefield .itp
    with open(ff_file, "r") as file_ff:
        for line in file_ff:
            if line.startswith(";"):
                continue
            if line.startswith("["):
                flag_ff = False
            if flag_ff and len(line.strip()) != 0:
                if ("martini" in forcefield) and line.split()[0] == line.split()[
                    1
                ]:
                    if forcefield == "martini3":
                        sigma[line.split()[0]] = float(line.split()[3])
                    else:
                        # sigma = (C12/C6)^(1/6)
                        C12 = float(line.split()[4])
                        C6 = float(line.split()[3])
                        if C6 >= 0.00001 or C6 <= -0.00001:
                            sigma[line.split()[0]] = (C12 / C6) ** (1 / 6)
                        else:
                            sigma[line.split()[0]] = 0.0
                elif forcefield == "charmm":
                    sigma[line.split()[0]] = float(line.split()[5])
            if (("martini" in forcefield) and "[ nonbond_params ]" in line) or (
                (forcefield == "charmm") and "[ atomtypes ]" in line
            ):
                flag_ff = True

    # Calculate the radii and prints the radii of each atom type
    with open(mol_file, "r") as file_mol, open(f"{output}.txt", "w") as file_out:
        for line in file_mol:
            if line.startswith(";"):
                continue
            if line.startswith("["):
                flag_mol = False
            if flag_mol and len(line.strip()) != 0:
                radii = ((math.pow(2, (1 / 6)) * sigma[line.split()[1]]) / 2) * 10
                res_name = line.split()[3]
                atom_name = line.split()[4]
                aliph, flag_aliph = get_aliphatic(
                    res_name, res_name_prev, atom_name, flag_aliph, forcefield
                )
                res_name_prev = res_name
                if "martini" in forcefield:
                    file_out.write(f"{res_name:<6s}{atom_name:4s} {radii:.3f} {aliph}\n")
                else:
                    file_out.write(f"{res_name:<6s}{atom_name:4s} {radii:.2f} {aliph}\n")
            if "[ atoms ]" in line or "[atoms]" in line:
                flag_mol = True


def main() -> None:
   
    args = get_args()

    launch(args.ff_file, args.mol_file, args.forcefield, args.output)


if __name__ == "__main__":
    main()
