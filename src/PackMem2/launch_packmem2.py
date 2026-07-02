import argparse
import sys
import subprocess

def get_arguments():
    """
    Get the arguments for the script and check that the inputfiles are valid.

    --------------------
    OUTPUT
    parser.parse_args
        Contains all the arguments for the script
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', action = 'store', dest = 'traj',
                        help = 'The trajectory file (.xtc)')
    parser.add_argument('-s', action = 'store', dest = 'topo',
                        help = 'The topology file (.gro)')
    parser.add_argument('-c', action = 'store', dest = 'cores',
                        help = 'The number of cores you want the script to work on',
                        default = 1, type = int)
    parser.add_argument('-fm', action = 'store', dest = 'frames',
                        help = 'The number of frames in your simulation',
                        type = int)
    parser.add_argument('-l', action='store', dest='lipid_name', 
                        help='The name of lipid(s). If multiple, seperate them with _')
    parser.add_argument('-p', action = 'store', dest = 'param',
                        help = 'The parameter file to be used. \nThe options are: param_Charmm.txt (default), param_Martini.txt',
                        default = "param_Charmm.txt")
    parser.add_argument('-r', action = 'store', dest = 'radii',
                        help = 'The radii file to be used.\nThe options are: vdw_radii_Charmm.txt (default), vdw_radii_Martini2.txt, vdw_radii_Martini2P.txt and vdw_radii_Martini3.txt',
                        default = "vdw_radii_Charmm.txt")
    parser.add_argument('-prot', action='store_true', dest='prot', 
                        help='Put if you want to see the packing defects close/far of the protein')
    args = parser.parse_args()

    return args

def main():
    """
    Launch the different python scripts for PackMem2
    """
    args = get_arguments()
    path = sys.argv[0]
    path = path.replace("launch_packmem2.py", "")

    # Prepare how much frames are going to be analysed per core
    step = int(args.frames / args.cores)
    intervals = list(range(0, args.frames+1, step))

    cmd = [
            "nice", "-19",
            sys.executable,  # ou "python3" si tu préfères
            f"{path}packmem2.py",
            "-f", args.traj,
            "-s", args.topo,
            "-l", args.lipid_name,
            "-r", f"{path}data/{args.radii}",
            "-p", f"{path}data/{args.param}",
            "-o", args.lipid_name
        ]
    if args.prot:
        cmd.append("-prot")

    processes = []
    out_file = open(f"OUT_packmem", "a")
    #If too few frames
    if args.frames < (args.cores):
        cmd += ["-b", "0", "-e", str(args.frames)]
        process = subprocess.Popen(cmd, stdout=out_file, stderr=subprocess.STDOUT)
        processes.append(process)
    else:
        # Loop to launch PackMem on several cores
        for i in range(len(intervals)-1):
            start = intervals[i]
            stop =  intervals[i+1]-1
            if i==(len(intervals)-2):
                stop = args.frames
            # Add commands
            cmd_i = cmd + ["-b", str(start), "-e", str(stop)]
            # Lancement du processus en arrière-plan
            process = subprocess.Popen(cmd_i, stdout=out_file, stderr=subprocess.STDOUT)
            processes.append(process)

    # Wait for the processes to finish
    for p in processes:
        p.wait()
    out_file.close()


    # Once everything is finished, launch concatenation
    cmd = [
            sys.executable,  # ou "python3" si tu préfères
            f"{path}concatenate.py",
            "-l", args.lipid_name,
            "-b", "0",
            "-e", str(args.frames)
        ]
    if args.prot:
        cmd.append("-prot")

    process = subprocess.Popen(cmd)
    process.wait()

    # Then launch packing defect analysis
    cmd = [
            sys.executable,  # ou "python3" si tu préfères
            f"{path}analysis.py"
        ]
    if args.prot:
        cmd.append("-prot")
    
    process = subprocess.Popen(cmd)

if __name__=="__main__":
    main()