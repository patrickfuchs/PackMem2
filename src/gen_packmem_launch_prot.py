import argparse
import sys

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
    parser.add_argument('-prot', action='store_true', dest='prot', 
                        help='Put if you want to see the packing defects close/far of the protein')
    args = parser.parse_args()

    return args

################################################################################

if __name__=="__main__":
    args = get_arguments()
    path = sys.argv[0]
    path = path.replace("src/gen_packmem_launch_prot.py", "")

    # Prepare how much frames are going to be analysed per core
    step = int(args.frames / args.cores)
    intervals = list(range(0, args.frames+1, step))

    print("\nPlease, copy paste these lines to launch PackMem:\n")

    #If too few frames
    if args.frames < (args.cores*2):
        if args.prot:
            print(f"nice -19 python {path}src/PackMem_prot.py -f {args.traj} -s {args.topo} -l {args.lipid_name} -r {path}data/vdw_radii_Charmm.txt -p {path}data/param_Charmm.txt -b 0 -e {args.frames} -prot -o {args.lipid_name} >& OUT_packmem0 &")
        else:
            print(f"nice -19 python {path}src/PackMem_prot.py -f {args.traj} -s {args.topo} -l {args.lipid_name} -r {path}data/vdw_radii_Charmm.txt -p {path}data/param_Charmm.txt -b 0 -e {args.frames} -o {args.lipid_name} >& OUT_packmem0 &")


    # Loop to create the prompt to launch PackMem on several cores
    for i in range(len(intervals)-1):
        start = intervals[i]
        stop =  intervals[i+1]-1
        if i==(len(intervals)-2):
            stop = args.frames
        if args.prot:
            print(f"nice -19 python {path}src/PackMem_prot.py -f {args.traj} -s {args.topo} -l {args.lipid_name} -r {path}data/vdw_radii_Charmm.txt -p {path}data/param_Charmm.txt -b {start} -e {stop} -prot -o {args.lipid_name} >& OUT_packmem0 &")
        else:
            print(f"nice -19 python {path}src/PackMem_prot.py -f {args.traj} -s {args.topo} -l {args.lipid_name} -r {path}data/vdw_radii_Charmm.txt -p {path}data/param_Charmm.txt -b {start} -e {stop} -o {args.lipid_name} >& OUT_packmem0 &")
