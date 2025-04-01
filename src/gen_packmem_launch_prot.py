import string
import argparse

# Getting the arguments
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
parser.add_argument('-n', action='store', dest='prefix', 
                    help='The type of lipid')
parser.add_argument('-prot', action='store_true', dest='prot', 
                    help='Put if you want to see the packing defects close/far of the protein')
args = parser.parse_args()

################################################################################

if __name__=="__main__":
    # Prepare how much frames are going to be analysed per core
    step = int(args.frames / args.cores)
    intervals = list(range(0, args.frames+1, step))

    #If too few frames
    if args.frames < (args.cores*2):
        if args.prot:
            print(f"nice -19 bash /home/maya/Documents/tools/New_PackMem/src/ScriptPackMem_prot.sh -f {args.traj} -s {args.topo} -n {args.prefix} -b 0 -e {args.frame} -p >& OUT_packmem0 &")
        else:
            print(f"nice -19 bash /home/maya/Documents/tools/New_PackMem/src/ScriptPackMem_prot.sh -f {args.traj} -s {args.topo} -n {args.prefix} -b 0 -e {args.frame} >& OUT_packmem0 &")


    # Loop to create the prompt to launch PackMem on several cores
    for i in range(len(intervals)-1):
        start = intervals[i]
        stop =  intervals[i+1]-1
        if i==(len(intervals)-2):
            stop = args.frames
        if args.prot:
            print(f"nice -19 bash /home/maya/Documents/tools/New_PackMem/src/ScriptPackMem_prot.sh -f {args.traj} -s {args.topo} -n {args.prefix} -b {start} -e {stop} -p >& OUT_packmem{i} &")
        else:
            print(f"nice -19 bash /home/maya/Documents/tools/New_PackMem/src/ScriptPackMem_prot.sh -f {args.traj} -s {args.topo} -n {args.prefix} -b {start} -e {stop} >& OUT_packmem{i} &")
