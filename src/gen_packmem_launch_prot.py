import string
import argparse

# Getting the arguments
parser = argparse.ArgumentParser()
parser.add_argument('-c', action = 'store', dest = 'cores',
                    help = 'The number of cores you want the script to work on',
                    default = 1, type = int)
parser.add_argument('-f', action = 'store', dest = 'frames',
                    help = 'The number of frames in your simulation',
                    type = int)
parser.add_argument('-p', action = 'store', dest = 'prefix',
                    help = 'The name you gave as prefix for the frames')
args = parser.parse_args()

################################################################################

if __name__=="__main__":
    # Prepare how much frames are going to be analysed per core
    step = int(args.frames / args.cores)
    intervals = list(range(0, args.frames+1, step))

    # Loop to create the prompt to launch PackMem on several cores
    for i in range(len(intervals)-1):
        start = intervals[i]
        stop =  intervals[i+1]-1
        if i==(len(intervals)-2):
            stop = args.frames
        print(f"nice -19 bash ../../New_PackMem/src/ScriptPackMem_prot.sh {args.prefix} {start} {stop} >& OUT_packmem{i} &")


