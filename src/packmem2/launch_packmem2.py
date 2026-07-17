from multiprocessing import Process

from packmem2.core import param as p
import packmem2.packmem2
import packmem2.concatenate
import packmem2.analysis

def launch(cores, topo, traj, lipid, start, end, paramFile, radiiFile, indexFile, output_dir, outputname, dist_suppl_Z, protein, pdbout, limx, limy, precision):
    """
    Launch the different python scripts for PackMem2
    """
    ######## PackMem2 #######
    # Prepare how much frames are going to be analysed per core
    step = int(end / cores)
    # Create a list of the minima / maxima frame of each interval
    intervals = list(range(start, end+1, step))
    processes = []

    # Launch PackMem2 on multiple processes
    for i in range(len(intervals)-1):
        # Select starting and ending frame (limits of the interval)
        start_i = intervals[i]
        end_i =  intervals[i+1]-1
        # If this is the last interval
        if i == (len(intervals)-2):
            end_i = end
        p = Process(target = packmem2.packmem2.launch,
                    args = (topo, traj, lipid, start_i, end_i, paramFile,\
                            radiiFile, indexFile, output_dir, outputname,\
                            dist_suppl_Z, protein, pdbout)
        )
        p.start()
        processes.append(p)
    
    # Wait for the processes to finish
    for p in processes:
        p.join()

    ######## Concatenate #######
    packmem2.concatenate.launch(output_dir, outputname, start, end, protein)

    ######## Analysis #######
    # Then launch packing defect analysis
    packmem2.analysis.launch(output_dir, outputname, protein,\
                            limx, limy, precision)

def main():
    args = p.get_args_launch_packmem2()

    launch(args.cores, args.topo, args.traj, args.lipid, args.start, args.end,\
           args.paramFile, args.radiiFile, args.indexFile,\
           args.output_dir, args.outputname, args.dist_suppl_Z, args.protein,\
           args.pdbout, args.limx, args.limy, args.precision)


if __name__=="__main__":
    main()