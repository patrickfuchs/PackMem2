import sys

from packmem2.core import param as p
import packmem2.packmem2
import packmem2.concatenate
import packmem2.analysis

def launch(cores, topo, traj, lipid, start, end, paramFile, radiiFile, indexFile, outputname, dist_suppl_Z, protein, pdbout, limx, limy, precision):
    """
    Launch the different python scripts for PackMem2
    """    
    ######## PackMem2 #######
    packmem2.packmem2.launch(topo, traj, lipid, start, end, paramFile,\
                             radiiFile, indexFile, outputname, dist_suppl_Z,\
                             protein, pdbout)

    ######## Concatenate #######
    packmem2.concatenate.launch(outputname, start, end, protein)

    ######## Analysis #######
    # Then launch packing defect analysis
    packmem2.analysis.launch(outputname, protein, limx, limy, precision)

def main():
    args = p.get_args_launch_packmem2()

    launch(args.cores, args.topo, args.traj, args.lipid, args.start, args.end,\
           args.paramFile, args.radiiFile, args.indexFile,\
           args.outputname, args.dist_suppl_Z, args.protein, args.pdbout,\
           args.limx, args.limy, args.precision)


if __name__=="__main__":
    main()