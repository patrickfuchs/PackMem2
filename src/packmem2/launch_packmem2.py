import sys

from packmem2.core import param as p
import packmem2.packmem2
import packmem2.concatenate
import packmem2.analysis

def main():
    """
    Launch the different python scripts for PackMem2
    """
    args = p.get_args_launch_packmem2()
    
    ######## PackMem2 #######
    packmem2.packmem2.launch(args.topo, args.traj, args.lipid, args.start,\
                             args.end, args.paramFile, args.radiiFile, args.indexFile,\
                             args.outputname, args.dist_suppl_Z, args.protein, args.pdbout)

    ######## Concatenate #######
    packmem2.concatenate.launch(args.outputname, args.start, args.end, args.protein)

    ######## Analysis #######
    # Then launch packing defect analysis
    packmem2.analysis.launch(args.outputname, args.protein, args.limx, args.limy, args.precision)
    

if __name__=="__main__":
    main()