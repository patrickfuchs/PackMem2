from pathlib import Path

import packmem2.launch_packmem2 as launch_packmem2

def test_launch():
    cores = 1
    topo = "tests/data/md.gro"
    traj = "tests/data/md_40ns.xtc"
    lipid = "DMPC"
    start = 0
    end = 400
    paramFile = "data/param_Charmm.txt"
    radiiFile = "data/vdw_radii_Charmm.txt"
    indexFile = None
    outputname = "DMPC"
    dist_suppl_Z = 1.0
    protein = False
    pdbout = False
    limx = 15
    limy = 1e-4
    precision = 2

    launch_packmem2.launch(cores, topo, traj, lipid, start, end, paramFile,\
                           radiiFile, indexFile, outputname, dist_suppl_Z,\
                           protein, pdbout, limx, limy, precision)
    
    # Check Total files
    expected_output_Deep = Path("Total_Deep.csv")
    expected_output_Shallow = Path("Total_Shallow.csv")
    expected_output_All = Path("Total_All.csv")
    with open("tests/data/Total_Deep.csv", 'r') as f_in:
        expected_content_Deep = f_in.read()

    assert expected_output_Deep.exists()
    assert expected_output_Shallow.exists()
    assert expected_output_All.exists ()
    assert expected_output_Deep.read_text() == expected_content_Deep

    # Check pdf
    expected_final_output = Path("DMPC.pdf")
    assert expected_final_output.exists()
