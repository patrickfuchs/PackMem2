from pathlib import Path

import packmem2.launch_packmem2 as launch_packmem2

def test_launch():
    ################# DMPC ###############
    cores = 1
    topo = "tests/data/end_to_end_DMPC/md.gro"
    traj = "tests/data/end_to_end_DMPC/md_10ns.xtc"
    lipid = "DMPC"
    start = 0
    end = 100
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
    with open("tests/data/end_to_end_DMPC/Total_Deep.csv", 'r') as f_in:
        expected_content_Deep = f_in.read()

    assert expected_output_Deep.exists()
    assert expected_output_Shallow.exists()
    assert expected_output_All.exists ()
    assert expected_output_Deep.read_text() == expected_content_Deep

    # Check pdf
    expected_final_output = Path("DMPC.pdf")
    assert expected_final_output.exists()


    ################# DMPC + protein ###############
    cores = 1
    topo = "tests/data/end_to_end_DMPC_prot/md.gro"
    traj = "tests/data/end_to_end_DMPC_prot/md_10ns.xtc"
    lipid = "DMPC"
    start = 0
    end = 100
    paramFile = "data/param_Charmm.txt"
    radiiFile = "data/vdw_radii_Charmm.txt"
    indexFile = None
    outputname = "DMPC"
    dist_suppl_Z = 1.0
    protein = True
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
    with open("tests/data/end_to_end_DMPC_prot/Total_Deep.csv", 'r') as f_in:
        expected_content_Deep = f_in.read()

    assert expected_output_Deep.exists()
    assert expected_output_Shallow.exists()
    assert expected_output_All.exists ()
    assert expected_output_Deep.read_text() == expected_content_Deep

    # Check Total files for protein
    expected_output_Deep = Path("Total_Up_Deep_prot.csv")
    expected_output_Shallow = Path("Total_Up_Shallow_prot.csv")
    expected_output_All = Path("Total_Up_All_prot.csv")
    with open("tests/data/end_to_end_DMPC_prot/Total_Up_Deep_prot.csv", 'r') as f_in:
        expected_content_Deep = f_in.read()

    assert expected_output_Deep.exists()
    assert expected_output_Shallow.exists()
    assert expected_output_All.exists ()
    assert expected_output_Deep.read_text() == expected_content_Deep

    # Check pdf
    expected_final_output = Path("DMPC.pdf")
    assert expected_final_output.exists()


################# DLPC (Martini) ###############
    cores = 1
    topo = "tests/data/end_to_end_DLPC/md.gro"
    traj = "tests/data/end_to_end_DLPC/md_10ns.xtc"
    lipid = "DLPC"
    start = 0
    end = 100
    paramFile = "data/param_Martini.txt"
    radiiFile = "data/vdw_radii_Martini_old.txt"
    indexFile = None
    outputname = "DLPC"
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
    with open("tests/data/end_to_end_DLPC/Total_Deep.csv", 'r') as f_in:
        expected_content_Deep = f_in.read()

    assert expected_output_Deep.exists()
    assert expected_output_Shallow.exists()
    assert expected_output_All.exists ()
    assert expected_output_Deep.read_text() == expected_content_Deep

    # Check pdf
    expected_final_output = Path("DLPC.pdf")
    assert expected_final_output.exists()