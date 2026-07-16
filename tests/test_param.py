import pytest
import sys
import numpy as np

import packmem2.core.param as p

def test_file_present():
    # Test if file present
    filename_OK = "src/packmem2/packmem2.py"
    p.file_present(str(filename_OK))

    # Test if file absent
    filename_NO = "piou.txt"
    with pytest.raises(FileNotFoundError):
        p.file_present(str(filename_NO))

def test_get_args_packmem2(monkeypatch, tmp_path):
    # Create tmp files
    traj = tmp_path / "traj.xtc"
    topo = tmp_path / "topo.gro"
    param = tmp_path / "param_Charmm.txt"
    radii = tmp_path / "vdw_radii_Charmm.txt"
    for f in [traj, topo, param, radii]:
        f.write_text("dummy content")

    ## TEST OK ##
    # Simulate the corect arguments
    monkeypatch.setattr(
        sys, "argv",
        [
            "prog",
            "-f", str(traj),
            "-s", str(topo),
            "-l", "DPPC",
            "-e", "10",
            "-p", str(param),
            "-r", str(radii),
            "-d", "1.0",
            "-prot"
        ]
    )
    args = p.get_args_packmem2()
    assert args.traj == str(traj)
    assert args.topo == str(topo)
    assert args.lipid == "DPPC"
    assert args.start == 0
    assert args.end == 10
    assert args.outputname == "output"
    assert args.dist_suppl_Z == 1.0
    assert args.pdbout is False
    assert args.protein is True

    ## TEST NO ##
    monkeypatch.setattr(
        sys, "argv",
        [
            "prog",
            "-f", str(traj),
            "-s", str(topo),
            "-l", "DPPC",
            "-e", "10",
            "-p", str(param),
            "-r", str(radii),
            "-d", "-1.0",
        ]
    )
    with pytest.raises(Exception) as exc_info:
        p.get_args_packmem2()
    assert "must be > 0.0" in str(exc_info.value)

def test_get_args_launch_packmem2(monkeypatch, tmp_path):
    # Create tmp files
    traj = tmp_path / "traj.xtc"
    topo = tmp_path / "topo.gro"
    param = tmp_path / "param_Charmm.txt"
    radii = tmp_path / "vdw_radii_Charmm.txt"
    for f in [traj, topo, param, radii]:
        f.write_text("dummy content")

    ## TEST OK ##
    # Simulate the corect arguments
    monkeypatch.setattr(
        sys, "argv",
        [
            "prog",
            "-f", str(traj),
            "-s", str(topo),
            "-l", "DPPC",
            "-e", "10",
            "-p", str(param),
            "-r", str(radii),
            "-d", "1.0",
            "-prot"
        ]
    )
    args = p.get_args_launch_packmem2()
    assert args.cores == 1
    assert args.precision == 2
    assert args.limx == 15
    assert args.limy == 1e-4

def test_dict_2columns():
    list_str = ["DOPC C2", "DMPG C2", "POPC C2"]
    tested_ouput = p.dict_2columns(list_str)
    wanted_output = {"DOPC": "C2", "DMPG": "C2", "POPC": "C2"}
    assert tested_ouput == wanted_output

def test_set_params(tmp_path):
    filename = tmp_path / "param.txt"
    for f in [filename]:
        f.write_text('DOP C2\nDOE C2\n')
    tested_output = p.set_params(filename)
    wanted_output = {'DOP': 'C2', 'DOE': 'C2'}
    assert tested_output == wanted_output

def test_dict_4columns():
    list_str = ['DOP N 1.89 n', 'DOE C 1.45 a']
    tested_ouput_1 = p.dict_4columns(list_str, 2)
    wanted_output_1 = {"DOP N": 1.89, "DOE C": 1.45}
    assert tested_ouput_1 == wanted_output_1
    tested_ouput_2 = p.dict_4columns(list_str, 3)
    wanted_output_2 = {"DOP N": 'n', "DOE C": 'a'}
    assert tested_ouput_2 == wanted_output_2

def test_set_rad_ali(tmp_path):
    filename = tmp_path / "file.txt"
    for f in [filename]:
        f.write_text('DOP N 1.89 n\nDOE C 1.45 a\n')
    output_rad, output_ali = p.set_rad_ali(filename)
    wanted_rad = {'DOP N': 1.89, 'DOE C': 1.45}
    wanted_ali = {'DOP N': 'n', 'DOE C': 'a'}
    assert output_rad == wanted_rad
    assert output_ali == wanted_ali

def test_read_ndx(tmp_path):
    indexfile = tmp_path / "param.txt"
    for f in [indexfile]:
        f.write_text('[ upper leaflet ]\n1 2 3 4 5 6 7 8 9 10\n[ lower leaflet ]\n11 12 13 14 15 16 17 18 19 20\n')
    tested_output = p.read_ndx(indexfile)
    wanted_outputUp = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    wanted_outputLo = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    np.testing.assert_array_equal(tested_output[0], wanted_outputUp)
    np.testing.assert_array_equal(tested_output[1], wanted_outputLo)