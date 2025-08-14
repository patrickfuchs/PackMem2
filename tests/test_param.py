import pytest
import sys
import os
import re

# Add the parent rep to the sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from bin import param as p

def test_file_present():
    # Test if file present
    filename_OK = "src/PackMem_prot.py"
    p.file_present(str(filename_OK))

    # Test if file absent
    filename_NO = "piou.txt"
    with pytest.raises(FileNotFoundError):
        p.file_present(str(filename_NO))

def test_get_args(monkeypatch, tmp_path):
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
            "-p", str(param),
            "-r", str(radii),
            "-d", "1.0",
            "-prot"
        ]
    )
    args = p.get_args()
    assert args.traj == str(traj)
    assert args.topo == str(topo)
    assert args.lipid == "DPPC"
    assert args.dist_suppl_Z == 1.0
    assert args.protein is True

    ## TEST NO ##
    monkeypatch.setattr(
        sys, "argv",
        [
            "prog",
            "-f", str(traj),
            "-s", str(topo),
            "-l", "DPPC",
            "-p", str(param),
            "-r", str(radii),
            "-d", "-1.0",
        ]
    )
    with pytest.raises(Exception) as exc_info:
        p.get_args()
    assert "must be > 0.0" in str(exc_info.value)

def test_dict_2columns():
    list_str = ["DOPC C2", "DMPG C2", "POPC C2"]
    tested_ouput = p.dict_2columns(list_str)
    wanted_output = {"DOPC": "C2", "DMPG": "C2", "POPC": "C2"}
    assert tested_ouput == wanted_output