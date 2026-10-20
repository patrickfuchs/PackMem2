from pathlib import Path
import packmem2.calculate_radii as calculate_radii


def test_calculate_radii_martini(tmp_path):
    """POPC_DMPS (Martini2)"""
    mol_file = "tests/data/calculate_radii_data/martini_v2.0_POPC_DMPS.itp"
    ff_file = "tests/data/calculate_radii_data/martini_v2.2.itp"
    output = str(tmp_path) + "/out"
    forcefield = "martini"
    

    calculate_radii.launch(
        ff_file,
        mol_file,
        forcefield,
        output
    )

    # Check pdf
    expected_final_output = Path(f"{output}.txt")
    with open("tests/data/calculate_radii_data/martini2_POPC_DMPS.txt", "r") as f_in:
        expected_content_summary = f_in.read()

    assert expected_final_output.exists()
    assert expected_final_output.read_text() == expected_content_summary


def test_calculate_radii_martini(tmp_path):
    """DOPC (Charmm)"""
    mol_file = "tests/data/calculate_radii_data/charmm36_DOPC.itp"
    ff_file = "tests/data/calculate_radii_data/charmm36.itp"
    output = str(tmp_path) + "/out"
    forcefield = "charmm"
    

    calculate_radii.launch(
        ff_file,
        mol_file,
        forcefield,
        output
    )

    # Check pdf
    expected_final_output = Path(f"{output}.txt")
    with open("tests/data/calculate_radii_data/charmm36_DOPC.txt", "r") as f_in:
        expected_content_summary = f_in.read()

    assert expected_final_output.exists()
    assert expected_final_output.read_text() == expected_content_summary