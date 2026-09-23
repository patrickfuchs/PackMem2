from pathlib import Path
import pytest
import packmem2.concatenate as concatenate

def test_concatenate_keep(tmp_path):
    """DMPG (Charmm)"""
    path = "tests/data/concatenate_data"
    output_dir = str(tmp_path)
    prefix = "DMPG"
    start = 0
    end = 10
    protein = False
    keep = True
    

    concatenate.launch(
        path,
        output_dir,
        prefix,
        start,
        end,
        protein,
        keep,
    )

    # Check content file
    expected_output_Deep = Path(f"{output_dir}/Total_Deep.csv")
    with open("tests/data/concatenate_data/Total_Deep.csv", "r") as f_in:
        expected_content_Deep = f_in.read()
    assert expected_output_Deep.exists()
    assert expected_output_Deep.read_text() == expected_content_Deep

    # Check check the we kept the input files
    expected_input_Deep0_Up = Path(f"{path}/DMPG0_Up_Deep_result.txt")
    expected_input_Shallow5_Up = Path(f"{path}/DMPG5_Up_Shallow_result.txt")
    expected_input_All10_Lo = Path(f"{path}/DMPG10_Lo_All_result.txt")
    assert expected_input_Deep0_Up.exists()
    assert expected_input_Shallow5_Up.exists()
    assert expected_input_All10_Lo.exists()
    


def test_concatenate_keep(tmp_path):
    """DMPG (Charmm)"""
    path = "tests/data/concatenate_data"
    output_dir = str(tmp_path)
    prefix = "DMPG"
    start = 0
    end = 10
    protein = False
    keep = True
    

    concatenate.launch(
        path,
        output_dir,
        prefix,
        start,
        end,
        protein,
        keep,
    )

    # Check Total files
    expected_output_Deep = Path(f"{output_dir}/Total_Deep.csv")
    expected_output_Shallow = Path(f"{output_dir}/Total_Shallow.csv")
    expected_output_All = Path(f"{output_dir}/Total_All.csv")
    with open("tests/data/concatenate_data/Total_Deep.csv", "r") as f_in:
        expected_content_Deep = f_in.read()
    assert expected_output_Deep.exists()
    assert expected_output_Shallow.exists()
    assert expected_output_All.exists()
    assert expected_output_Deep.read_text() == expected_content_Deep

    # Check Total_Up files
    expected_output_Deep_Up = Path(f"{output_dir}/Total_Up_Deep.csv")
    expected_output_Shallow_Up = Path(f"{output_dir}/Total_Up_Shallow.csv")
    expected_output_All_Up = Path(f"{output_dir}/Total_Up_All.csv")
    assert expected_output_Deep_Up.exists()
    assert expected_output_Shallow_Up.exists()
    assert expected_output_All_Up.exists()

    # Check Total_Lo files
    expected_output_Deep_Lo = Path(f"{output_dir}/Total_Lo_Deep.csv")
    expected_output_Shallow_Lo = Path(f"{output_dir}/Total_Lo_Shallow.csv")
    expected_output_All_Lo = Path(f"{output_dir}/Total_Lo_All.csv")
    assert expected_output_Deep_Lo.exists()
    assert expected_output_Shallow_Lo.exists()
    assert expected_output_All_Lo.exists()
