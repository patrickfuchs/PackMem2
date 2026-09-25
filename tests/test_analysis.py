from pathlib import Path
import packmem2.analysis as analysis


def test_analysis(tmp_path):
    """DMPG_prot (Charmm)"""
    input_path = "tests/data/analysis_data"
    output_dir = str(tmp_path)
    output_name = "DMPG"
    protein = True
    nb_block = 3
    limx = 15
    limy = 1e-4
    precision = 2
    

    analysis.launch(
        input_path,
        output_dir,
        output_name,
        protein,
        nb_block,
        limx,
        limy,
        precision,
    )

    # Check pdf
    expected_final_output = Path(f"{output_dir}/{output_name}.pdf")
    expected_final_summary = Path(f"{output_dir}/{output_name}.csv")
    with open("tests/data/analysis_data/DMPG.csv", "r") as f_in:
        expected_content_summary = f_in.read()

    assert expected_final_output.exists()
    assert expected_final_summary.exists()
    assert expected_final_summary.read_text() == expected_content_summary
