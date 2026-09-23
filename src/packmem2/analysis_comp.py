# M. Zygadlo 2025

import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def get_arguments() -> argparse.Namespace:
    """
    Get the arguments for the script and check that the inputfiles are valid.

    --------------------
    OUTPUT
    parser.parse_args
    """
    # Getting the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f1', action = 'store', dest = 'file1',
        help = 'The path towards the first PackMem results')
    parser.add_argument('-f2', action = 'store', dest = 'file2',
        help = 'The path towards the second PackMem results')
    parser.add_argument('-f3', action = 'store', dest = 'file3',
        default=None,
        help = 'The path towards the third PackMem results')
    parser.add_argument('-lg', action = 'store', dest = 'legend',
        help = 'The legend of the graph. Separated by :')
    parser.add_argument("-prot", action="store_true", dest="prot",
        help="If there is a protein")
    parser.add_argument("-w_wo", action="store_true", dest="w_wo",
            help="If there is a file with a prot and the other without")
    parser.add_argument("-o", action="store", dest="output",
        default="Res_comp",
        help="The name of the output .pdf file. Default = Res_comp")
    parser.add_argument("-od", action="store", dest="output_dir",
        default="./",
        help="Name for output directory (default: ./)")
    args = parser.parse_args()

    return args


def read_file(filename: str, prot: bool = False) -> pd.DataFrame:
    """
    Read the Packing defect constants file

    --------------------
    INPUT
    filename : string
        The name of the packing defect area file.
    prot : boolean
        If the file is the one with the protein packing defect area
    
    --------------------
    OUTPUT
    Pandas DataFrame
        Contains the packing defect constants
    """
    if Path(filename).is_file():
        if filename[-3:] == "csv":
            def_csts = pd.read_csv(filename, index_col=0)
        elif filename[-3:] == "txt":
            csts = []
            with open(filename, "r") as file_in:
                for line in file_in:
                    if "Using" in line:
                        csts.append(float(line.split(' ')[4]))
                    elif "Mean" in line:
                        csts.append(float(line.split(' ')[7]))
                        csts.append(float(line.split(' ')[9]))
            # Reshape to form 3 columns of csts
            def_csts = pd.DataFrame(np.array(csts).reshape(3, 6)).T
            def_csts.columns = ["Deep_Total", "Shallow_Total", "All_Total"]
            def_csts.index = ["cst_global", "cst_block1", "cst_block2", "cst_block3", "cst_mean_blocks", "error_mean_blocks"]
    else:
        sys.exit(f"file {filename} not found")

    return def_csts


def change_colnames(
        legend: str,
        name: str, packdef_csts: pd.DataFrame, errors: pd.Series,
) -> tuple[pd.DataFrame, pd.Series]:
    """
    Change the column names of the dtf to have the wanted legend

    --------------------
    INPUT
    legend: string
        The legend fo the plot given by the user
        Of the form XX:YY:ZZ 
    name: string
        the type to analyse: Total/Total_Up/Total_Lo
    packdef_csts: pandas DataFrame
        Contains the statistics of the defects
    errors: pandas Series
        Contains the errors of each defect

    --------------------
    OUTPUT
    tuple
        pandas DataFrame
            Has new columns names
        pandas Series
            Has new columns names
    """
    # Create the dataframe needed for the graph
    # A mix of the first and the second dtf ordered by defect type
    cst_packing = pd.DataFrame(
        packdef_csts[[f"Deep_{name}", f"Shallow_{name}", f"All_{name}"]])
    cst_packing.columns = [f"Deep\n{legend}", f"Shallow\n{legend}", f"All\n{legend}"]
    errors.index = [f"Deep\n{legend}", f"Shallow\n{legend}", f"All\n{legend}"]

    return cst_packing, errors


def plot_comp_defect_constants(
    pdf: PdfPages, legend: str,
    name1: str, packdef_csts1: pd.DataFrame, errors1: pd.Series,
    name2: str, packdef_csts2: pd.DataFrame, errors2: pd.Series,
    name3: str = None, packdef_csts3: pd.DataFrame = None, errors3: pd.Series = None,
    w_wo: bool = False,
) -> None:
    """
    Plot the defects constants for each defect.

    --------------------
    INPUT
    pdf: matplotlib.backends.backend_pdf.PdfPages
        Contains the figures in the final pdf
    legend: string
        The legend fo the plot given by the user
        Of the form XX:YY:ZZ 
    name{1..3}: string
        the type to analyse: Total/Total_Up/Total_Lo
    packdef_csts{1..3}: pandas DataFrame
        Contains the statistics of the defects
    errors{1..3}: pandas Series
        Contains the errors of each defect
    w_wo: boolean
        If the user wants to compare a case with and without a protein
    """
    # Get the legend
    legend = legend.split(':')
    for i in range(len(legend)):
        legend[i] = legend[i].replace('_', '\n')
    if w_wo:
        legend.append(legend[1]+"\nfar")
        legend[1] = legend[1]+"\nclose"


    cst_packing1, errors1 = change_colnames(legend[0], name1, packdef_csts1, errors1)
    cst_packing2, errors2 = change_colnames(legend[1], name2, packdef_csts2, errors2)

    if packdef_csts3 is not None:
        cst_packing3, errors3 = change_colnames(legend[2], name3, packdef_csts3, errors3)

        # Concatenate the dtf and errors
        cst_packing_comp = pd.concat([cst_packing1, cst_packing2, cst_packing3], axis = 1).iloc[:, [0, 3, 6, 1, 4, 7, 2, 5, 8]]
        errors_comp = pd.concat([errors1, errors2, errors3]).iloc[[0, 3, 6, 1, 4, 7, 2, 5, 8]]

        # Get the colour, text position and bars to compare
        colour = ["darkred", "firebrick", "lightcoral", "darkgreen", "forestgreen", "darkseagreen", "royalblue", "cornflowerblue", "lightsteelblue"]
        text_pos = [0.0, 1.01, 2.01, 3.01, 4.01, 5.01, 6.01, 7.01, 8.01]
        comparison_bar = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    else:
        cst_packing_comp = pd.concat([cst_packing1, cst_packing2], axis = 1).iloc[:, [0, 3, 1, 4, 2, 5]]
        errors_comp = pd.concat([errors1, errors2]).iloc[[0, 3, 1, 4, 2, 5]]

        colour = ["darkred", "firebrick", "darkgreen", "forestgreen", "royalblue", "cornflowerblue"]
        text_pos = [0.0, 1.01, 2.01, 3.01, 4.01, 5.01]
        comparison_bar = [[0, 1], [2, 3], [4, 5]]


    # Create the figure
    fig, ax = plt.subplots()

    # Get the variables needed for the graph
    x_index = cst_packing_comp.columns
    packdef_values = cst_packing_comp.loc["cst_mean_blocks"]

    ax.bar(x_index, packdef_values, color=colour, width=0.5, yerr=errors_comp, error_kw=dict(ecolor='darkgrey'), capsize=3)

    if packdef_csts3 is not None:
        # Add the text for the packing constant and the percentage:
        for i, defect in enumerate(packdef_values.index):
            if i < 3:
                j=0
            elif i < 6:
                j=1
            else:
                j=2
            ax.text(text_pos[i],
                packdef_values.loc[f"{defect}"] + (max(packdef_values)/100),
                f"{packdef_values.loc[f'{defect}'] / max(packdef_values.iloc[comparison_bar[j]]) * 100:.1f}%",
                verticalalignment="bottom", horizontalalignment="center",
                fontsize=8)
            ax.text(
                text_pos[i],
                packdef_values.loc[f"{defect}"] + (7*max(packdef_values)/100),
                f"{packdef_values.loc[f'{defect}']:.1f} $Å^2$",
                verticalalignment="bottom", horizontalalignment="center",
                fontsize=8)
    else:
        for i, defect in enumerate(packdef_values.index):
            if i < 2:
                j=0
            elif i < 4:
                j=1
            else:
                j=2
            ax.text(text_pos[i],
                packdef_values.loc[f"{defect}"] + (max(packdef_values)/100),
                f"{packdef_values.loc[f'{defect}'] / max(packdef_values.iloc[comparison_bar[j]]) * 100:.1f}%",
                verticalalignment="bottom", horizontalalignment="center",
                fontsize=8)
            ax.text(text_pos[i],
                packdef_values.loc[f"{defect}"] + (7*max(packdef_values)/100),
                f"{packdef_values.loc[f'{defect}']:.1f} $Å^2$",
                verticalalignment="bottom", horizontalalignment="center",
                fontsize=8)
    ax.set_ylabel("Defect size constant ${Å^2}$")
    ax.set_ylim(0, int(max(packdef_values)) + (20*max(packdef_values)/100))
    ax.set_title(f"{name1} Packing defect constants computed by block averaging")
    plt.tight_layout()
    pdf.savefig()  # Save the current figure to the PDF
    plt.close()


def launch(
    file1: str, file2: str, file3: str, legend: str,
    output_dir: str, output: str, prot: bool,
    w_wo: bool,
) -> None:
    """
    Analyse the defect size computed by PackMem2
    """
    # Open a pdf device
    # Create a single PDF file
    pdf = PdfPages(f"{output_dir}/{output}.pdf")

    packdef_csts_1 = read_file(file1)
    packdef_csts_2 = read_file(file2)
    if packdef_csts_1.shape[1] == 3 or packdef_csts_2.shape[1] == 3:
        list_name = ["Total"]
    else:
        list_name = ["Total", "Total_Up", "Total_Lo"]

    if file3 is not None:
        packdef_csts_3 = read_file(file3)
        if packdef_csts_1.shape[1] == 3 or \
            packdef_csts_2.shape[1] == 3 or \
            packdef_csts_3.shape[1] == 3:
            list_name = ["Total"]
        else:
            list_name = ["Total", "Total_Up", "Total_Lo"]
    
    for name in list_name:
        # Plot the packing defects constants + errors
        # (for each packdef) on a barplot
        errors1 = packdef_csts_1.loc["error_mean_blocks", [f"Deep_{name}", f"Shallow_{name}", f"All_{name}"]]
        errors2 = packdef_csts_2.loc["error_mean_blocks", [f"Deep_{name}", f"Shallow_{name}", f"All_{name}"]]

        if file3 is not None:
            errors3 = packdef_csts_3.loc["error_mean_blocks", [f"Deep_{name}", f"Shallow_{name}", f"All_{name}"]]
            plot_comp_defect_constants(pdf, legend,
                name, packdef_csts_1, errors1,
                name, packdef_csts_2, errors2,
                name, packdef_csts_3, errors3)
        else:
            plot_comp_defect_constants(pdf, legend,
                name, packdef_csts_1, errors1,
                name, packdef_csts_2, errors2)

    if prot:
        for name in ["Total_Up", "Total_Lo"]:
            if (packdef_csts_1[
                    [f"Deep_{name}_close", f"Shallow_{name}_close", f"All_{name}_close",
                    f"Deep_{name}_far", f"Shallow_{name}_far", f"All_{name}_far"]]
                    .isna().all().all()):
                continue
            elif (packdef_csts_2[
                    [f"Deep_{name}_close", f"Shallow_{name}_close", f"All_{name}_close",
                    f"Deep_{name}_far", f"Shallow_{name}_far", f"All_{name}_far"]]
                    .isna().all().all()):
                continue
            # Plot the final decays + errors computed with block averaging
            # (for each packdef) on a barplot
            errors1 = packdef_csts_1.loc["error_mean_blocks"]
            errors2 = packdef_csts_2.loc["error_mean_blocks"]

            if file3 is not None:
                if (packdef_csts_3[
                        [f"Deep_{name}_close", f"Shallow_{name}_close", f"All_{name}_close",
                        f"Deep_{name}_far", f"Shallow_{name}_far", f"All_{name}_far"]]
                        .isna().all().all()):
                    continue
                errors3 = packdef_csts_3.loc["error_mean_blocks"]
                plot_comp_defect_constants(pdf, legend,
                    f"{name}_close", packdef_csts_1, errors1.loc[[f"Deep_{name}_close", f"Shallow_{name}_close", f"All_{name}_close"]],
                    f"{name}_close", packdef_csts_2, errors2.loc[[f"Deep_{name}_close", f"Shallow_{name}_close", f"All_{name}_close"]],
                    f"{name}_close", packdef_csts_3, errors3.loc[[f"Deep_{name}_close", f"Shallow_{name}_close", f"All_{name}_close"]])
                plot_comp_defect_constants(pdf, legend,
                    f"{name}_far", packdef_csts_1, errors1.loc[[f"Deep_{name}_far", f"Shallow_{name}_far", f"All_{name}_far"]],
                    f"{name}_far", packdef_csts_2, errors2.loc[[f"Deep_{name}_far", f"Shallow_{name}_far", f"All_{name}_far"]],
                    f"{name}_far", packdef_csts_3, errors3.loc[[f"Deep_{name}_far", f"Shallow_{name}_far", f"All_{name}_far"]])
            else:
                plot_comp_defect_constants(pdf, legend,
                    f"{name}_close", packdef_csts_1, errors1.loc[[f"Deep_{name}_close", f"Shallow_{name}_close", f"All_{name}_close"]],
                    f"{name}_close", packdef_csts_2, errors2.loc[[f"Deep_{name}_close", f"Shallow_{name}_close", f"All_{name}_close"]])
                plot_comp_defect_constants(pdf, legend,
                    f"{name}_far", packdef_csts_1, errors1.loc[[f"Deep_{name}_far", f"Shallow_{name}_far", f"All_{name}_far"]],
                    f"{name}_far", packdef_csts_2, errors2.loc[[f"Deep_{name}_far", f"Shallow_{name}_far", f"All_{name}_far"]])

    if w_wo:
        # if the first file given is the one with the protein
        if packdef_csts_1.shape[1] == 21:
            PD_csts1 = packdef_csts_2
            PD_csts2 = packdef_csts_1
        elif packdef_csts_2.shape[1] == 21:
            PD_csts1 = packdef_csts_1
            PD_csts2 = packdef_csts_2

        errors1 = PD_csts1.loc["error_mean_blocks"]
        errors2 = PD_csts2.loc["error_mean_blocks"]
        for name in ["Total_Up", "Total_Lo"]:
            if (PD_csts2[
                    [f"Deep_{name}_close", f"Shallow_{name}_close", f"All_{name}_close",
                    f"Deep_{name}_far", f"Shallow_{name}_far", f"All_{name}_far"]]
                    .isna().all().all()):
                continue

            plot_comp_defect_constants(pdf, legend,
                name, PD_csts1, errors1.loc[['Deep_Total_Up', 'Shallow_Total_Up', 'All_Total_Up']],
                f"{name}_close", PD_csts2, errors2.loc[[f"Deep_{name}_close", f"Shallow_{name}_close", f"All_{name}_close"]],
                f"{name}_far", PD_csts2, errors2.loc[[f"Deep_{name}_far", f"Shallow_{name}_far", f"All_{name}_far"]],
                w_wo)
    
    # Close the PDF file
    pdf.close()


def main() -> None:
    # Get arguments
    args = get_arguments()

    launch(args.file1, args.file2, args.file3, args.legend,
        args.output_dir, args.output, args.prot,
        args.w_wo)


if __name__ == "__main__":
    main()
