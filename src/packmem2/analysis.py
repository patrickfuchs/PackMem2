# M. Zygadlo 2025

import argparse
import os
import pandas as pd
import numpy as np
import math
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
    parser.add_argument('-prot', action = 'store_true', dest = 'prot',
                    help = 'If there is a protein')
    parser.add_argument('-p', action = 'store', dest = 'precision',
                    type=int, default=2,
                    help = 'The precision for writing packdef constants (nb of decimals) in the output. Default = 2')
    parser.add_argument('-lx', action = 'store', dest = 'limx',
                    type=int, default=15,
                    help = 'The lowest defect area used for the fit (we recommand not to touch to this value). Default = 15')
    parser.add_argument('-ly', action = 'store', dest = 'limy',
                    type=float, default=1e-4,
                    help = 'The lowest probability used for the fit (we recommand not to touch to this value). Default = 1e-4')
    parser.add_argument('-o', action = 'store', dest = 'output',
                    default="Res_membrane",
                    help = 'The name of the output .pdf file. Default = Res_membrane')
    parser.add_argument('-od', action='store', dest='output_dir',
                        default = './',
                        help = 'Name for output directory (default: ./)')
    args = parser.parse_args()

    return args

def log(
    y_list: list
    ) -> list:
    """
    Calculate the logarithm of a list.
    
    --------------------
    INPUT
    y_list : list
        The y values.
    
    --------------------
    OIUTPUT
    list
        The y values that have been through the logarithm.
    """
    log_y = []
    for y_ind in y_list:
        if y_ind > 0:
            log_y.append(math.log(y_ind))
        else:
            log_y.append(0)
    return log_y

def fit_decay(
    x: np.array,
    y: np.array,
    limx: int | float,
    limy: float
    ) -> np.array:
    """
    Function does a linear fit.
    
    The linear fit is made on the probability of having a certain packing area.

    --------------------
    INPUT
    x : numpy array
        Packing area.
    y : numpy array
        probability of having a certain packing area.
    limx : int
        The lowest defect area used for the fit.
    limy : float
        The lowest probability used for the fit
    
    --------------------
    OUTPUT
    numpy array
        A linear fit of x and y.
    """
    # fit with defects above LIMX nm and proba > LIMY
    y = y[x >= limx]
    x = x[x >= limx]
    x = x[y >= limy]
    y = y[y >= limy]
    FIT = np.polyfit(x, log(y), 1)
    return (FIT)

def block_averaging(
    vect: pd.DataFrame,
    nb_block: int,
    limx: int | float,
    limy: float
    ) -> list[float, float, float]:
    """
    Divide the packing data into n blocks.
    
    --------------------
    INPUT
    vect : pandas dataframe
        The size of the packing defect.
    nb_block : int
        The number of blocks we want to have.
    limx : int
        The lowest defect area used for the fit.
    limy : float
        The lowest probability used for the fit.
    
    --------------------
    OUTPUT
    list
        A vector of 3 decays.
    """
    bornes = [int((len(vect)/3)*nb) for nb in range(nb_block+1)]
    decays = []
    for i in range(0,nb_block):
        subvect = vect[bornes[i]:bornes[i+1]]
        H = plt.hist(subvect, bins = np.arange(0.5, max(vect)+0.5))
        x = H[1]+0.5
        x = x[:len(x)-1]
        y = H[0]/sum(H[0])
        FIT = fit_decay(x, y, limx, limy)
        decays.append(abs(1/FIT[0]))
    return decays

def plot_defect_fit(
    type: str,
    defect: str,
    packdef_data: pd.DataFrame,
    packdef_constants: pd.DataFrame,
    limx: int | float,
    limy: float,
    precision: int,
    pdf: matplotlib.backends.backend_pdf.PdfPages
    ) -> pd.DataFrame:
    """
    Compute and plot the fit of the defects distribution.

    --------------------
    INPUT
    type: string
        The type to analyse: Total/Total_Up/Total_Lo
    defect: string
        The type  of  defect: Deep/Shallow/All
    packdef_data: pandas DataFrame
        Contains the sizes of the defects
    packdef_constants: pandas DataFrame
        Will contains the statistics of the defects found in this function
    limx : int
        The lowest defect area used for the fit
    limy : float
        The lowest probability used for the fit
    precision : int
        The precision for writing packdef constants
    pdf : matplotlib.backends.backend_pdf.PdfPages
        Contains the figures in the final pdf

    --------------------
    OUTPUT
    pandas DataFrame
        Contains the statistics of the defects found in this function
    """  
    # Compute PackDef distributions (on the whole set)
    H = plt.hist(packdef_data, bins=np.arange(0.5,max(packdef_data)+0.5))
    # Length of the defects
    x = np.arange(0, max(packdef_data)-1)
    # nb of observations of a certain defect length
    y = H[0]/sum(H[0])
    FIT = fit_decay(x,y, limx, limy)
    fit_function = np.poly1d(FIT)

    ### Plot for defect fit
    plt.clf()
    plt.scatter(x, log(y), marker='o', facecolor="none", edgecolor="black")
    plt.xlim(-2, 102)
    plt.ylim(-10, -3.5)
    plt.yticks([math.log(1e-4), math.log(1e-3), math.log(1e-2)], labels=[str(1e-4), str(1e-3), str(1e-2)])
    plt.ylabel("Probability")
    plt.xlabel("Defect area ${A^2}$")
    plt.title(f"{defect} {type}")
    plt.axvline(limx, color='gray', linestyle='--')
    plt.axhline(math.log(limy), color='gray', linestyle='-')
    plt.plot(x, fit_function(x), color='red', label='Fit')
    pdf.savefig()  # Save the current figure to the PDF
    plt.close()

    global_inv_decay = abs(1/FIT[0])
    
    # print global results
    print("")
    print(f"Results for {defect} {type} defects")
    
    # Compute PackDef distributions with block averaging method, and print results
    FITS_3blocks=block_averaging(packdef_data, 3, limx, limy)
    for i in range(3):
        print(f"Using block {i+1}: {round(FITS_3blocks[i], precision)} Å\u00b2")
    # Compute the mean of the 3 blocks averages
    inv_decay_block = sum(FITS_3blocks)/len(FITS_3blocks)
    # Compute the standard deviation of the 3 values
    error_inv_decay_block = np.std(FITS_3blocks)
    print(f"Mean +/- sd on 3 blocks: {round(inv_decay_block, precision)} ± {round(error_inv_decay_block, precision)} Å\u00b2")

    # store all the results in the data frame
    packdef_constants.loc["PackDef_cst_global", f"{defect}_{type}"] = global_inv_decay
    packdef_constants.loc["PackDef_cst_block1", f"{defect}_{type}"] = FITS_3blocks[0]
    packdef_constants.loc["PackDef_cst_block2", f"{defect}_{type}"] = FITS_3blocks[1]
    packdef_constants.loc["PackDef_cst_block3", f"{defect}_{type}"] = FITS_3blocks[2]
    packdef_constants.loc["PackDef_cst_all_blocks", f"{defect}_{type}"] = inv_decay_block
    packdef_constants.loc["error_all_blocks", f"{defect}_{type}"] = error_inv_decay_block
    
    return packdef_constants

def plot_defect_constants_blocks(
    type: str,
    packdef_constants: pd.DataFrame,
    errors: pd.Series,
    pdf: matplotlib.backends.backend_pdf.PdfPages
    ) -> None:
    """
    Plot the defect constants of all block and their average.

    --------------------
    INPUT
    type: string
        the type to analyse: Total/Total_Up/Total_Lo
    packdef_constants: pandas DataFrame
        Will contains the statistics of the defects found in this function
    errors: pandas Series
        Contains the errors of each defect
    pdf : matplotlib.backends.backend_pdf.PdfPages
        Contains the figures in the final pdf
    """    
    # Plot and save the second figure (bar plot for just one row)
    packdef_constants[[f"Deep_{type}", f"Shallow_{type}", f"All_{type}"]][:5].T.plot.bar(color=["darkred", "firebrick", "indianred", "lightcoral", "rosybrown"], yerr=errors, capsize=3, rot=0)
    plt.ylabel("Defect size constant ${A^2}$")
    plt.ylim(0, int(max(packdef_constants[[f"Deep_{type}", f"Shallow_{type}", f"All_{type}"]][:5].max()))+2)
    plt.title("Packing defect constants")
    pdf.savefig()  # Save the current figure to the PDF
    plt.close()

def get_outliers(
    dtf_packdef_values: pd.DataFrame,
    type: str
    ) -> pd.DataFrame:
    """
    Get the outliers in a dataframe

    --------------------
    INPUT
    dtf_packdef_values: pandas DataFrame
        Contains the informations on the packing defects (mean, data, error)
    type: string
        the type to analyse: Total/Total_Up/Total_Lo

    --------------------
    OUTPUT
    pandas DataFrame
        Contains the x and y values of the outliers
    """
    # Prepare empty dtf
    outliers = pd.DataFrame(columns = ['outliers', 'x_index'],)
    list_outliers = []
    list_index = []

    for defect in ['Deep', 'Shallow', 'All']:
        # Compute the maximum and minimum values for the standard deviation
        mean_top = dtf_packdef_values.loc['PackDef_cst_all_blocks', f'{defect}_{type}'] + dtf_packdef_values.loc['error_all_blocks', f'{defect}_{type}']
        mean_bot = dtf_packdef_values.loc['PackDef_cst_all_blocks', f'{defect}_{type}'] - dtf_packdef_values.loc['error_all_blocks', f'{defect}_{type}']

        sub_dtf = dtf_packdef_values.loc[['PackDef_cst_block1', 'PackDef_cst_block2', 'PackDef_cst_block3'], f'{defect}_{type}']
        list_outliers += list(sub_dtf[(sub_dtf < mean_bot) | (sub_dtf > mean_top)])

        if defect == 'Deep':
            list_index+= [0] * len(sub_dtf[(sub_dtf < mean_bot) | (sub_dtf > mean_top)])
        elif defect == 'Shallow':
            list_index += [1] * len(sub_dtf[(sub_dtf < mean_bot) | (sub_dtf > mean_top)])
        else:
            list_index+= [2] * len(sub_dtf[(sub_dtf < mean_bot) | (sub_dtf > mean_top)])

    # Complete the final dtf
    outliers['outliers'] = np.array(list_outliers)
    outliers['x_index'] = np.array(list_index)

    return outliers

def plot_defect_constants(
    type: str,
    packdef_constants: pd.DataFrame,
    errors: pd.Series,
    pdf: matplotlib.backends.backend_pdf.PdfPages
    ) -> None:
    """
    Plot the defects constants for each defect.

    --------------------
    INPUT
    type: string
        the type to analyse: Total/Total_Up/Total_Lo
    packdef_constants: pandas DataFrame
        Will contains the statistics of the defects found in this function
    errors: pandas Series
        Contains the errors of each defect
    pdf : matplotlib.backends.backend_pdf.PdfPages
        Contains the figures in the final pdf
    """
    cst_packing = pd.DataFrame(packdef_constants[[f"Deep_{type}", f"Shallow_{type}", f"All_{type}"]])

    # Create the figure
    fig, ax = plt.subplots()

    # Get the variables needed for the graph
    x_index = cst_packing.columns
    packdef_values = cst_packing.loc['PackDef_cst_all_blocks']
    colour = ["firebrick", "forestgreen", "royalblue"]
    outliers = get_outliers(cst_packing, type)
    text_pos = [0.0, 1.01, 2.01]

    ax.bar(x_index, packdef_values, color = colour, width = 0.5, yerr = errors, capsize=3)
    plt.scatter(x= outliers['x_index'], y=outliers['outliers'], s=15, c='black', alpha=0.5)
    # Add the text for the packing constant and the percentage
    for i, defect in enumerate(['Deep', 'Shallow', 'All']):
        ax.text(text_pos[i], 0.25, f'{packdef_values.loc[f"{defect}_{type}"]/max(packdef_values)*100:.1f}%',
            verticalalignment='bottom', horizontalalignment='center',
            fontsize=12)
        ax.text(text_pos[i], 1.4, f'{packdef_values.loc[f"{defect}_{type}"]:.1f} $A^2$',
            verticalalignment='bottom', horizontalalignment='center',
            fontsize=12)
    ax.set_ylabel("Defect size constant ${A^2}$")
    ax.set_ylim(0, int(max(packdef_values))+2)
    ax.set_title("Packing defect constants computed by block averaging")
    pdf.savefig()  # Save the current figure to the PDF
    plt.close()

def launch(
    output_dir: str,
    output: str,
    prot: bool,
    limx: int | float,
    limy: float,
    precision: int
    ) -> None:
    """
    Analyse the defect size computed by PackMem2
    """
    # Open a pdf device
    # Create a single PDF file
    pdf = PdfPages(f'{output_dir}/{output}.pdf')

    columns_dtf =  ['Deep_Total', 'Shallow_Total', 'All_Total', 'Deep_Total_Up', 'Shallow_Total_Up', 'All_Total_Up', 'Deep_Total_Lo', 'Shallow_Total_Lo', 'All_Total_Lo']
    index_dtf = ["PackDef_cst_global", "PackDef_cst_block1", "PackDef_cst_block2", "PackDef_cst_block3","PackDef_cst_all_blocks","error_all_blocks"]
    if prot:
        columns_dtf += ['Deep_Total_Up_close', 'Shallow_Total_Up_close', 'All_Total_Up_close', 'Deep_Total_Up_far', 'Shallow_Total_Up_far', 'All_Total_Up_far', 'Deep_Total_Lo_close', 'Shallow_Total_Lo_close', 'All_Total_Lo_close', 'Deep_Total_Lo_far', 'Shallow_Total_Lo_far', 'All_Total_Lo_far']
    # Initialize a data frame to store packdef constants + errors
    packdef_constants = pd.DataFrame(columns = columns_dtf, index = index_dtf)

    for name in ["Total", "Total_Up", "Total_Lo"]:
        # Now loop over the three default types
        for defect in ["Deep","Shallow","All"]:
            filename = f"{output_dir}/{name}_{defect}.csv"
            def_area = pd.read_csv(filename, header=None)[1]

            # Plot the fit of the defects distribution
            packdef_constants = plot_defect_fit(name, defect, def_area, packdef_constants, limx, limy, precision, pdf)
    
        # Plot all packdef constants on a single barplot
        # Allows to estimate the relative convergence of the simulation
        errors=packdef_constants.T.error_all_blocks.to_frame('PackDef_cst_all_blocks')
        plot_defect_constants_blocks(name, packdef_constants, errors, pdf)
        
        # Plot the final decays + errors computed with block averaging
        # (for each packdef) on a barplot
        plot_defect_constants(name, packdef_constants, errors.loc[[f'Deep_{name}', f'Shallow_{name}', f'All_{name}'], 'PackDef_cst_all_blocks'], pdf)
    
    if prot:
        for name in ["Total_Up", "Total_Lo"]:
            # Now loop over the three default types
            for defect in ["Deep","Shallow","All"]:
                # Load PackMem data
                filename = f"{output_dir}/{name}_{defect}_prot.csv"
                if not Path(filename).is_file():
                    continue
                def_area_prot = pd.read_csv(filename, header=None).iloc[:,1:]
                def_area_prot_close = def_area_prot[def_area_prot[1] == "close"][2]
                def_area_prot_far = def_area_prot[def_area_prot[1] == "far"][2]
                
                # Plot the fit of the defects distribution
                packdef_constants = plot_defect_fit(f'{name}_close', defect, def_area_prot_close, packdef_constants, limx, limy, precision, pdf)
                packdef_constants = plot_defect_fit(f'{name}_far', defect, def_area_prot_far, packdef_constants, limx, limy, precision, pdf)
            if packdef_constants[[f'Deep_{name}_close', f'Shallow_{name}_close', f'All_{name}_close', f'Deep_{name}_far', f'Shallow_{name}_far', f'All_{name}_far']].isna().all().all():
                continue
            # Plot all packdef constants on a single barplot
            # Allows to estimate the relative convergence of the simulation
            errors=packdef_constants.T.error_all_blocks.to_frame('PackDef_cst_all_blocks')
            plot_defect_constants_blocks(f'{name}_close', packdef_constants, errors, pdf)
            plot_defect_constants_blocks(f'{name}_far', packdef_constants, errors, pdf)
            
            # Plot the final decays + errors computed with block averaging
            # (for each packdef) on a barplot
            plot_defect_constants(f'{name}_close', packdef_constants, errors.loc[[f'Deep_{name}_close', f'Shallow_{name}_close', f'All_{name}_close'], 'PackDef_cst_all_blocks'], pdf)
            plot_defect_constants(f'{name}_far', packdef_constants, errors.loc[[f'Deep_{name}_far', f'Shallow_{name}_far', f'All_{name}_far'], 'PackDef_cst_all_blocks'], pdf)

    # Close the PDF file
    pdf.close()

def main() -> None:
    # Get arguments
    args = get_arguments()

    launch(args.output_dir, args.output, args.prot, args.limx, args.limy, args.precision)

if __name__=="__main__":
    main()
