# M. Zygadlo 2025

import argparse
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def get_arguments():
    """
    Get the arguments for the script and check that the inputfiles are valid.

    --------------------
    OUTPUT
    parser.parse_args
    """
    # Getting the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', action = 'store', dest = 'lipid_name',
                    help = 'The lipid name(s). If multiple, separate them with _')
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
    args = parser.parse_args()

    return args

def log(y_list):
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

def fit_decay(x,y, LIMX, LIMY):
    """
    Function does a linear fit.
    
    The linear fit is made on the probability of having a certain packing area.
    
    --------------------
    INPUT
    x : numpy array
        Packing area.
    y : numpy array
        probability of having a certain packing area.
    LIMX : int
        The lowest defect area used for the fit.
    LIMY : int
        The lowest probability used for the fit
    
    --------------------
    OUTPUT
    numpy array
        A linear fit of x and y.
    """
    # fit with defects above LIMX nm and proba > LIMY
    y = y[x >= LIMX]
    x = x[x >= LIMX]
    x = x[y >= LIMY]
    y = y[y >= LIMY]
    FIT = np.polyfit(x, log(y), 1)
    return (FIT)

def block_averaging(vect, nb_block):
    """
    Divide the packing data into n blocks.
    
    --------------------
    INPUT
    vect : pandas dataframe
        The size of the packing defect.
    nb_block : int
        The number of blocks we want to have.
    
    --------------------
    OUTPUT
    list
        Avector of 3 decays.
    """
    bornes = [int((len(vect)/3)*nb) for nb in range(nb_block+1)]
    decays = []
    for i in range(0,nb_block):
        subvect = vect[bornes[i]:bornes[i+1]]
        H = plt.hist(subvect, bins = np.arange(0.5, max(vect)+0.5))
        x = H[1]+0.5
        x = x[:len(x)-1]
        y = H[0]/sum(H[0])
        FIT = fit_decay(x, y, args.limx, args.limy)
        decays.append(abs(1/FIT[0]))
    return decays

def plot_defect_fit(type, packdef_constants):
    """
    Compute and plot the fit of the defects distribution.

    --------------------
    INPUT
    type: string
        the type to analyse: Total/Total_Up/Total_Lo
    packdef_constants: pandas DataFrame
        Will contains the statistics of the defects found in this function
    
    --------------------
    OUTPUT
    pandas DataFrame
        Contains the statistics of the defects found in this function
    """
    # Now loop over the three default types
    for defect in ["Deep","Shallow","All"]:
        # load PackMem data
        filename = f"{type}_{args.lipid_name}_{defect}_clean.txt"
        packdef_data = pd.read_csv(filename, header=None, sep=r"\s+")[1]
        
        # Compute PackDef distributions (on the whole set)
        H = plt.hist(packdef_data, bins=np.arange(0.5,max(packdef_data)+0.5))
        # Length of the defects
        x = np.arange(0, max(packdef_data)-1)
        # nb of observations of a certain defect length
        y = H[0]/sum(H[0])
        FIT = fit_decay(x,y, args.limx, args.limy)
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
        plt.axvline(args.limx, color='gray', linestyle='--')
        plt.axhline(math.log(args.limy), color='gray', linestyle='-')
        plt.plot(x, fit_function(x), color='red', label='Fit')
        pdf.savefig()  # Save the current figure to the PDF
        plt.close()
    
        global_inv_decay = abs(1/FIT[0])
        
        # print global results
        print("")
        print(f"Results on {args.lipid_name} for {defect} {type} defects")
        
        # Compute PackDef distributions with block averaging method, and print results
        FITS_3blocks=block_averaging(packdef_data, 3)
        for i in range(3):
            print(f"Using block {i+1}: {round(FITS_3blocks[i],args.precision)} A^2")
        inv_decay_block = sum(FITS_3blocks)/len(FITS_3blocks)
        error_inv_decay_block = np.std(FITS_3blocks)
        print(f"Mean +/- sd on 3 blocks: {round(inv_decay_block,args.precision)} +/- {round(error_inv_decay_block,args.precision)} A^2")

        # store all the results in the data frame
        packdef_constants.loc["PackDef_cst_global", f"{defect}_{type}"] = global_inv_decay
        packdef_constants.loc["PackDef_cst_block1", f"{defect}_{type}"] = FITS_3blocks[0]
        packdef_constants.loc["PackDef_cst_block2", f"{defect}_{type}"] = FITS_3blocks[1]
        packdef_constants.loc["PackDef_cst_block3", f"{defect}_{type}"] = FITS_3blocks[2]
        packdef_constants.loc["PackDef_cst_all_blocks", f"{defect}_{type}"] = inv_decay_block
        packdef_constants.loc["error_all_blocks", f"{defect}_{type}"] = error_inv_decay_block
    
    return packdef_constants

def plot_defect_constants_blocks(type, packdef_constants, errors):
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
    """    
    # Plot and save the second figure (bar plot for just one row)
    packdef_constants[[f"Deep_{type}", f"Shallow_{type}", f"All_{type}"]][:5].T.plot.bar(color=["darkred", "firebrick", "indianred", "lightcoral", "rosybrown"], yerr=errors, capsize=3, rot=0)
    plt.ylabel("Defect size constant ${A^2}$")
    plt.title("Packing defect constants")
    pdf.savefig()  # Save the current figure to the PDF
    plt.close()

def plot_defect_constants(type, packdef_constants, errors):
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
    """
    cst_packing = pd.DataFrame(packdef_constants.loc["PackDef_cst_all_blocks", [f"Deep_{type}", f"Shallow_{type}", f"All_{type}"]])
    df_reset = cst_packing.reset_index()
    
    br = df_reset.plot.bar(x='index', y='PackDef_cst_all_blocks', color=["firebrick", "forestgreen", "royalblue"], rot=0, yerr=errors, capsize=3, legend=False)
    for p in br.patches:
        br.annotate(f'{p.get_height()/max(packdef_constants.loc["PackDef_cst_all_blocks"])*100:.1f}%', (p.get_x() + p.get_width() / 2, p.get_height()),
                    ha='center', va='top', xytext=(0, 10), textcoords='offset points')
        br.annotate(f'{p.get_height():.1f} $A^2$', (p.get_x() + p.get_width() / 2, p.get_height()),
                    ha='center', va='top', xytext=(0, -10), textcoords='offset points')
    plt.ylabel("Defect size constant ${A^2}$")
    plt.title("Packing defect constants computed by block averaging")
    pdf.savefig()  # Save the current figure to the PDF
    plt.close()

if __name__=="__main__":
    # Get arguments
    args = get_arguments()
    # Open a pdf device
    # Create a single PDF file
    pdf = PdfPages(f'{args.output}.pdf')
    # Initialize a data frame to store packdef constants + errors
    packdef_constants = pd.DataFrame(columns = ['Deep_Total', 'Shallow_Total', 'All_Total', 'Deep_Total_Up', 'Shallow_Total_Up', 'All_Total_Up', 'Deep_Total_Lo', 'Shallow_Total_Lo', 'All_Total_Lo'], 
                       index = ["PackDef_cst_global", "PackDef_cst_block1", "PackDef_cst_block2", 
                                "PackDef_cst_block3","PackDef_cst_all_blocks","error_all_blocks"])

    for type in ["Total", "Total_Up", "Total_Lo"]:
        # Plot the fit of the defects distribution
        packdef_constants = plot_defect_fit(type, packdef_constants)
    
        # Plot all packdef constants on a single barplot
        # Allows to estimate the relative convergence of the simulation
        errors=packdef_constants.T.error_all_blocks.to_frame('PackDef_cst_all_blocks')
        plot_defect_constants_blocks(type, packdef_constants, errors)
        
        # Plot the final decays + errors computed with block averaging
        # (for each packdef) on a barplot
        plot_defect_constants(type, packdef_constants, errors)

    
    # Close the PDF file
    pdf.close()
