import argparse
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from matplotlib.backends.backend_pdf import PdfPages

# Getting the arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', action = 'store', dest = 'filename',
                help = 'The first clean .txt file to analyse')
parser.add_argument('-f2', action = 'store', dest = 'filename2',
                help = 'The second clean .txt file to analyse')
parser.add_argument('-p', action = 'store', dest = 'precision',
                help = 'The precision for writing packdef cosntants (nb of decimals) in the output. Default = 2', type=int, default=2)
parser.add_argument('-lx', action = 'store', dest = 'limx',
                help = 'The lowest defect area used for the fit (we recommand not to touch to this value). Default = 15', type=int, default=15)
parser.add_argument('-ly', action = 'store', dest = 'limy',
                help = 'The lowest probability used for the fit (we recommand not to touch to this value). Default = 1e-4', type=float, default=1e-4)
parser.add_argument('-o', action = 'store', dest = 'output',
                help = 'The name of the output .pdf file. Default = Res_membrane', default="Res_membrane")
args = parser.parse_args()

################################################################################

def log(y_list):
    """Calculate the logarithm of a list.
    
    Parameters
    ----------
    y_list : list
        The y values.
    
    Returns
    -------
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
    """Function does a linear fit.
    
    The linear fit is made on the probability of having a certain packing area.
    
    Parameters
    ----------
    x : numpy array
        Packing area.
    y : numpy array
        probability of having a certain packing area.
    LIMX : int
        The lowest defect area used for the fit.
    LIMY : int
        The lowest probability used for the fit
    Returns
    -------
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
    """Divide the packing data into n blocks.
    
    Parameters
    ----------
    vect : pandas dataframe
        The size of the packing defect.
    nb_block : int
        The number of blocks we want to have.
    
    Returns
    -------
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


if __name__=="__main__":
    # Open a pdf device
    # Create a single PDF file
    pdf = PdfPages(f'{args.output}.pdf')
    # Initialize a data frame to store packdef constants + errors
    packdef_constants1 = pd.DataFrame(columns = ['all'], 
                       index = ["PackDef_cst_global", "PackDef_cst_block1", "PackDef_cst_block2", 
                                "PackDef_cst_block3","PackDef_cst_all_blocks","error_all_blocks"])
    packdef_constants2 = pd.DataFrame(columns = ['all'], 
                       index = ["PackDef_cst_global", "PackDef_cst_block1", "PackDef_cst_block2", 
                                "PackDef_cst_block3","PackDef_cst_all_blocks","error_all_blocks"])
    
    # load PackMem data
    filename1 = args.filename
    packdef_data1 = pd.read_csv(filename1, header=None, delim_whitespace=True)[1]
    filename2 = args.filename2
    packdef_data2 = pd.read_csv(filename2, header=None, delim_whitespace=True)[1]
    
    # Compute PackDef distributions (on the whole set)
    H = plt.hist(packdef_data1, bins=np.arange(0.5,max(packdef_data1)+0.5))
    # Length of the defects
    x1 = np.arange(0, max(packdef_data1)-1)
    # nb of observations of a certain defect length
    y1 = H[0]/sum(H[0])
    FIT1 = fit_decay(x1,y1, args.limx, args.limy)
    fit_function1 = np.poly1d(FIT1)

    # Compute PackDef distributions (on the whole set)
    H = plt.hist(packdef_data2, bins=np.arange(0.5,max(packdef_data2)+0.5))
    # Length of the defects
    x2 = np.arange(0, max(packdef_data2)-1)
    # nb of observations of a certain defect length
    y2 = H[0]/sum(H[0])
    FIT2 = fit_decay(x2,y2, args.limx, args.limy)
    fit_function2 = np.poly1d(FIT2)

    ### Plot for defect fit
    plt.clf()
    plt.scatter(x1, log(y1), marker='o', facecolor="none", edgecolor="orange")
    plt.xlim(-2, 102)
    plt.ylim(-10, -3.5)
    plt.yticks([math.log(1e-4), math.log(1e-3), math.log(1e-2)], labels=[str(1e-4), str(1e-3), str(1e-2)])
    plt.ylabel("Probability")
    plt.xlabel("Defect area ${A^2}$")
    plt.title("all")
    plt.axvline(args.limx, color='gray', linestyle='--')
    plt.axhline(math.log(args.limy), color='gray', linestyle='-')
    plt.plot(x1, fit_function1(x1), color='orange', label='Fit')
    plt.scatter(x2, log(y2), marker='o', facecolor="none", edgecolor="seagreen")
    plt.plot(x2, fit_function2(x2), color='seagreen', label='Fit')
    pdf.savefig()  # Save the current figure to the PDF
    plt.close()

    global_inv_decay1 = abs(1/FIT1[0])
    global_inv_decay2 = abs(1/FIT2[0])
    
    # print global results
    print("")
    print(f"Results for all defects")
    print(f"Total number of defects = {len(packdef_data1)}")
    print(f"Using all data: {round(global_inv_decay1, args.precision)} A^2")
    
    # Compute PackDef distributions with block averaging method, and print results
    FITS_3blocks1=block_averaging(packdef_data1, 3)
    for i in range(3):
        print(f"Using block {i+1}: {round(FITS_3blocks1[i],args.precision)} A^2")
    inv_decay_block1 = sum(FITS_3blocks1)/len(FITS_3blocks1)
    error_inv_decay_block1 = np.std(FITS_3blocks1)
    print(f"Mean +/- sd on 3 blocks: {round(inv_decay_block1,args.precision)} +/- {round(error_inv_decay_block1,args.precision)} A^2")

    FITS_3blocks2=block_averaging(packdef_data2, 3)
    for i in range(3):
        print(f"Using block {i+1}: {round(FITS_3blocks2[i],args.precision)} A^2")
    inv_decay_block2 = sum(FITS_3blocks2)/len(FITS_3blocks2)
    error_inv_decay_block2 = np.std(FITS_3blocks2)
    print(f"Mean +/- sd on 3 blocks: {round(inv_decay_block2,args.precision)} +/- {round(error_inv_decay_block2,args.precision)} A^2")

    # store all the results in the data frame
    packdef_constants1["all"]["PackDef_cst_global"] = global_inv_decay1
    packdef_constants1["all"]["PackDef_cst_block1"] = FITS_3blocks1[0]
    packdef_constants1["all"]["PackDef_cst_block2"] = FITS_3blocks1[1]
    packdef_constants1["all"]["PackDef_cst_block3"] = FITS_3blocks1[2]
    packdef_constants1["all"]["PackDef_cst_all_blocks"] = inv_decay_block1
    packdef_constants1["all"]["error_all_blocks"] = error_inv_decay_block1

    packdef_constants2["all"]["PackDef_cst_global"] = global_inv_decay2
    packdef_constants2["all"]["PackDef_cst_block1"] = FITS_3blocks2[0]
    packdef_constants2["all"]["PackDef_cst_block2"] = FITS_3blocks2[1]
    packdef_constants2["all"]["PackDef_cst_block3"] = FITS_3blocks2[2]
    packdef_constants2["all"]["PackDef_cst_all_blocks"] = inv_decay_block2
    packdef_constants2["all"]["error_all_blocks"] = error_inv_decay_block2

    errors1=packdef_constants1.T.error_all_blocks.to_frame('PackDef_cst_all_blocks')
    errors2=packdef_constants2.T.error_all_blocks.to_frame('PackDef_cst_all_blocks')
    
    
    ####
    # FINAL PLOT
    # Now we plot the final decays + errors computed with block averaging
    # (for each packdef) on a barplot
    cst_packing1 = pd.DataFrame(packdef_constants1.loc["PackDef_cst_all_blocks"])
    df_reset1 = cst_packing1.reset_index()
    cst_packing2 = pd.DataFrame(packdef_constants2.loc["PackDef_cst_all_blocks"])
    df_reset2 = cst_packing2.reset_index()
    
    br1 = df_reset1.plot.bar(x='index', y='PackDef_cst_all_blocks', color=["orange", "seagreen"], rot=0, yerr=errors1, capsize=3, legend=False)
    br2 = df_reset2.plot.bar(x='index', y='PackDef_cst_all_blocks', color=["orange", "seagreen"], rot=0, yerr=errors2, capsize=3, legend=False)
    for p in br1.patches:
        br1.annotate(f'{p.get_height()/max(packdef_constants1.loc["PackDef_cst_all_blocks"])*100:.1f}%', (p.get_x() + p.get_width() / 2, p.get_height()),
                    ha='center', va='top', xytext=(0, 10), textcoords='offset points')
        br1.annotate(f'{p.get_height():.1f} $A^2$', (p.get_x() + p.get_width() / 2, p.get_height()),
                    ha='center', va='top', xytext=(0, -10), textcoords='offset points')
    for p in br2.patches:
        br2.annotate(f'{p.get_height()/max(packdef_constants2.loc["PackDef_cst_all_blocks"])*100:.1f}%', (p.get_x() + p.get_width() / 2, p.get_height()),
                    ha='center', va='top', xytext=(0, 10), textcoords='offset points')
        br2.annotate(f'{p.get_height():.1f} $A^2$', (p.get_x() + p.get_width() / 2, p.get_height()),
                    ha='center', va='top', xytext=(0, -10), textcoords='offset points')
    plt.ylabel("Defect size constant ${A^2}$")
    plt.title("Packing defect constants computed by block averaging")
    pdf.savefig()  # Save the current figure to the PDF
    plt.close()
    
    # Close the PDF file
    pdf.close()