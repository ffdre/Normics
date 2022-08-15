### (c) Franz F. Dressler 2022
### Please cite https://doi.org/10.1016/j.mcpro.2022.100269


### import modules
import urllib
import Tkinter, Tkconstants, tkFileDialog
from PIL import ImageTk, Image
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import pylab as plt
import pandas as pd
import multiprocessing as mp
import subprocess
import sys, os, inspect
import scipy.stats
import math


plt.close('all')
directory = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
os.chdir(directory)
colors = ["#0072C5","#DC4438"]

################## prompt input ##################
def Exit():
    root.destroy()
    sys.exit()

def Start():
    root.destroy()

def name_file():
    global data_file
    data_file = tkFileDialog.askopenfilename(initialdir=directory,
                                             title="Select input data file",
                                             filetypes=(("excel files", "*.xlsx"),("all files", "*.*"))) #
    button_text = "Selected input data file: ..." + data_file[-20:]
    Tkinter.Label(root, text=button_text).grid(row=3, column=1, columnspan=2)
    pass

root = Tkinter.Tk()
root.title("Normics v.3.5.10.1")
#Image
logo = ImageTk.PhotoImage(Image.open("S2_NORMICS.jpg"))
Tkinter.Label(root, image=logo).grid(
    row=0, column=0, columnspan=4)


Tkinter.Label(root, text="").grid(
    row=1, column=1, columnspan=2, sticky="NSEW")
Tkinter.Button(root, text="Select input data file", highlightbackground="#0072C5", command=name_file).grid(
    row=2, column=1, columnspan=2, sticky="NSEW")

min_data_values = Tkinter.StringVar(value = 0.8)
Tkinter.Label(root, text='Minimum share of samples per protein').grid(
    row=4, column=1, columnspan=1, sticky="W")
Tkinter.Entry(root, textvariable=min_data_values).grid(
    row=4, column=2, columnspan=1, sticky="NSEW")

reduce_correlation_by = Tkinter.StringVar(value = 1)
Tkinter.Label(root, text='Reduce correlation by [1 = No reduction]').grid(
    row=5, column=1, columnspan=1, sticky="W")
Tkinter.Entry(root, textvariable=reduce_correlation_by).grid(
    row=5, column=2, columnspan=1, sticky="NSEW")

cores_to_spare = Tkinter.StringVar(value = 1)
Tkinter.Label(root, text='# of cores not to use for computation').grid(
    row=6, column=1, columnspan=1, sticky="W")
Tkinter.Entry(root, textvariable=cores_to_spare).grid(
    row=6, column=2, columnspan=1, sticky="NSEW")

start_from_longlist = Tkinter.BooleanVar(value = False)
Tkinter.Label(root, text='Re-run with existing computation').grid(
    row=7, column=1, columnspan=1, sticky="W")
Tkinter.Checkbutton(root,variable=start_from_longlist, onvalue=True, offvalue=False).grid(
    row=7, column=2, columnspan=1, sticky="NSEW")

TMT_ratio = Tkinter.BooleanVar(value = False)
Tkinter.Label(root, text='Ratio-reported data (e.g. TMT)').grid(
    row=8, column=1, columnspan=1, sticky="W")
Tkinter.Checkbutton(root,variable=TMT_ratio, onvalue=True, offvalue=False).grid(
    row=8, column=2, columnspan=1, sticky="NSEW")


Tkinter.Label(root, text='\n____________________________________ Additional settings ____________________________________').grid(
    row=10, column=1, columnspan=2, sticky="NSEW")

handle_zero_as_nan = Tkinter.BooleanVar(value = True)
Tkinter.Label(root, text='Handle zero intensities as None').grid(
    row=11, column=1, columnspan=1, sticky="W")
Tkinter.Checkbutton(root,variable=handle_zero_as_nan, onvalue=True, offvalue=False).grid(
    row=11, column=2, columnspan=1, sticky="NSEW")

reverse_log2 = Tkinter.BooleanVar(value = False)
Tkinter.Label(root, text='Reverse log2').grid(
    row=12, column=1, columnspan=1, sticky="W")
Tkinter.Checkbutton(root,variable=reverse_log2, onvalue=True, offvalue=False).grid(
    row=12, column=2, columnspan=1, sticky="NSEW")

reverse_log = Tkinter.BooleanVar(value = False)
Tkinter.Label(root, text='Reverse log').grid(
    row=13, column=1, columnspan=1, sticky="W")
Tkinter.Checkbutton(root,variable=reverse_log, onvalue=True, offvalue=False).grid(
    row=13, column=2, columnspan=1, sticky="NSEW")

reduce_data_by = Tkinter.StringVar(value = 1)
Tkinter.Label(root, text='Reduction of input data [1 = No reduction]').grid(
    row=14, column=1, columnspan=1, sticky="W")
Tkinter.Entry(root, textvariable=reduce_data_by).grid(
    row=14, column=2, columnspan=1, sticky="NSEW")

VSN_quantile = Tkinter.StringVar(value = 0.5)
Tkinter.Label(root, text='VSN_quantile').grid(
    row=15, column=1, columnspan=1, sticky="W")
Tkinter.Entry(root, textvariable=VSN_quantile).grid(
    row=15, column=2, columnspan=1, sticky="NSEW")

NormicsVSN_quantile = Tkinter.StringVar(value = 0.8)
Tkinter.Label(root, text='VSN_quantile for Normics subset').grid(
    row=16, column=1, columnspan=1, sticky="W")
Tkinter.Entry(root, textvariable=NormicsVSN_quantile).grid(
    row=16, column=2, columnspan=1, sticky="NSEW")


Tkinter.Label(root, text="").grid(
    row=25, column=1, columnspan=2, sticky="NSEW")
Tkinter.Button(root, text="Start", highlightbackground="#DC4438", command=Start).grid(
    row=26, column=1,columnspan=2, sticky="NSEW")
Tkinter.Button(root, text="Cancel", command=Exit).grid(
    row=28, column=1, columnspan=2, sticky="NSEW")
Tkinter.Label(root, text="").grid(
    row=29, column=1, columnspan=2, sticky="NSEW")

root.mainloop()

##### basic settings #####
# copy this file to the same directory as a folder containing your data
# set the data_folder name
data_folder = "/".join(data_file.split("/")[:-1])
# set the data file name (samples as columns; proteins as rows; row 0 := sample names; column 0 := protein names)
data_file = data_file
# set the share of minimum sample values per protein to be included. Default = 0.8
min_data_values = float(min_data_values.get())
# if proteins exceed 2000, consider reducing data by integer steps (1, 2, 3, ...). Default = 1
reduce_correlation_by = int(reduce_correlation_by.get())
# select how many cores should be left for other tasks. Default = 1
cores_to_spare = int(cores_to_spare.get())
# already computed a longlist?
start_from_longlist = start_from_longlist.get()
##########################

### advanced settings ####
handle_zero_as_nan = handle_zero_as_nan.get()
reverse_log2 = reverse_log2.get()
reverse_log = reverse_log.get()
reduce_data_by = int(reduce_data_by.get())
VSN_quantile = float(VSN_quantile.get())
NormicsVSN_quantile = float(NormicsVSN_quantile.get())
TMT_ratio = TMT_ratio.get()
##########################

### create results and figure folders
try: os.mkdir(data_folder + "/NORMICS_results")
except: pass
try: os.mkdir(data_folder + "/NORMICS_figures")
except: pass
try: os.mkdir(data_folder + "/NORMICS_R-files")
except: pass

### save log file
with open(data_folder + "/NORMICS_results/_log-file.txt", 'w') as log_file:
    log_file.write("".join(["_______ Selected parameters: _______\n",
                   "Data folder\t%s\n" %(data_folder),
                    "Data_file\t%s\n" %(data_file),
                    "Min. data values\t%s\n" %(min_data_values),
                    "Correlation reduction\t%s\n" %(reduce_correlation_by),
                    "Cores to spare\t%s\n" %(cores_to_spare),
                    "Ratio-reported\t%s\n" %(TMT_ratio),
                    "Re-run only\t%s\n" %(start_from_longlist),
                    "Reverse log2\t%s\n" %(reverse_log2),
                    "Reverse log\t%s\n" %(reverse_log),
                    "Data reduction\t%s\n" %(reduce_data_by),
                    "VSN_quantile\t%s\n" %(VSN_quantile),
                    "VSN_quantile Normics\t%s\n" %(NormicsVSN_quantile),
                   "___________________________________"]))
print(
"\n_______ Selected parameters: _______")
print("Data folder\t%s" %data_folder)
print("Data_file\t%s" %data_file)
print("Min. data values\t%s" %min_data_values)
print("Correlation reduction\t%s" %reduce_correlation_by)
print("Cores to spare\t%s" %cores_to_spare)
print("Ratio-reported\t%s" %TMT_ratio)
print("Re-run only\t%s" %start_from_longlist)
print("")
print("Reverse log2\t%s" %reverse_log2)
print("Reverse log\t%s" %reverse_log)
print("Data reduction\t%s" %reduce_data_by)
print("VSN_quantile\t%s" %VSN_quantile)
print("VSN_quantile Normics\t%s" %NormicsVSN_quantile)
print(
"___________________________________\n")

### load data
print(
">>> Start ...")

data_all = pd.read_excel(data_file).iloc[::reduce_data_by,:]
num_original = len(data_all.iloc[:,0])


### handle missing values
print(
"... removing all proteins with less than %s %% values" %np.round(min_data_values * 100, 0))
if handle_zero_as_nan == True:
    data_all = data_all.replace([0],[np.nan])
else:
    data_min = np.array(data_all.min(skipna=True, numeric_only=True))[0]
    data_all = data_all.replace([0], [10**-6])
    print("... zeros replaced by %s"%10**-6)
data_full = data_all.copy()
for i, d in data_all.iterrows():
    if sum(list(d.isnull().astype(int))) > int((1 - min_data_values) * (len(data_all.columns) - 1)):
        data_all = data_all.drop([i])
print(
"... %s proteins removed (-%s %%)"%((num_original - len(data_all.iloc[:,0])), np.round((1.- len(data_all.iloc[:,0])/
                                                                                        float(num_original)) * 100, 2)))

### rearrange data
data_proteins = data_all.iloc[:,0]
data_proteins_full = data_full.iloc[:,0]
data = data_all.iloc[:,1:]
data_full = data_full.iloc[:,1:]
if reverse_log2 == True:
    data = 2. ** data
    data_full = 2. ** data_full
elif reverse_log == True:
    data = np.exp(data)
    data_full = np.exp(data_full)


### check data to be non log-transformed
if (data < 0).any().any() == True:
	print("... WARNING: data likely log transformed => set reverse_log2/log to True and restart")
	sys.exit()

data_proteins.to_excel(data_folder + "/NORMICS_results/RES_data_proteins.xlsx", index=False)


### reduce data if chosen
if reduce_correlation_by > 1:
    data_complete = data.copy()
    data_proteins_complete = data_proteins.copy()
    data = data.iloc[::reduce_correlation_by, :]
    data_proteins = data_proteins.iloc[::reduce_correlation_by]
data_corr = data.applymap(np.log2)
data_corr.to_excel(data_folder + "/NORMICS_results/RES_data_corrected.xlsx")


### start calculation or import longlist
cols = ["Protein_ID", "Rank_sum(RS)", "Mean_correlation(MC)", "Coefficient_of_variation(CV)", "Variance_correlation(VC)"]
if start_from_longlist != True:
    longlist = pd.DataFrame(columns=cols)
    print(
    "... start parallel correlation calculation")

    pairs = range(len(data_proteins))

    def calculate_protein(i):
        global data, data_proteins, pairs
        global TMT_ratio
        temp = []
        for j, P in enumerate(data_proteins):
            if i != j:
                A, B = [], []
                for a, b, a_nan, b_nan in zip(data.iloc[i, :], data.iloc[j, :], data.iloc[i, :].isnull(),
                                              data.iloc[j, :].isnull()):
                    if a_nan == False and b_nan == False:
                        A.append(a)
                        B.append(b)
                temp.append(scipy.stats.spearmanr(A, B)[0])
        if TMT_ratio == False:
            CV = data.iloc[i, :].var() ** .5 / data.iloc[i, :].mean()
        elif TMT_ratio == True:
            CV = data.iloc[i, :].var() ** .5
        if np.mod(i, 10) == 0:
            percent_done = np.round(i / float(len(pairs)) * 100, 2)
            sys.stdout.write("... |%s%s|    %s %% done                           \r" % (
            "#" * int(percent_done / 5.), "." * (20 - int(percent_done / 5.)), percent_done))
            sys.stdout.flush()
        return [data_proteins.iloc[i], None, np.nanmean(temp), CV, np.nanvar(temp)]

    pool = mp.Pool(mp.cpu_count() - cores_to_spare)
    res = pool.map(calculate_protein, pairs, chunksize=1)
    pool.close()
    sys.stdout.write("... |%s%s|    %s %% done                           \r" % (
    "#" * int(100 / 5.), "." * (20 - int(100 / 5.)), 100))
    sys.stdout.flush()
    print("")
    for r in res:
        longlist = longlist.append(pd.DataFrame([r], columns=cols))

    # sort and determine rank sum
    print(
    "... sorting longlist")
    longlist = longlist.sort_values(by=cols[3], ascending=True)
    longlist = longlist.reset_index(drop=True)
    for i in range(len(data_proteins)):
        longlist.loc[i, cols[1]] = i
    longlist = longlist.sort_values(by=cols[2], ascending=False)
    longlist = longlist.reset_index(drop=True)
    for i in range(len(data_proteins)):
        longlist.loc[i, cols[1]] = longlist.loc[i, cols[1]] + i
    longlist = longlist.sort_values(by=cols[1], ascending=True)
    longlist = longlist.reset_index(drop=True)
    longlist.to_excel(data_folder + "/NORMICS_results/RES_longlist.xlsx", index=False)

else:
    longlist = pd.read_excel(data_folder + "/NORMICS_results/RES_longlist.xlsx")
print(
"... sorted results:")
print(longlist)


### visualize longlist
top = 50
final = False
while final == False:
    plt.close()
    fig = plt.figure("Sorting results", [8, 8])
    fig.subplots_adjust(hspace=.3)
    fig.subplots_adjust(wspace=.3)
    # subplot 2
    ax = fig.add_subplot(223)
    ax.set_title("Mean vs. SD")
    tops = list(longlist.loc[:top, cols[0]])
    for i, p in enumerate(data_proteins):
        if p in tops:
            ax.plot(np.log2(np.mean(data.iloc[i, :])), np.log2(np.var(data.iloc[i, :]) ** .5), ".", markeredgewidth=0., alpha=.75,
                    color=colors[0])
        else:
            ax.plot(np.log2(np.mean(data.iloc[i, :])), np.log2(np.var(data.iloc[i, :]) ** .5), ".", markeredgewidth=0., alpha=.25,
                    color=colors[1])
    ax.set_xlabel("Mean [log2]")
    ax.set_ylabel("Standard deviation [log2]")
    ax.grid(False)
    # subplot 1
    ax = fig.add_subplot(211)
    ax.set_title("Mean correlation vs. CV")
    ax.set_xlabel("Mean correlation [1]")
    ax.set_ylabel("Coefficient of variation [log2]")
    for i, p in enumerate(longlist.loc[:top-1, cols[0]]):
        ax.plot(longlist.loc[i, cols[2]], np.log2(longlist.loc[i, cols[3]]), ".", markeredgewidth=0., alpha=.75, color=colors[0])
    for i, p in enumerate(longlist.loc[top:, cols[0]]):
        i += top
        ax.plot(longlist.loc[i, cols[2]], np.log2(longlist.loc[i, cols[3]]), ".", markeredgewidth=0., alpha=.25, color=colors[1])
    ax.plot([], [], ".", markeredgewidth=0., alpha=.75, color=colors[0], label="selected n=%s"%top)
    ax.plot([], [], ".", markeredgewidth=0., alpha=.75, color=colors[1], label="remaining n=%s"%(len(data_proteins)-top))
    ax.legend(fancybox=False)
    ax.grid(False)
    # subplot 3
    ax = fig.add_subplot(224)
    ax.set_title("Sample coverage")
    tops = list(longlist.loc[:top, cols[0]])
    heatmap = []
    for i, p in enumerate(data_proteins):
        if p in tops:
            heatmap.append(data.iloc[i, :].isna())
    ax.imshow(np.array(heatmap), interpolation='none', aspect='auto', cmap=plt.get_cmap("Blues_r"))
    ax.set_xlabel("Samples [1]")
    ax.set_ylabel("Selected proteins [1]")
    ax.grid(False)
    plt.show()
    fig.savefig(data_folder + "/NORMICS_figures/PLOT_data-longlisted.png", dpi = 500)
    # select
    top = int(raw_input("... choose top X proteins for normalization: (integer) "))
    final = raw_input("... visualize again? (y/n) ")
    if final == "y":
        final = False
    else:
        final = True
shortlist = longlist.iloc[:top, :]


### calculate NORMICSmedian ratios
params_ratios = pd.DataFrame(columns = data.columns)
data_proteins_annotated = data_proteins.copy()
for i, protein in enumerate(data_proteins):
    if protein in list(shortlist.loc[:, cols[0]]):
        data_proteins_annotated.update(pd.Series([str(protein) + "_norm"], index=[i]))
        params_ratios = params_ratios.append(data.iloc[i, :])
params_ratios.to_excel(data_folder + "/NORMICS_results/RES_params_medians-raw.xlsx", index=False)
ratios = pd.DataFrame(columns = data.columns)
temp = []
for j, r in enumerate(params_ratios.columns):
    temp.append(params_ratios.iloc[:, j].median())
ratios = ratios.append(pd.DataFrame([temp], columns = data.columns))
ratios.to_excel(data_folder + "/NORMICS_results/RES_params_medians-ratios.xlsx", index=False)


### create R files
data_R_shortlist = pd.DataFrame(columns = data.columns)
for i, protein in enumerate(data_proteins):
    if protein in list(shortlist.loc[:, cols[0]]):
        data_R_shortlist = data_R_shortlist.append(data.iloc[i, :])

data_R_shortlist.to_csv(data_folder + "/NORMICS_R-files/data_R_shortlist.csv", sep=';', index=False)
if reduce_correlation_by > 1:
    data = data_complete
    data_proteins = data_proteins_complete

data.to_csv(data_folder + "/NORMICS_R-files/data_R.csv", sep=';', index=False)
with open(data_folder + "/NORMICS_R-files/normalizeVSN.R", "w") as F:
    F.write('library(readxl)\n'
            'library(limma)\n'
            'library(openxlsx)\n'
            'library(NormalyzerDE)\n'
            'library(vsn)\n')
    F.write('dir <- "%s"\n' % (data_folder + "/NORMICS_R-files"))
    F.write('data <- as.matrix(read.table(file.path(dir,"data_R_shortlist.csv"), header=TRUE, sep=";", as.is=TRUE))\n'

            'data_all <- as.matrix(read.table(file.path(dir,"data_R.csv"), header=TRUE, sep=";", as.is=TRUE))\n'

            'data_median <- medianNormalization(data_all)\n'  # requires raw data
            'df_m = data.frame(data_median)\n'
            'write.xlsx(df_m,file.path(dir,"RES_data_normalizedMEDIAN_log2.xlsx"))\n'

            'data_quantile <- normalizeQuantiles(log2(data_all))\n'
            'df_q = data.frame(data_quantile)\n'
            'write.xlsx(df_q,file.path(dir,"RES_data_normalizedQUANTILE_log2.xlsx"))\n'

            'data_loess <- normalizeCyclicLoess(log2(data_all))\n'
            'df_l = data.frame(data_loess)\n'
            'write.xlsx(df_l,file.path(dir,"RES_data_normalizedLOESS_log2.xlsx"))\n'

            'data_all <- ExpressionSet(assayData = data_all)\n'

            'png(file=file.path(dir,"PLOT_before.png"),width=1000, height=1000)\n'
            'meanSdPlot(data)\n')
    F.write('data_new <- normalizeVSN(data, lts.quantile = %s)\n' % NormicsVSN_quantile)
    F.write('png(file=file.path(dir,"PLOT_after.png"),width=1000, height=1000)\n'
            'meanSdPlot(data_new)\n'

            'data <- ExpressionSet(assayData = data)\n')
    F.write('fit <- vsn2(data, lts.quantile = %s)\n' % NormicsVSN_quantile)
    F.write('params <- coef(fit)[1,,]\n')
    F.write('fit_all <- vsn2(data_all, lts.quantile = %s)\n' % VSN_quantile)
    F.write('params_VSN <- coef(fit_all)[1,,]\n'
            'data_all_new <- predict(fit_all, data_all, log2scale = TRUE)\n'
            'df_p = data.frame(params)\n'
            'df_VSN = data.frame(params_VSN)\n'
            'df = data.frame(data_all_new)\n'
            'write.xlsx(df,file.path(dir,"RES_data_VSN-normalized.xlsx"))\n'
            'write.xlsx(df_p,file.path(dir,"RES_parameters_VSN-shortlist.xlsx"))\n'
            'write.xlsx(df_VSN,file.path(dir,"RES_parameters_VSN.xlsx"))')


### run Normics + VSN + comparative normalization algorithms
subprocess.call('Rscript "%s"' % (data_folder + "/NORMICS_R-files/normalizeVSN.R"), shell=True)
# load VSN parameters
params = pd.read_excel(data_folder + "/NORMICS_R-files/RES_parameters_VSN-shortlist.xlsx")
params_VSN = pd.read_excel(data_folder + "/NORMICS_R-files/RES_parameters_VSN.xlsx")
# transform data to glog
def h_transform(Y, ai, bi, B):
    ''' https://www.bioconductor.org/packages/release/bioc/vignettes/vsn/inst/doc/A-vsn.html#sec:calib
    '''
    return np.arcsinh(np.exp(bi) * Y + ai) / np.log(2) - np.log2(2 * np.exp(B))
data_new = data.copy()
data_new_VSN = data.copy()
data_new_full = data_full.copy()
data_new_VSN_full = data_full.copy()
B = np.mean(params.loc[:, "X2"])
B_VSN = np.mean(params_VSN.loc[:, "X2"])
for j, c in enumerate(data.columns):
    ai = params.loc[j, "X1"]
    bi = params.loc[j, "X2"]
    for i, d in enumerate(data.iloc[:, j]):
        data_new.iat[i, j] = h_transform(data.iloc[i, j], ai, bi, B)
    for i, d in enumerate(data_full.iloc[:, j]):
        data_new_full.iat[i, j] = h_transform(data_full.iloc[i, j], ai, bi, B)
    ai_VSN = params_VSN.loc[j, "X1"]
    bi_VSN = params_VSN.loc[j, "X2"]
    for i, d in enumerate(data.iloc[:, j]):
        data_new_VSN.iat[i, j] = h_transform(data.iloc[i, j], ai_VSN, bi_VSN, B_VSN)
    for i, d in enumerate(data_full.iloc[:, j]):
        data_new_VSN_full.iat[i, j] = h_transform(data_full.iloc[i, j], ai_VSN, bi_VSN, B_VSN)
        
# transform data by NORMICSmedians
data_new_ratios = data.copy()
data_new_ratios_full = data_full.copy()
for j, c in enumerate(data.columns):
    ri = ratios.iloc[0, j]
    data_new_ratios.iloc[:, j] = data.iloc[:, j] / float(ri) * ratios.iloc[0, :].mean()
    data_new_ratios_full.iloc[:, j] = data_full.iloc[:, j] / float(ri) * ratios.iloc[0, :].mean()
data_new_ratios = data_new_ratios.applymap(np.log2)
data_new_ratios_full = data_new_ratios_full.applymap(np.log2)

# exclude proteins used for normalization in separate dataset
data_new_noNorm = pd.DataFrame(columns=data.columns)
data_new_med_noNorm = pd.DataFrame(columns=data.columns)
for i, p in enumerate(data_proteins_annotated):
    if "norm" not in str(p):
        data_new_noNorm = data_new_noNorm.append(data_new.iloc[i, :])
        data_new_med_noNorm = data_new_med_noNorm.append(data_new.iloc[i, :])

# transform data
for d in ["QUANTILE", "MEDIAN", "LOESS"]:
    D = pd.read_excel(data_folder + "/NORMICS_R-files/RES_data_normalized%s_log2.xlsx"%d)
    D.columns = list(data.columns)
    D.to_excel(data_folder + "/NORMICS_results/RES_data_normalized%s-log2.xlsx"%d, index=False)
data_new_ratios.to_excel(data_folder + "/NORMICS_results/RES_data_normalizedNORMICSmed-log2.xlsx", index=False)
data_new_ratios_full.to_excel(data_folder + "/NORMICS_results/RES_data_normalizedNORMICSmed-log2_COMPLETE-DATA.xlsx", index=False)
data_new_VSN.to_excel(data_folder + "/NORMICS_results/RES_data_normalizedVSN-glog.xlsx", index=False)
data_new_VSN_full.to_excel(data_folder + "/NORMICS_results/RES_data_normalizedVSN-glog_COMPLETE-DATA.xlsx", index=False)
data_new_noNorm.to_excel(data_folder + "/NORMICS_results/RES_data_normalizedNORMICS-glog_noNorm.xlsx", index=False)
data_new_med_noNorm.to_excel(data_folder + "/NORMICS_results/RES_data_normalizedNORMICSmed-glog_noNorm.xlsx", index=False)
data_new.to_excel(data_folder + "/NORMICS_results/RES_data_normalizedNORMICS-glog.xlsx", index=False)
data_new_full.to_excel(data_folder + "/NORMICS_results/RES_data_normalizedNORMICS-glog_COMPLETE-DATA.xlsx", index=False)
data_proteins.to_excel(data_folder + "/NORMICS_results/RES_data_proteins.xlsx", index=False)
data_proteins_full.to_excel(data_folder + "/NORMICS_results/RES_data_proteins_COMPLETE-DATA.xlsx", index=False)
data_proteins_annotated.to_excel(data_folder + "/NORMICS_results/RES_data_proteins_annotated.xlsx", index=False)
data = data.applymap(np.log2)
data.to_excel(data_folder + "/NORMICS_results/RES_data_non-normalized-log2.xlsx", index=False)

print(
"... data transformed and saved")
# save shortlist
shortlist.to_excel(data_folder + "/NORMICS_results/RES_shortlist.xlsx", index=False)
print(
"... shortlist saved")
print(
">>> ... all done. <<<")