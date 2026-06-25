#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 13:56:07 2026

@author: andrea
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.ticker as mticker
import pandas as pd
import seaborn as sns 
import scipy.stats as st 
import os

# set working directory
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# results path
path_results_script = '../results/growthcurve_glucose_gradients/'

plt.rcParams["font.family"] = "Arial"

# =========================
# confidence interval
# =========================
def confidence_interval_95(data, confidence=0.95):
    data = np.array(data)
    m = np.mean(data)
    sem = st.sem(data)
    h = sem * st.t.ppf((1 + confidence) / 2., len(data)-1)
    return m, m - h, m + h


# glucose concentrations used
glc_concentration = [1.6]
for i in range(15):
    glc_concentration.append(glc_concentration[-1]/2)

# convert to C-moles (as in your script)
glc_moles = [c*0.00666/0.02 for c in glc_concentration]

# palette
color_palette = ['#322929','#b18971','#906e60','#6e6a64','#8d94b3','#73879c','#b6c5d5']

# strains info
strains_name = ['PA','KT','P1','P2','P3', 'Salmonella', 'Serratia']
starting_ix = [11,9,9,9,9,9,9]

i_blank = [23,24,35,36,47,48,59,60,71,72,83,84]

for si, sp in enumerate(strains_name):

    # load data
    df = pd.read_excel('../data/growthcurve_glucose_gradients/' + sp + '_gradient.xlsx')
    
    # time in hours
    time_hs = df.iloc[starting_ix[si]:starting_ix[si]+144,1] / (60*60)
    
    i_sp = 15
    
    od24_mean = []
    od24_std = []
    fitness_mean = []
    fitness_std = []
    
    f = plt.figure()
    ax = plt.subplot(111)
    
    for ic, conc in enumerate(glc_concentration):

        # blank
        df_blank = df.iloc[starting_ix[si]:starting_ix[si]+144, i_blank].mean(axis=1)

        # replicates
        data_here = df.iloc[starting_ix[si]:starting_ix[si]+144, [i_sp, i_sp+24, i_sp+48]]

        # update index
        if conc == 0.0125:
            i_sp = 27
        else:
            i_sp += 1

        df_strain_mean = data_here.mean(axis=1) - df_blank
        df_strain_std = data_here.std(axis=1)

        od24_mean.append(df_strain_mean.iloc[-1])
        od24_std.append(df_strain_std.iloc[-1])

        fitness = []
        for i in range(len(data_here.columns)):
            final_od = data_here.iloc[-1, i] - df_blank.iloc[-1]
            initial_od = data_here.iloc[0, i] - df_blank.iloc[0]
            fitness_here = np.log(final_od / initial_od)
            fitness.append(fitness_here)

        fitness_mean.append(np.mean(fitness))
        fitness_std.append(np.std(fitness))

        # growth curves
        plt.errorbar(time_hs, df_strain_mean, df_strain_std, label=str(conc)+'%')

        plt.ylabel('Fitness (F)', size=12, family='Arial')
        plt.xlabel('Time (h)', size=12, family='Arial')
        plt.xticks(size=12, family='Arial')
        plt.yticks(size=12, family='Arial')

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()
    f.savefig(path_results_script + "growthcurve_lower_glc_" + sp + ".pdf", bbox_inches='tight')


    # =========================
    # OD24 vs glucose %
    # =========================
    f = plt.figure()
    plt.errorbar(glc_concentration, od24_mean, yerr=od24_std, fmt='ko')
    plt.ylabel('Fitness (F)', size=12)
    plt.xlabel('Glucose concentration (%)', size=12)
    plt.title(sp)
    plt.show()
    f.savefig(path_results_script + "OD24_glc_" + sp + ".pdf", bbox_inches='tight')

    f = plt.figure()
    ax = plt.axes()
    ax.set_xscale("log")

    plt.errorbar(glc_moles, od24_mean, yerr=od24_std, fmt='ko')

    # reference lines with legend
    ax.axvline(x=0.003, linestyle='--', color='#76AD1F', linewidth=1.5, label='0.003 C-mol/L')
    ax.axvline(x=0.024, linestyle='--', color='#175717', linewidth=1.5, label='0.024 C-mol/L')

    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.4f'))

    plt.ylabel('Fitness (F)', size=12)
    plt.xlabel('C-mol/L', size=12)
    plt.title(sp)

    ax.legend()

    plt.show()
    f.savefig(path_results_script + "OD24_molesCperL_loglog_" + sp + ".pdf", bbox_inches='tight')


    # =========================
    # fitness vs C-moles
    # =========================
    f = plt.figure()
    ax = plt.axes()
    ax.set_xscale("log")

    plt.errorbar(glc_moles, fitness_mean, yerr=fitness_std, fmt='ko')

    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.4f'))

    plt.ylabel('Log(final OD/initial OD)', size=12)
    plt.xlabel('C-mol/L per L', size=12)
    plt.title(sp)

    plt.show()
    f.savefig(path_results_script + "fitness_molesCperL_loglog_" + sp + ".pdf", bbox_inches='tight')