"""
SCRIPT TO SHOW THE VARIANCE DECOMPOSITION BY INTERACTIONS ORDERS
"""

############################
#    IMPORT LIBRARIES      #
############################

import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator, MultipleLocator
from math import comb

############################
#     CONFIGURE PATHS      #
############################

# Run from the scripts/ directory so relative ../results paths resolve.
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Folder containing the Epistasia data tree (GlobalEpistasis/...).
# Set the EPISTASIA_DATA environment variable, or edit the default below.
base_path = os.environ.get("EPISTASIA_DATA", os.path.expanduser("~/epistasia_data"))

############################
#    IMPORT EPISTASIA      #
############################

# Epistasia must be installed (pip) or available on PYTHONPATH.
# (Package maintained by Camacho-Mateu; see README.)
# If you keep the library inside base_path, uncomment the next line:
# sys.path.insert(1, base_path)
import epistasia as ep

############################################
#          INPUT / OUTPUT DATASETS         #
############################################

# INPUT DIRECTORY
input_dir = os.path.join(base_path, "GlobalEpistasis/Datasets_Clean/Metalandscape_biomass/Resources/")

# LIST ALL DATASETS IN FOLDER
dataset_list = sorted([
    f for f in os.listdir(input_dir)
    if os.path.isfile(os.path.join(input_dir, f)) and not f.startswith(".")
])

# OUTPUT DIRECTORY
output_dir = "../results"
os.makedirs(output_dir, exist_ok=True)

# OUTPUT FILE
output_file = os.path.join(output_dir, "variance_decomposition_by_order.svg")

# SHORT NAMES (mapping of dataset file name -> pretty label).
# Defaults to a names_map.csv next to the input data; override with NAMES_MAP.
names_map_path = os.environ.get("NAMES_MAP", os.path.join(input_dir, "names_map.csv"))
name_df = pd.read_csv(names_map_path)
name_map = dict(zip(name_df["file"], name_df["pretty_name"]))

#############################################
#               DECLARE FIGURE              #
#############################################

N = len(dataset_list)

nrows = 3  
ncols = max(1, math.ceil(N / nrows))

fig, axes = plt.subplots(nrows , ncols, figsize=(9,8), sharex=True, sharey=True)

plt.subplots_adjust(hspace=0.5,wspace=0.3)

#############################################
#             LOOP OVER DATASETS            #
#############################################

for k,dataset in enumerate(dataset_list):
    
    print("Ploting V(k) vs k for "+f"{dataset}")

    # FILE PATH
    file = os.path.join(input_dir,dataset)    

    # LOAD EPISTASIA OBJETC 
    L = ep.landscape_from_file(file)

    ############################
    #  COMPUTE DATA (NO PLOT)  #
    ############################

    fig_tmp, axes_tmp, data = ep.plot_variance_and_amplitude(
        L,
        B_uncertainty=1000,
        B_null=1000,
        as_fraction=True,
        show_components=False,    
        show=False,
        rng=np.random.default_rng(125),
        return_data=True,
        ci_method_uncertainty="bca",
    )

    plt.close(fig_tmp)
    
    #############################
    #        SELECT DATA        #
    #############################
    
    df_variance = data["variance"]
    
    orders = df_variance["Order"].values
    Vk = df_variance["Fraction of variance V(S)/sum_T V(T)"].values
    ci_low = df_variance["CI low"].values
    ci_high = df_variance["CI high"].values
    null_low = df_variance["Null CI low"].values
    null_high = df_variance["Null CI high"].values

    v1 = Vk[orders == 1][0] if np.any(orders == 1) else np.nan
    v2 = Vk[orders == 2][0] if np.any(orders == 2) else np.nan

    yerr = np.vstack([np.clip(Vk - ci_low, 0, None),np.clip(ci_high - Vk, 0, None)])
    
    #############################
    #            PLOT           #
    #############################
    
    i = k // ncols
    j = k % ncols
    ax = axes[i, j] if nrows > 1 and ncols > 1 else axes[max(i, j)]

    ax.fill_between(orders, null_low, null_high, color="gray", alpha=0.25, label="Null CI")
    ax.bar(orders, Vk, color="C0", label="Observed")
    ax.errorbar(orders, Vk, yerr=yerr, fmt="none", ecolor="k", elinewidth=1, capsize=3, zorder=3)
    
    pretty_name = name_map.get(dataset, dataset)
    ax.set_title(f"{pretty_name}", fontsize=10)

    if j == 0:
        ax.set_ylabel(r"$V(k)$", fontsize=10)
    if i == 2:
        ax.set_xlabel("Interaction order k")

    ax.text(
        0.6, 0.6,
        f"V(1) = {100*v1:.1f}%\nV(2) = {100*v2:.1f}%",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,
    )

    if k == 0:
        ax.legend(loc="upper right", frameon=False, fontsize=9)
        
    ax.xaxis.set_major_locator(MultipleLocator(2))

#############################################
#          HIDE EMPTY SUBPLOTS              #
#############################################

for k in range(len(dataset_list), nrows * ncols):
    i = k // ncols
    j = k % ncols
    axes[i, j].axis("off")
    
############################################
#               SAVE FIGURE                #
############################################

plt.tight_layout()
print(f"Saving figure to:\n  {output_file}")
fig.savefig(output_file, bbox_inches="tight", format="svg")
plt.show()