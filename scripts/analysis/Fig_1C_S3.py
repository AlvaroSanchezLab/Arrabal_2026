#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os
import numpy as np

# =========================================
# WORKING DIRECTORY
# =========================================
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# =========================================
# STRAINS
# =========================================
strains = ['P1', 'P2', 'P3', 'PA', 'KT', 'Salmonella', 'Serratia']

# =========================================
# COLORS
# =========================================
light_blue = '#c2c2c2'
dark_blue = '#3d3d3d'

all_strains_data = []

# =========================================
# DATA PROCESSING
# =========================================
for st in strains:
    print('Processing data for:', st)

    path_save_results = f'../results/{st}'
    os.makedirs(path_save_results, exist_ok=True)

    df = pd.read_csv(
        f'../data/3_clean_data_module/{st}_L8_resources_module.csv'
    )

    resource_cols = [
        'starch','glycerol','glucose','fructose',
        'citrate','butyrate','ascorbic','acetate'
    ]

    diversity = df[resource_cols].sum(axis=1)
    function_mean = df['mean']
    function_sd = df['sd']

    total_function_mean = []
    total_function_sd = []

    for div in range(9):
        vals = function_mean[diversity == div]
        total_function_mean.append(np.mean(vals) if len(vals) > 0 else 0)
        total_function_sd.append(np.std(vals) if len(vals) > 0 else 0)

    all_strains_data.append({
        'strain': st,
        'diversity': diversity.values,
        'function_mean': function_mean.values,
        'function_sd': function_sd.values,
        'total_mean': np.array(total_function_mean),
        'total_sd': np.array(total_function_sd)
    })

    # ======================================================
    # INDIVIDUAL FIGURE (FIXED LAYERING)
    # ======================================================
    fig_ind = plt.figure(figsize=(7,5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1], wspace=0.05)

    ax_left = fig_ind.add_subplot(gs[0, 0])

    # --- puntos individuales (BACKGROUND) ---
    ax_left.errorbar(
        diversity,
        function_mean,
        yerr=function_sd,
        fmt='o',
        markersize=3.5,
        markerfacecolor='none',
        markeredgecolor=light_blue,
        ecolor=light_blue,
        linestyle='None',
        alpha=0.5,
        zorder=1,
        capsize=0,
        elinewidth=1
    )

    # --- barras de error media (sin marker) ---
    ax_left.errorbar(
        range(9),
        total_function_mean,
        yerr=total_function_sd,
        fmt='none',
        ecolor=dark_blue,
        elinewidth=1.3,
        capsize=3,
        zorder=2
    )

    # --- puntos media (ENCIMA SIEMPRE) ---
    ax_left.plot(
        range(9),
        total_function_mean,
        'o',
        markersize=5.5,
        markerfacecolor=dark_blue,
        markeredgecolor=dark_blue,
        zorder=10
    )

    ax_left.set_title(st, size=13, weight='bold')
    ax_left.set_xticks(range(9))
    ax_left.set_xlabel('Number of C-sources')
    ax_left.set_ylabel('Fitness (F)')

    for s in ['top','right']:
        ax_left.spines[s].set_visible(False)

    ax_hist = fig_ind.add_subplot(gs[0, 1], sharey=ax_left)

    ax_hist.hist(
        function_mean,
        bins=25,
        orientation='horizontal',
        color=dark_blue,
        edgecolor='white',
        linewidth=0.3,
        alpha=0.85
    )

    for s in ['top','right','bottom']:
        ax_hist.spines[s].set_visible(False)

    ax_hist.tick_params(
        left=False,
        labelleft=False,
        bottom=False,
        labelbottom=False
    )

    # SAVE INDIVIDUAL
    strain_dir = f'../results/{st}'
    os.makedirs(strain_dir, exist_ok=True)

    fig_ind.savefig(
        os.path.join(strain_dir, "DiversityFunction.pdf"),
        bbox_inches='tight'
    )

    plt.close(fig_ind)

# =========================================
# GLOBAL FIGURE
# =========================================
print("\nGenerando figura global...")

fig = plt.figure(figsize=(20,14))

ncols = 3
nrows = 3
outer_gs = gridspec.GridSpec(nrows, ncols, hspace=0.35, wspace=0.28)

all_means = np.concatenate([d['function_mean'] for d in all_strains_data])
y_max = np.max(all_means) * 1.15

for idx, data in enumerate(all_strains_data):

    st_name = data['strain']

    inner_gs = gridspec.GridSpecFromSubplotSpec(
        1, 2,
        subplot_spec=outer_gs[idx],
        width_ratios=[4, 1],
        wspace=0.03
    )

    ax_left = fig.add_subplot(inner_gs[0, 0])

    # --- puntos individuales (BACKGROUND) ---
    ax_left.errorbar(
        data['diversity'],
        data['function_mean'],
        yerr=data['function_sd'],
        fmt='o',
        markersize=3.5,
        markerfacecolor='none',
        markeredgecolor=light_blue,
        ecolor=light_blue,
        linestyle='None',
        alpha=0.5,
        zorder=1,
        capsize=0,
        elinewidth=1
    )

    # --- barras media ---
    ax_left.errorbar(
        range(9),
        data['total_mean'],
        yerr=data['total_sd'],
        fmt='none',
        ecolor=dark_blue,
        elinewidth=1.3,
        capsize=3,
        zorder=2
    )

    # --- puntos media (TOP LAYER) ---
    ax_left.plot(
        range(9),
        data['total_mean'],
        'o',
        markersize=5.5,
        markerfacecolor=dark_blue,
        markeredgecolor=dark_blue,
        zorder=10
    )

    ax_left.set_title(st_name, size=13, weight='bold')
    ax_left.set_xticks(range(9))
    ax_left.set_ylim(0, y_max)

    for s in ['top','right']:
        ax_left.spines[s].set_visible(False)

    if idx >= 4:
        ax_left.set_xlabel('Number of C-sources')

    if idx % 3 == 0:
        ax_left.set_ylabel('Fitness (F)')

    ax_hist = fig.add_subplot(inner_gs[0, 1], sharey=ax_left)

    ax_hist.hist(
        data['function_mean'],
        bins=25,
        orientation='horizontal',
        color=dark_blue,
        edgecolor='white',
        linewidth=0.3,
        alpha=0.85
    )

    for s in ['top','right','bottom']:
        ax_hist.spines[s].set_visible(False)

    ax_hist.tick_params(
        left=False,
        labelleft=False,
        bottom=False,
        labelbottom=False
    )

# limpiar vacíos del grid 3x3
fig.delaxes(fig.add_subplot(outer_gs[2, 1]))
fig.delaxes(fig.add_subplot(outer_gs[2, 2]))

plt.tight_layout()

global_save_path = '../results/S3_FunctionDiversity_3col_GLOBAL.pdf'
os.makedirs('../results', exist_ok=True)

fig.savefig(global_save_path, bbox_inches='tight')
plt.close(fig)

print(f"✔ Figura global guardada en: {global_save_path}")
print(f"✔ Figuras individuales en: ../results/<strain>/DiversityFunction.pdf")