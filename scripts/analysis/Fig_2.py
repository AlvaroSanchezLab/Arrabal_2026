#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 11:18:02 2026

@author: andrea
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

resources = [
    'starch','glycerol','glucose','fructose',
    'citrate','butyrate','ascorbic','acetate'
]

# color gris solicitado
hist_color = "#c2c2c2"

for st in ['PA','KT','P1','P2','P3','Salmonella','Serratia']:

    df = pd.read_csv(f'../results/{st}/fitness_effects.csv')

    # nueva carpeta
    path_out = f'../results/{st}/histograms_grey'
    os.makedirs(path_out, exist_ok=True)

    # =========================
    # GLOBAL LIMITS
    # =========================
    all_effects = df['Fitness effect mean']

    min_fe = all_effects.min()
    max_fe = all_effects.max()

    # bins comunes
    bins = np.linspace(min_fe, max_fe, 40)

    # mismo Y máximo para todos
    y_max_global = 0
    hist_data = {}

    # calcular histogramas
    for r in resources:

        sub = df[df['Resource'] == r]['Fitness effect mean']

        counts, _ = np.histogram(sub, bins=bins)

        hist_data[r] = sub

        y_max_global = max(y_max_global, counts.max())

    # =========================
    # PLOT POR RECURSO
    # =========================
    for r in resources:

        data = hist_data[r]

        fig, ax = plt.subplots(figsize=(5,4))

        ax.hist(
            data,
            bins=bins,
            color=hist_color,
            edgecolor="white"
        )

        ax.axvline(
            0,
            linestyle='--',
            color='grey',
            linewidth=1
        )

        # mismas escalas
        ax.set_xlim(min_fe, max_fe)
        ax.set_ylim(0, y_max_global * 1.1)

        ax.set_xlabel("Fitness effect")
        ax.set_ylabel("")
        ax.set_title(f"{st} - {r}")

        ax.grid(False)

        plt.tight_layout()

        plt.savefig(
            f'{path_out}/hist_{r}_grey.pdf',
            bbox_inches='tight'
        )

        plt.close()