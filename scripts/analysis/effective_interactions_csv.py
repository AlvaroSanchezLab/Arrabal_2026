#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 13:14:22 2026

@author: andrea
"""

import pandas as pd
import numpy as np
import os

# =========================
# WORKING DIRECTORY
# =========================
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# =========================
# RESOURCES (NEW FORMAT CONSISTENCY)
# =========================
resources_1to8 = [
    'acetate','ascorbic','butyrate','citrate',
    'fructose','glucose','glycerol','starch'
]
resources_8to1 = resources_1to8[::-1]

# =========================
# STRAINS
# =========================
strains = ['PA','KT','P1','P2','P3', 'Salmonella', 'Serratia']

# =========================
# LOOP
# =========================
for st in strains:

    print(f"Processing: {st}")

    path_save_results = f'../results/{st}'

    # =========================
    # LOAD DATA (NEW FORMAT)
    # =========================
    df_fe = pd.read_csv(
        f'../results/{st}/fitness_effects.csv',
        dtype={'Background environment': object}
    )

    df_e = pd.read_csv(
        f'../results/{st}/interaction_pairwise.csv',
        dtype={'Background environment': object}
    )

    # =========================
    # OUTPUT FILE
    # =========================
    fil = open(f'{path_save_results}/effective_interaction.csv', 'w')
    fil.write('Resource i,Resource j,Effective interaction\n')

    # =========================
    # COMPUTE EFFECTIVE INTERACTION
    # =========================
    for ri in resources_8to1:

        for rj in resources_8to1:

            if ri == rj:
                continue

            # -------------------------
            # mean pairwise interaction
            # -------------------------
            mask_ij = (
                ((df_e['Resource_i'] == ri) & (df_e['Resource_j'] == rj)) |
                ((df_e['Resource_i'] == rj) & (df_e['Resource_j'] == ri))
            )

            mean_interaction_ri_rj = df_e.loc[mask_ij, 'Interaction'].mean()

            if np.isnan(mean_interaction_ri_rj):
                continue

            # -------------------------
            # fitness effect of rj when ri absent
            # -------------------------
            mean_fitness_effect_rj = df_fe.loc[
                (df_fe['Resource'] == rj) & (df_fe[ri] == 0),
                'Fitness effect mean'
            ].mean()

            # -------------------------
            # normalization term
            # -------------------------
            norm_terms = []

            for rk in resources_8to1:
                if rk == ri:
                    continue

                val = df_fe.loc[
                    (df_fe['Resource'] == rk) & (df_fe[ri] == 0),
                    'Fitness effect mean'
                ].mean()

                if not np.isnan(val):
                    norm_terms.append(val ** 2)

            denom = np.sum(norm_terms)

            if denom == 0 or np.isnan(denom):
                effective_interaction = np.nan
            else:
                effective_interaction = (
                    mean_interaction_ri_rj *
                    mean_fitness_effect_rj /
                    denom
                )

            # -------------------------
            # SAVE
            # -------------------------
            fil.write(f"{ri},{rj},{effective_interaction}\n")

    fil.close()