#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 11:18:29 2024

@author: magdalena
"""

import pandas as pd
import numpy as np
import os
from scipy.stats import ttest_ind_from_stats
import module_andalena as mym

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# resources
resources_1to8 = [
    'starch','glycerol','glucose','fructose',
    'citrate','butyrate','ascorbic','acetate'  #es el orden en el que aparecen en el dataset
]

#resources_8to1 = resources_1to8[::-1]

for st in ['PA','KT','P1','P2','P3','Serratia','Salmonella']:

    path_save_results = f'../results/{st}'
    os.makedirs(path_save_results, exist_ok=True)

    fil = open(f'{path_save_results}/fitness_effects.csv','w')

    fil.write(
        'Resource,Background environment,' +
        ','.join(resources_1to8) +
        ',Fitness background mean,Fitness background sd,' +
        'Fitness b+ri mean,Fitness b+ri sd,' +
        'Fitness effect mean,Fitness effect sd,p-value\n'
    )

    df = pd.read_csv(
        f'../data/3_clean_data_module/{st}_L8_resources_module.csv'
    )
    resource_cols = resources_1to8
    df['Resources'] = df[resource_cols].astype(str).agg(''.join, axis=1)

    # identify replicate columns dynamically
    rep_cols = [c for c in df.columns if 'rep' in c.lower()]

    # function to get mean/sd per environment
    def get_stats(env):
        subset = df[df['Resources'] == env]
        return subset['mean'].values[0], subset['sd'].values[0]

    every_environment = mym.possible_environments(8)

    for i, ri in enumerate(resources_1to8):

        envs_with_ri = [env for env in every_environment if env[i] == '1']

        for env in envs_with_ri:

            env_list = list(env)

            env_background_list = env_list[:]
            env_background_list[i] = '0'

            env_background = ''.join(env_background_list)
            env_ri = ''.join(env_list)

            # skip missing data
            if env_background not in df['Resources'].values:
                continue
            if env_ri not in df['Resources'].values:
                continue

            # fitness values from mean/sd columns
            f_background, f_background_std = get_stats(env_background)
            f_with_ri, f_with_ri_std = get_stats(env_ri)

            # effect
            f_effect = f_with_ri - f_background

            # approximate propagation of error
            f_effect_std = np.sqrt(f_background_std**2 + f_with_ri_std**2)

            # t-test (using stats)
            t_stat, p_value = ttest_ind_from_stats(
                mean1=f_background,
                std1=f_background_std,
                nobs1=len(rep_cols),
                mean2=f_with_ri,
                std2=f_with_ri_std,
                nobs2=len(rep_cols)
            )

            # write output
            fil.write(
                ri + ',' + env_background + ',' +
                ','.join(env_background_list) + ',' +
                str(f_background) + ',' + str(f_background_std) + ',' +
                str(f_with_ri) + ',' + str(f_with_ri_std) + ',' +
                str(f_effect) + ',' + str(f_effect_std) + ',' +
                str(p_value) + '\n'
            )

    fil.close()