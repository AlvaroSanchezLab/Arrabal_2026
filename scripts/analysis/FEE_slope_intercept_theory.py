
import pandas as pd
from matplotlib import pyplot as plt
import os

# Run from the scripts/ directory regardless of where this is launched from,
# so the relative ../results and ../data paths resolve correctly.
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

base_path = '../results'


strains = ['P1', 'P2', 'P3', 'KT', 'PA', 'Salmonella', 'Serratia']

resources = ['starch','glycerol','glucose','fructose','citrate','butyrate','ascorbic','acetate']

resource_colors = {
    'acetate':  '#264653',
    'ascorbic': '#2a9d8f',
    'butyrate': '#8AB17D',
    'citrate':  '#E9C46A',
    'fructose': '#F4A261',
    'glucose':  '#E36040',
    'glycerol': '#BC6B85',
    'starch':   '#9576C9'
}

for st in strains:

    print(f"Processing {st}...")
    
    # =========================
    # OUTPUT FILE
    # =========================
    
    path_save_results = f'../results/{st}'
    
    fil = open(f'{path_save_results}/FEEs_slope_intercept_theory.csv', 'w')
    fil.write('Resource ,slope (b),intercept (a)\n')
    
    # =========================
    # OPEN DATA
    # =========================

    df_ei = pd.read_csv(f'{base_path}/{st}/effective_interaction.csv')
    df_fe = pd.read_csv(f'{base_path}/{st}/fitness_effects.csv')

    # =========================
    # CALCULATE SLOPE AND INTERCEPT
    # =========================
    for ri in resources:

        # slope theory
        eff = df_ei.loc[df_ei['Resource i'] == ri, 'Effective interaction']
        slope_theory = eff.sum()
        
        # intercept theory
        mean_fitness_effect_r = df_fe.loc[df_fe['Resource'] == ri, 'Fitness effect mean'].mean()
        mean_fitness_background = df_fe.loc[df_fe['Resource'] == ri, 'Fitness background mean'].mean()
        intercept_theory = mean_fitness_effect_r - slope_theory*mean_fitness_background
        
        # write to file
        fil.write(ri + ',' + str(slope_theory) + ',' + str(intercept_theory) + '\n')
        
    fil.close()