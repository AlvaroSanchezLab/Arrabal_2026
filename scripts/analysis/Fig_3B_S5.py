
"""
Created on Fri May 15 13:24:46 2026

@author: andrea
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

# =========================
# WORKDIR
# =========================
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# =========================
# RESOURCES
# =========================
resources = [
    'acetate',
    'ascorbic',
    'butyrate',
    'citrate',
    'fructose',
    'glucose',
    'glycerol',
    'starch'
]

# =========================
# COLORS
# =========================
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

# =========================
# FUNCTION
# =========================
def plot_global_epistasis(
    ax_l,
    ax_r,
    df_fe,
    df_ei,
    df_FEEs_theory,
    resource,
    min_background,
    max_background,
    min_fitness_effect,
    max_fitness_effect,
    min_ei,
    max_ei
):

    # =========================
    # FITNESS EFFECTS
    # =========================
    sub = df_fe[df_fe['Resource'] == resource]

    f_background = sub['Fitness background mean'].values
    f_background_sd = sub['Fitness background sd'].values

    f_effect = sub['Fitness effect mean'].values
    f_effect_sd = sub['Fitness effect sd'].values

    col = resource_colors[resource]

    # =========================
    # LEFT PANEL
    # =========================
    ax_l.errorbar(
        f_background,
        f_effect,
        xerr=f_background_sd,
        yerr=f_effect_sd,
        fmt='o',
        markerfacecolor=col,
        markeredgecolor='black',
        markeredgewidth=0.01,
        markersize=12,
        ecolor='#D3D3D3',
        elinewidth=0.7,
        capsize=1.5,
        alpha=0.9
    )

    # =========================
    # LINEAR FIT
    # =========================
    slope, intercept = np.polyfit(f_background, f_effect, 1)

    y_pred = slope * f_background + intercept

    r2 = r2_score(f_effect, y_pred)

    x_vals = np.array([min_background, max_background])

    ax_l.plot(
        x_vals,
        x_vals * slope + intercept,
        color='black',
        linewidth=1.5
    )
    
    # =========================
    # LINE THEORY
    # =========================
    slope_theory = df_FEEs_theory.loc[df_FEEs_theory['Resource '] == resource, 'slope (b)']
    intercept_theory = df_FEEs_theory.loc[df_FEEs_theory['Resource '] == resource, 'intercept (a)']

    ax_l.plot(
        x_vals,
        [x * slope_theory + intercept_theory for x in x_vals],
        linestyle='--',
        color='grey',
        linewidth=1.5
    )
    
    # =========================
    # ZERO LINE
    # =========================
    ax_l.axhline(
        0,
        linestyle=':',
        color='grey',
        linewidth=1
    )

    # =========================
    # LIMITS
    # =========================
    ax_l.set_xlim(min_background, max_background)
    ax_l.set_ylim(min_fitness_effect, max_fitness_effect)

    # =========================
    # TITLE
    # =========================
    ax_l.set_title(
        resource,
        fontsize=12,
        color=col,
        fontweight='bold'
    )

    # =========================
    # R²
    # =========================
    ax_l.text(
        0.05,
        0.90,
        r'$r^2=$' + str(round(r2, 2)),
        transform=ax_l.transAxes,
        fontsize=18,
        color='black'
    )

    # aesthetics
    for s in ['top', 'right']:
        ax_l.spines[s].set_visible(False)

    ax_l.spines['left'].set_linewidth(2)
    ax_l.spines['bottom'].set_linewidth(2)

    ax_l.tick_params(labelsize=18)

    # =========================
    # RIGHT PANEL:
    # EFFECTIVE INTERACTIONS
    # =========================
    sub_ei = df_ei[df_ei['Resource i'] == resource]

    effective_interaction = sub_ei['Effective interaction']
    names = sub_ei['Resource j']

    color_bar = []

    for ei in effective_interaction:

        if ei > 0:
            color_bar.append('#348bd4ff') 

        elif ei < 0:
            color_bar.append('#cc0b0eff')

        else:
            color_bar.append('#987d82ff')

    ax_r.barh(
        names,
        effective_interaction,
        color=color_bar,
        alpha=0.4
    )

    ax_r.axvline(
        0,
        linestyle='--',
        linewidth=1,
        color='grey'
    )

    ax_r.set_xlim(min_ei, max_ei)

    ax_r.set_xlabel(
        'Effective interaction',
        fontsize=16
    )

    # aesthetics
    for s in ['top', 'bottom', 'left', 'right']:
        ax_r.spines[s].set_visible(False)

    ax_r.tick_params(labelsize=18)

    ax_r.grid(
        visible=True,
        color='grey',
        linestyle='-.',
        linewidth=0.5,
        alpha=0.2
    )

    return slope, intercept, r2


# =========================
# STRAINS
# =========================
for st in [
    'PA',
    'KT',
    #'P1',
    'P2',
    'P3',
    'Serratia',
    'Salmonella'
]:

    print(f'\nProcessing: {st}')

    path_save_results = f'../results/{st}'
    path_save_results_ge = f'{path_save_results}/global_epistasis'

    os.makedirs(path_save_results_ge, exist_ok=True)

    # =========================
    # LOAD DATA
    # =========================
    df_fe = pd.read_csv(
        f'../results/{st}/fitness_effects.csv',
        dtype={'Background environment': object}
    )

    df_ei = pd.read_csv(
        f'../results/{st}/effective_interaction.csv'
    )
    
    df_FEEs_theory = pd.read_csv(
        f'../results/{st}/FEEs_slope_intercept_theory.csv'
    )

    # =========================
    # GLOBAL LIMITS
    # =========================
    min_fitness_effect = (
        df_fe['Fitness effect mean']
        - df_fe['Fitness effect sd']
    ).min()

    max_fitness_effect = (
        df_fe['Fitness effect mean']
        + df_fe['Fitness effect sd']
    ).max()

    min_background = (
        df_fe['Fitness background mean']
        - df_fe['Fitness background sd']
    ).min()

    max_background = (
        df_fe['Fitness background mean']
        + df_fe['Fitness background sd']
    ).max()

    # =========================
    # EI LIMITS
    # =========================
    min_ei = df_ei['Effective interaction'].min()
    max_ei = df_ei['Effective interaction'].max()

    # =========================
    # OUTPUT FILE
    # =========================
    fil = open(f'{path_save_results}/FEEs.csv', 'w')
    fil.write('Resource,Intercept,Slope,r2\n')

    # =====================================================
    # INDIVIDUAL PLOTS
    # =====================================================
    for resource in resources:

        fig, (ax_l, ax_r) = plt.subplots(
            1,
            2,
            figsize=(8, 4),
            gridspec_kw={'width_ratios': [4, 1]}
        )

        slope, intercept, r2 = plot_global_epistasis(
            ax_l,
            ax_r,
            df_fe,
            df_ei,
            df_FEEs_theory,
            resource,
            min_background,
            max_background,
            min_fitness_effect,
            max_fitness_effect,
            min_ei,
            max_ei
        )

        # =========================
        # DEBUG
        # =========================
        sub = df_fe[df_fe['Resource'] == resource]

        f_background = sub['Fitness background mean'].values
        f_effect = sub['Fitness effect mean'].values

        print("\n", resource)
        print("corr Pearson:", np.corrcoef(f_background, f_effect)[0, 1])
        print("R² sklearn:", r2)

        # =========================
        # LABELS
        # =========================
        ax_l.set_xlabel(
            'Fitness in background environment',
            fontsize=18
        )

        ax_l.set_ylabel(
            'Fitness effect',
            fontsize=18
        )

        # =========================
        # SAVE RESULTS
        # =========================
        fil.write(
            f'{resource},{intercept},{slope},{r2}\n'
        )

        plt.tight_layout()

        fig.savefig(
            f'{path_save_results_ge}/{resource}_global_epistasis_plot.pdf',
            bbox_inches='tight'
        )
        plt.show()
        plt.close()

    fil.close()

