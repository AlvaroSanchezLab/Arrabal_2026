import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import os
import module_andalena as mym
from scipy.stats import ttest_ind_from_stats



# =========================
# WORKDIR
# =========================
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# =========================
# RESOURCES
# =========================
resources = [
    'starch',
    'glycerol',
    'glucose',
    'fructose',
    'citrate',
    'butyrate',
    'ascorbic',
    'acetate'
]

# =========================
# STRAINS
# =========================
strains = ['PA', 'KT', 'P1', 'P2', 'P3', 'Salmonella', 'Serratia']

# =========================
# LOOP STRAINS
# =========================
for st in strains:

    print(f'\nProcessing {st}')

    # =========================
    # PATHS
    # =========================
    path_save = f'../results/{st}'
    path_pairwise = f'{path_save}/pairwise_interaction'

    os.makedirs(path_pairwise, exist_ok=True)

    # =========================
    # LOAD DATA
    # =========================
    df = pd.read_csv(
        f'../data/3_clean_data_module/{st}_L8_resources_module.csv'
    )

    # ensure binary columns are integers
    df[resources] = df[resources].astype(int)

    # =========================
    # OUTPUT FILE
    # =========================
    fil = open(f'{path_save}/interaction_pairwise.csv', 'w')

    fil.write(
        'Resource_i,Resource_j,Background_size,Background,Interaction,pvalue\n'
    )

    # =========================
    # STORAGE
    # =========================
    e_positive = []
    e_negative = []
    pairs = []

    lowest_interaction = 0
    lowest_pair = ''

    # =========================
    # LOOP RESOURCE PAIRS
    # =========================
    for i in range(len(resources)):

        ri = resources[i]

        for j in range(i + 1, len(resources)):

            rj = resources[j]

            pair_name = f'{ri}-{rj}'
            pairs.append(pair_name)

            print(f'  Pair: {pair_name}')

            eps_list = []
            bg_sizes = []
            pvalues = []

            # =========================
            # SELECT ENVIRONMENTS
            # only rows where both ri and rj are present
            # =========================
            df_pair = df[
                (df[ri] == 1) &
                (df[rj] == 1)
            ]

            # =========================
            # LOOP ENVIRONMENTS
            # =========================
            for idx, row in df_pair.iterrows():

                # ---------------------------------
                # BACKGROUND
                # ---------------------------------
                other_resources = [
                    r for r in resources
                    if r not in [ri, rj]
                ]

                bg_vals = row[other_resources]

                bg_size = int(bg_vals.sum())

                # ---------------------------------
                # SAME BACKGROUND MASK
                # ---------------------------------
                mask_bg = (
                    df[other_resources] == bg_vals.values
                ).all(axis=1)

                # ---------------------------------
                # B
                # ---------------------------------
                mask_0 = (
                    mask_bg &
                    (df[ri] == 0) &
                    (df[rj] == 0)
                )

                # ---------------------------------
                # Bi
                # ---------------------------------
                mask_i = (
                    mask_bg &
                    (df[ri] == 1) &
                    (df[rj] == 0)
                )

                # ---------------------------------
                # Bj
                # ---------------------------------
                mask_j = (
                    mask_bg &
                    (df[ri] == 0) &
                    (df[rj] == 1)
                )

                # ---------------------------------
                # Bij
                # current row
                # ---------------------------------
                mask_ij = (
                    mask_bg &
                    (df[ri] == 1) &
                    (df[rj] == 1)
                )

                # ---------------------------------
                # CHECK ALL EXIST
                # ---------------------------------
                if (
                    mask_0.sum() == 0 or
                    mask_i.sum() == 0 or
                    mask_j.sum() == 0 or
                    mask_ij.sum() == 0
                ):
                    continue

                # ---------------------------------
                # EXTRACT VALUES
                # ---------------------------------
                g_0 = df.loc[mask_0, 'mean'].values[0]
                g_i = df.loc[mask_i, 'mean'].values[0]
                g_j = df.loc[mask_j, 'mean'].values[0]
                g_obs = df.loc[mask_ij, 'mean'].values[0]

                # standard deviations
                sd_0 = df.loc[mask_0, 'sd'].values[0]
                sd_i = df.loc[mask_i, 'sd'].values[0]
                sd_j = df.loc[mask_j, 'sd'].values[0]
                sd_obs = df.loc[mask_ij, 'sd'].values[0]

                # ---------------------------------
                # ADDITIVE EXPECTATION
                # epsilon = gij - gi - gj + g0
                # ---------------------------------
                g_exp = g_i + g_j - g_0

                eps_val = g_obs - g_exp

                # expected variance propagation
                sd_exp = np.sqrt(
                    sd_i**2 +
                    sd_j**2 +
                    sd_0**2
                )

                # ---------------------------------
                # T-TEST
                # observed vs expected
                # ---------------------------------
                t_stat, p_val = ttest_ind_from_stats(
                    mean1=g_obs,
                    std1=sd_obs,
                    nobs1=3,
                    mean2=g_exp,
                    std2=sd_exp,
                    nobs2=3
                )

                # ---------------------------------
                # STORE
                # ---------------------------------
                eps_list.append(eps_val)
                bg_sizes.append(bg_size)
                pvalues.append(p_val)

                # background string
                bg_string = ''.join(
                    [str(int(v)) for v in bg_vals.values]
                )

                fil.write(
                    f'{ri},{rj},{bg_size},{bg_string},{eps_val},{p_val}\n'
                )

            # =========================
            # SKIP EMPTY
            # =========================
            if len(eps_list) == 0:
                e_positive.append(np.nan)
                e_negative.append(np.nan)
                continue

            # =========================
            # POSITIVE / NEGATIVE
            # =========================
            eps_array = np.array(eps_list)

            pos_frac = np.mean(eps_array > 0)
            neg_frac = np.mean(eps_array < 0)

            e_positive.append(pos_frac)
            e_negative.append(neg_frac)

            # =========================
            # MOST NEGATIVE
            # =========================
            if np.min(eps_array) < lowest_interaction:

                lowest_interaction = np.min(eps_array)
                lowest_pair = pair_name

            # =========================
            # COLORS
            # =========================
            colors = []

            for e, p in zip(eps_array, pvalues):

                if e > 0 and p <= 0.05:
                    colors.append('#743C8A')

                elif e < 0 and p <= 0.05:
                    colors.append('#B5933C')

                else:
                    colors.append('#987d82')

            # =========================
            # PLOT
            # =========================
            fig, ax = plt.subplots(figsize=(5, 4))

            ax.scatter(
                bg_sizes,
                eps_array,
                c=colors,
                alpha=0.8
            )

            ax.axhline(
                0,
                linestyle=':',
                color='grey'
            )

            for s in ['top', 'right']:
                ax.spines[s].set_visible(False)

            ax.set_xlabel('Resources in background')
            ax.set_ylabel('Observed - expected')

            plt.tight_layout()
            plt.show()

            fig.savefig(
                f'{path_pairwise}/interaction_{ri}_{rj}.pdf',
                bbox_inches='tight'
            )

            plt.close()

    # =========================
    # CLOSE FILE
    # =========================
    fil.close()

    # =========================
    # PRINT
    # =========================
    print('\nMost negative interaction:')
    print(lowest_pair, lowest_interaction)

    # =========================
    # SUMMARY BARPLOT
    # =========================
    fig, ax = plt.subplots(figsize=(8, 3))

    x = np.arange(len(e_positive))

    ax.bar(
        x,
        e_negative,
        color='#B5933C'
    )

    ax.bar(
        x,
        e_positive,
        bottom=e_negative,
        color='#743C8A'
    )

    ax.set_xticks(x)

    ax.set_xticklabels(
        pairs,
        rotation=90,
        fontsize=8
    )

    for s in ['top', 'right', 'left', 'bottom']:
        ax.spines[s].set_visible(False)

    ax.grid(
        visible=True,
        linestyle='-.',
        linewidth=0.5,
        alpha=0.2
    )

    ax.set_ylabel('Fraction interactions')

    plt.tight_layout()
    plt.show()

    fig.savefig(
        f'{path_save}/pairwise_pos_neg_interaction.pdf',
        bbox_inches='tight'
    )

    plt.close()