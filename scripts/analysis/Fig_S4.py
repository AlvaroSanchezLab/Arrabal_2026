import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import random

# =========================================
# WORKDIR
# =========================================
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# =========================================
# STRAINS
# =========================================
strains = ['P1','P2','P3','PA','KT','Salmonella','Serratia']



resources_ordered = ['starch', 'glycerol', 'glucose', 'fructose', 'citrate', 'butyrate', 'ascorbic', 'acetate']

# =========================================
# GLOBAL X LIMITS (FITNESS EFFECT)
# =========================================
all_vals = []

for st in strains:

    file_path = f'../results/{st}/fitness_effects.csv'

    if not os.path.exists(file_path):
        continue

    df_tmp = pd.read_csv(file_path,
                          dtype={'Background environment': object})

    all_vals.extend(df_tmp['Fitness effect mean'].values)

x_min = np.min(all_vals)
x_max = np.max(all_vals)

pad = 0.05 * (x_max - x_min)
x_min -= pad
x_max += pad

# =========================================
# COUNT POSITIVE / NEGATIVE FITNESS EFFECTS
# =========================================

for st in strains:

    file_path = f'../results/{st}/fitness_effects.csv'

    if not os.path.exists(file_path):
        continue

    df = pd.read_csv(file_path,
                     dtype={'Background environment': object})

    print(f"\n{'='*50}")
    print(f"STRAIN: {st}")
    print(f"{'='*50}")

    for resource in resources_ordered:

        vals = df.loc[df['Resource'] == resource,
                      'Fitness effect mean']

        n_pos = (vals > 0).sum()
        n_neg = (vals < 0).sum()
        n_zero = (np.abs(vals) < 1e-9).sum()

        print(
            f"{resource:10s} | "
            f">0: {n_pos:3d} | "
            f"<0: {n_neg:3d} | "
            f"=0: {n_zero:3d}"
        )

# =========================================
# FIGURE
# =========================================
fig, axes = plt.subplots(1, 7, figsize=(16, 4), sharey=True)
axes = axes.flatten()

# =========================================
# PLOT FUNCTION
# =========================================
def plot_strain(ax, st):

    file_path = f'../results/{st}/fitness_effects.csv'

    if not os.path.exists(file_path):
        ax.set_title(st, fontsize=11)
        ax.axis("off")
        return

    df = pd.read_csv(file_path,
                     dtype={'Background environment': object})

    ax.axvline(x=0, linestyle='--', linewidth=1, color='grey')

    for ir, resource in enumerate(resources_ordered):

        vals = df.loc[df['Resource'] == resource]['Fitness effect mean']

        if len(vals) == 0:
            continue

        fe_pos = vals[vals > 0]
        fe_neg = vals[vals < 0]
        fe_neu = vals[np.abs(vals) < 1e-9]

        y_pos = ir + np.random.uniform(-0.25, 0.25, len(fe_pos))
        y_neg = ir + np.random.uniform(-0.25, 0.25, len(fe_neg))
        y_neu = ir + np.random.uniform(-0.25, 0.25, len(fe_neu))

        ax.scatter(fe_pos, y_pos, c="#348BD4", alpha=0.4, s=25)
        ax.scatter(fe_neg, y_neg, c="#CC0B0E", alpha=0.4, s=25)
        ax.scatter(fe_neu, y_neu, c="grey", alpha=0.4, s=25)

        ax.errorbar(
            np.mean(vals),
            ir,
            xerr=np.std(vals),
            fmt='o',
            color='black',
            ecolor='black',
            alpha=0.6,
            capsize=3,
            markersize=4
        )

    ax.set_title(st, fontsize=11)

    ax.set_yticks(range(len(resources_ordered)))
    ax.set_yticklabels(resources_ordered, fontsize=10)

    ax.set_xlim(x_min, x_max)

    ax.tick_params(axis='x', labelsize=10)

    for s in ['top','right','left','bottom']:
        ax.spines[s].set_visible(False)

    ax.grid(True, linestyle='--', alpha=0.2)


# =========================================
# LOOP
# =========================================
for ax, st in zip(axes, strains):
    plot_strain(ax, st)

axes[0].set_xlabel('Fitness effect')

plt.tight_layout()

# =========================================
# SAVE
# =========================================
output = '../results/fitness_effect_stripplot.pdf'
fig.savefig(output, bbox_inches='tight', dpi=300)
plt.close(fig)

print(f"Saved: {output}")