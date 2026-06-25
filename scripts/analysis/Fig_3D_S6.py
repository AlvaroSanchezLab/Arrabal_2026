import pandas as pd
from matplotlib import pyplot as plt
import os

# Run from the scripts/ directory regardless of where this is launched from,
# so the relative ../results and ../data paths resolve correctly.
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

base_path = '../results'

strains = ['P1', 'P2', 'P3', 'KT', 'PA', 'Salmonella', 'Serratia']

resources = ['starch','glycerol','glucose','fructose',
             'citrate','butyrate','ascorbic','acetate']

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

    df_ei = pd.read_csv(f'{base_path}/{st}/effective_interaction.csv')
    df_fees = pd.read_csv(f'{base_path}/{st}/FEEs.csv')

    # Figura cuadrada
    fig, ax = plt.subplots(figsize=(5, 5))

    for ri in resources:

        # Theory
        eff = df_ei.loc[
            df_ei['Resource i'] == ri,
            'Effective interaction'
        ]
        slope_theory = eff.sum()

        # Data
        slope_series = df_fees.loc[
            df_fees['Resource'] == ri,
            'Slope'
        ]

        if slope_series.empty:
            print(f"  Missing {ri} in {st}")
            continue

        slope_data = slope_series.iloc[0]

        ax.plot(
            slope_theory,
            slope_data,
            marker='o',
            linestyle='None',
            markersize=16,
            color=resource_colors.get(ri, 'black'),
            label=ri
        )

    # Obtener límites actuales
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    # Mismos límites en ambos ejes
    lim_min = min(xmin, ymin)
    lim_max = max(xmax, ymax)

    # Añadir un pequeño margen
    pad = 0.05 * (lim_max - lim_min)
    lim_min -= pad
    lim_max += pad

    ax.set_xlim(lim_min, lim_max)
    ax.set_ylim(lim_min, lim_max)

    # Línea y = x
    ax.plot(
        [lim_min, lim_max],
        [lim_min, lim_max],
        color='grey',
        linestyle=':',
        linewidth=1
    )

    # Aspecto 1:1
    ax.set_aspect('equal', adjustable='box')

    # Styling
    ax.legend(
        prop={"size": 9, "family": "Arial"},
        frameon=False
    )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xlabel(
        'FEE slope from theory',
        size=12,
        family='Arial'
    )

    ax.set_ylabel(
        'FEE slope from data',
        size=12,
        family='Arial'
    )

    plt.xticks(size=12, family='Arial')
    plt.yticks(size=12, family='Arial')

    plt.tight_layout()

    # Save inside each strain folder
    outpath = f'{base_path}/{st}/FEE_slope_theory_vs_data.pdf'

    fig.savefig(
        outpath,
        dpi=300,
        bbox_inches='tight'
    )

    plt.close(fig)

    print(f"  Saved: {outpath}")