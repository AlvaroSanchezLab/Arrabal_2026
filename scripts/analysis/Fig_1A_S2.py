import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# =========================================
# STRAINS
# =========================================
strains = ['KT', 'P1', 'P2', 'P3', 'Serratia', 'Salmonella', 'PA']

# =========================================
# RESOURCES
# =========================================
resources_order = [
    'starch','glycerol','glucose','fructose',
    'citrate','butyrate','ascorbic','acetate'
]

resource_colors = {
    'acetate':  '#264653',
    'ascorbic': '#2a9d8f',
    'butyrate': '#8AB17D',
    'citrate':  '#E9C46A',
    'fructose': '#F4A261',
    'glucose':  '#E36040',
    'glycerol': '#BC6B85',
    'starch':   '#9576C9',
    'Base medium': 'lightgrey'
}

DETECTION_LIMIT = 0.001

# =========================================
# LOOP
# =========================================
for st in strains:

    print(f"\nProcessing {st}")

    file_path = f'../data/3_clean_data_module/{st}_L8_resources_module.csv'
    if not os.path.exists(file_path):
        print("File not found")
        continue

    df = pd.read_csv(file_path)
    rep_cols = [c for c in df.columns if 'rep' in c]

    # =========================================
    # BASE MEDIUM (00000000)
    # =========================================
    base_mask = df[resources_order].sum(axis=1) == 0

    base_values = df.loc[base_mask, rep_cols].values.flatten()
    base_values = base_values.astype(float)
    base_values = base_values[~np.isnan(base_values)]

    base_mean = np.mean(base_values)
    base_std = np.std(base_values, ddof=1)

    # =========================================
    # SINGLE CARBON SOURCES ONLY
    # =========================================
    df_single = df[df[resources_order].sum(axis=1) == 1]

    results = []

    for r in resources_order:

        sub = df_single[df_single[r] == 1]

        if sub.empty:
            continue

        values = sub[rep_cols].values.flatten()
        values = values.astype(float)
        values = values[~np.isnan(values)]

        if len(values) < 3:
            continue

        mean = np.mean(values)
        std = np.std(values, ddof=1)
        sem = std / np.sqrt(len(values))

        # =========================================
        # DETECTION vs NOISE (CI METHOD)
        # =========================================
        ci_low = mean - 1.96 * sem
        detected = ci_low > DETECTION_LIMIT

        results.append({
            "resource": r,
            "values": values,
            "mean": mean,
            "std": std,
            "detected": detected
        })

    # sort by mean
    results = sorted(results, key=lambda x: x["mean"])

    # =========================================
    # ADD BASE MEDIUM AS FIRST ELEMENT
    # =========================================
    results = [{
        "resource": "Base medium",
        "values": base_values,
        "mean": base_mean,
        "std": base_std,
        "detected": False
    }] + results

    # =========================================
    # PLOT
    # =========================================
    fig, ax = plt.subplots(figsize=(20, 10))

    x_pos = np.arange(len(results))

    for i, r in enumerate(results):

        vals = np.array(r["values"])
        x = np.random.normal(i, 0.05, len(vals))

        ax.scatter(
            x,
            vals,
            s=1000,
            color=resource_colors.get(r["resource"], 'black'),
            alpha=0.85
        )

        ax.errorbar(
            i,
            r["mean"],
            yerr=r["std"],
            fmt='o',
            color='black',
            capsize=6,
            linewidth=2
        )

        

    # =========================================
    # AXES
    # =========================================
    ax.set_xticks(x_pos)
    
    ax.set_xticklabels(
        [r["resource"] for r in results],
        rotation=45,
        ha='right',
        fontsize=30
    )
    
    ax.set_ylabel("Fitness (F)", fontsize=36)
    
    # Zero line
    ax.axhline(
        0,
        color='black',
        linewidth=2,
        zorder=0
    )
    
    # Detection limit
    ax.axhline(
        DETECTION_LIMIT,
        color='black',
        linestyle='--',
        linewidth=2,
        alpha=0.8,
        zorder=0
    )
    

    
    ax.set_xlim(-0.6, len(results) - 0.4)
    
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    
    ax.tick_params(axis='y', labelsize=28)
    ax.tick_params(axis='x', length=0)
    
    plt.tight_layout()

    # =========================================
    # SAVE
    # =========================================
    outdir = f'../results/{st}'
    os.makedirs(outdir, exist_ok=True)

    plt.savefig(
        f"{outdir}/growth_detection_with_base_medium.pdf",
        bbox_inches='tight'
    )

    plt.close()

    print("saved")