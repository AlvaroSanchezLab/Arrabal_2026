# Code and data for *[Statistical learning of bacterial growth in combinatorially constructed environments]*

This repository contains the analysis code that reproduces the figures of the
paper *[Statistical learning of bacterial growth in combinatorially constructed environments, Andrea Arrabal, Magdalena San Román, Juan Diaz-Colunga, Alvaro Sanchez, preprint 2026]*.

The **data** are archived separately on Zenodo:
**[10.5281/zenodo.20848258]**

To run the pipeline, download the `data/` folder from Zenodo and place it next
to `scripts/` so the layout is:

```
.
├── data/        # downloaded from Zenodo
│   ├── 0_raw_data/                  # raw OD600 plate readouts (.xlsx), one folder per strain
│   ├── 1_csv_clean_data/            # cleaned CSVs (output of step 1)
│   ├── 2_csv_clean_batch_effect/    # batch-corrected CSVs (output of step 2)
│   ├── 3_clean_data_module/         # final per-strain landscapes used by all analyses
│   ├── growthcurve_glucose_gradients/   # glucose-gradient growth curves (Fig. S1)
│   ├── d_plateWell_resourceCombinationBinary.json
│   └── d_resourceCombinationBinary_plateWell.json
├── scripts/
│   ├── clean_data/   # data-cleaning pipeline (steps 1–3)
│   └── analysis/     # figure-generating scripts
└── results/          # created automatically when scripts run (git-ignored)
```

All scripts locate their inputs **relative to their own location**, so they can
be run from any working directory without editing paths.

## Requirements

- Python 3.10+ — `pip install -r requirements.txt`
- R (for the single `.R` cleaning step) with: `tidyverse`, `readxl`, `readr`,
  `tidyr`, `ggplot2`, `ggh4x`, `scales`, `reshape2`, `colorspace`, `stringr`,
  `grid`.
- **Epistasia** (Bayesian batch-effect correction), maintained by
  Camacho-Mateu. Required only for `clean_data/2_clean_batch_effect.py` and
  `analysis/Fig_3C_S7.py`. Install it separately and point the `EPISTASIA_DATA`
  environment variable at its data folder.

## Strains

`KT`, `Salmonella`, `Serratia`, `P1`, `P2`, `P3`, `PA`.

## How to reproduce

### A. Data cleaning (`scripts/clean_data/`)

Only needed to regenerate the cleaned data from raw; the cleaned outputs are
already provided on Zenodo.

1. **`1_clean_data.py`** — processes the raw OD600 `.xlsx` files in
   `0_raw_data/`, extracts blanks, corrects measurements and writes the cleaned
   datasets to `1_csv_clean_data/`. Requires `module_andalena.py` and the two
   `d_*.json` files.
2. **`2_clean_batch_effect.py`** — Bayesian batch-effect correction with the
   **Epistasia** package; produces `2_csv_clean_batch_effect/`.
3. **`3_to_generate_clean_data_module.R`** — reformats column headers/structure
   (no biological/numerical values are changed) and writes the final
   `3_clean_data_module/` files used by every downstream analysis.

(Supplementary Fig. 1 uses the separate `growthcurve_glucose_gradients/` data
via `analysis/Fig_S1.py`.)

### B. Analysis (`scripts/analysis/`)

> `module_andalena.py` contains shared helper functions and must be importable
> (it lives in the same folder, so running the scripts as-is is enough).

Some scripts consume the CSVs that earlier scripts write into `results/`, so run
them in this order:

| Order | Script | Produces / Figure |
|------|--------|-------------------|
| 1 | `Fig_1A_S2.py` | Fig. 1A, S2 |
| 2 | `Fig_1C_S3.py` | Fig. 1C, S3 |
| 3 | `fitness_effect_calculation_every_strain.py` | `results/<strain>/fitness_effects.csv` |
| 4 | `Fig_2.py` | Fig. 2 |
| 5 | `Fig_S4.py` | Fig. S4 |
| 6 | `interaction_pairwise_csv.py` | `results/<strain>/interaction_pairwise.csv` |
| 7 | `effective_interactions_csv.py` | `results/<strain>/effective_interaction.csv` |
| 8 | `FEE_slope_intercept_theory.py` | `results/<strain>/FEEs_slope_intercept_theory.csv` |
| 9 | `Fig_3B_S5.py` | Fig. 3B, S5 (writes `FEEs.csv`) |
| 10 | `Fig_3C_S7.py` | Fig. 3C, S7 (requires Epistasia) |
| 11 | `Fig_3D_S6.py` | Fig. 3D, S6 |
| — | `Fig_S1.py` | Fig. S1 (independent; glucose gradients) |

Outputs are written to `results/` (created automatically, git-ignored).

## License

- **Code**: MIT (see `LICENSE`).
- **Data** (Zenodo): CC-BY-4.0.

## Citation

If you use this code or data, please cite the paper and the Zenodo archive.
