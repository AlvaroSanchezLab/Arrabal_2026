"""
SCRIPT FOR AUTOMATING BATCH-EFFECT CORRECTION ACROSS MULTIPLE DATASETS.

This script scans a folder of raw CSV datasets, converts each into an epistasia
Landscape, applies Bayesian batch-effect correction, and saves the cleaned
datasets into a separate output directory (created if it does not exist).

The saved CSV files include the binary state columns (the N features) at the
front, keeping the original feature names stored inside the Landscape object.
"""

############################
#    IMPORT LIBRARIES      #
############################

import os
import sys
import pandas as pd

############################
#     CONFIGURE PATHS      #
############################

# Folder containing the Epistasia data tree (GlobalEpistasis/...).
# Set the EPISTASIA_DATA environment variable, or edit the default below,
# to point to your local copy.
base_path = os.environ.get("EPISTASIA_DATA", os.path.expanduser("~/epistasia_data"))

# Path where the raw datasets are stored
raw_dir = os.path.join(base_path, "GlobalEpistasis", "Datasets_Raw")

# Path where the batch-corrected datasets will be saved
clean_dir = os.path.join(base_path, "GlobalEpistasis", "Datasets_Clean")

# Create the output directory if it does not exist
os.makedirs(clean_dir, exist_ok=True)

##### PROJECTS #####
proyects = ["Sabool/RB2_RB3"]

############################
#    IMPORT EPISTASIA      #
############################

# Epistasia must be installed (pip) or available on PYTHONPATH.
# (Package maintained by Camacho-Mateu; see README.)
# If you keep the library inside base_path, uncomment the next line:
# sys.path.insert(1, base_path)

import epistasia as ep

############################
#    LOOP OVER DATASETS    #
############################

for proyect in proyects:
    
    raw_dir_project = os.path.join(raw_dir, proyect)
    clean_dir_project = os.path.join(clean_dir, proyect)
    
    # Create sub directory for every project
    os.makedirs(clean_dir_project , exist_ok=True)
    
    print(f"\n{'='*40}")
    print(f"ENTERING PROJECT: {proyect}")
    print(f"{'='*40}")

    for fname in os.listdir(raw_dir_project):
        # Process only CSV files (skip other files/folders)
        if not fname.lower().endswith(".csv"):
            continue
    
        print(f"\nProcessing dataset: {fname}")
    
        # Full input and output paths
        data_path = os.path.join(raw_dir_project, fname)
        out_path = os.path.join(clean_dir_project, fname)
    
        ############################
        #      READ RAW DATA       #
        ############################
        try:
            # prefer package helper (handles csv/tsv and simple parsing)
            df = getattr(ep, "read_table", None)
            if callable(df):
                df = ep.read_table(data_path)
            else:
                # fallback to pandas
                df = pd.read_csv(data_path, index_col=False)

            if df is None or df.empty:
                print(f"Skipping empty or unreadable file: {data_path}")
                continue
        except Exception as e:
            print(f"Error reading {data_path}: {e}")
            continue

        # Convert to Landscape object using the core API
        try:
            L_raw = ep.core.Landscape.from_dataframe(df)
        except Exception as e:
            print(f"Error creating Landscape from {data_path}: {e}")
            continue
    
        ############################
        #   CORRECT BATCH EFFECTS  #
        ############################
        L_clean, post = ep.correct_batch_effect(
            L_raw,
            return_posteriors=True,
            seed=12345678,
            chains=3,
            iter_warmup=1000,
            iter_sampling=1000,
        )
    
        ############################
        #  BUILD CLEAN OUTPUT DF   #
        ############################
        # Binary states (M × N)
        states = pd.DataFrame(
            L_clean.states,
            columns=L_clean.feature_names  # use original feature names
        )
    
        # Values (M × R)
        R = L_clean.values.shape[1]
        replicate_names = [f"rep_{i+1}" for i in range(R)]
        
        values = pd.DataFrame(
            L_clean.values,
            columns=replicate_names
        )
    
        # Combine states + values
        df_clean = pd.concat([states, values], axis=1)
    
        ############################
        #    SAVE CLEANED DATA     #
        ############################
        df_clean.to_csv(out_path, index=False)
    
        print(f"Saved batch-corrected dataset to: {out_path}")

############################
#       END OF SCRIPT      #
############################

print("\nAll datasets processed.")
