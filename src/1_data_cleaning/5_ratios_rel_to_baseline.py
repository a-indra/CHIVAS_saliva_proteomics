# This script transforms the abundance data such that it is 
# expressed as a ratio, relative to baseline samples.
# These ratios are then log2-transformed.

### [1] SET UP WORKSPACE ###
import pandas as pd
import numpy as np
from pathlib import Path

# Set project root directory and paths
project_root = Path(__file__).parent.parent.parent
input_dir = project_root / "data" / "processed"
figures_dir = project_root / "data" / "figures" / "percept"
results_dir = project_root / "data" / "processed" / "results"

# Load in data
abundance_df= pd.read_excel(input_dir/"4_pg_chivas_filtered.xlsx")

### [2] CONVERT ABUNDANCE DATA TO RATIOS ###
# Define baseline suffixes
baseline_suffixes = ['_screen', '_prior', '_1mth', '_3mth']

# Identify unique patient identifiers
patient_ids = set(col.split('_')[0] for col in abundance_df.columns if col.split('_')[0].isdigit())

# Initialize a dictionary to store new baseline columns
new_baseline_cols = {}

# Calculate the mean for each patient and each baseline suffix
for patient_id in patient_ids:
    baseline_cols = [f"{patient_id}{suffix}" for suffix in baseline_suffixes if f"{patient_id}{suffix}" in abundance_df.columns]
    if baseline_cols:
        new_baseline_cols[f'{patient_id}_base'] = abundance_df[baseline_cols].mean(axis=1, skipna=True)

# Create a DataFrame from the new baseline columns
baseline_means_df = pd.DataFrame(new_baseline_cols)

# Merge the new baseline means DataFrame with the original DataFrame
abundance_with_baseline_df = pd.concat([abundance_df, baseline_means_df], axis=1)

# Compute the ratios between 24h, 48h, discharge and base
ratios = {}
for patient_id in patient_ids:
    base_col = f'{patient_id}_base'
    if base_col in abundance_with_baseline_df.columns:
        if f'{patient_id}_24h' in abundance_with_baseline_df.columns:
            ratios[f'{patient_id}_24h_base_ratio'] = abundance_with_baseline_df[f'{patient_id}_24h'] / abundance_with_baseline_df[base_col]
        if f'{patient_id}_48h' in abundance_with_baseline_df.columns:
            ratios[f'{patient_id}_48h_base_ratio'] = abundance_with_baseline_df[f'{patient_id}_48h'] / abundance_with_baseline_df[base_col]
        if f'{patient_id}_disch' in abundance_with_baseline_df.columns:
            ratios[f'{patient_id}_disch_base_ratio'] = abundance_with_baseline_df[f'{patient_id}_disch'] / abundance_with_baseline_df[base_col]

# Create a DataFrame for the ratios
ratios_df = pd.DataFrame(ratios)

# Combine the ratios DataFrame with the protein name and ID columns
final_ratios_df = pd.concat([abundance_with_baseline_df[['Protein.Names', 'Protein.Ids','Genes']], ratios_df], axis=1)
final_ratios_df.to_excel(input_dir/"5A_protein_abundance_ratios.xlsx", index=False)

### [3] LOG-TRANSFORM VALUES ###
# Convert ratios to log2 fold changes
log2_fc_df = final_ratios_df.copy()
ratio_cols = [col for col in log2_fc_df.columns if col.endswith('_base_ratio')]
log2_fc_df[ratio_cols] = np.log2(log2_fc_df[ratio_cols])

log2_fc_df.to_excel(input_dir/"5B_log2fc_abundance.xlsx", index=False)