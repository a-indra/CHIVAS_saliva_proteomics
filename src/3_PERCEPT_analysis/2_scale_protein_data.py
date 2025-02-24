# This script generates means and p-values with a 1-sample t-test
# in preparation for PERCEPT-scaling 
### [1] SET UP WORKSPACE ###
import pandas as pd
import numpy as np
from scipy.stats import ttest_1samp
from pathlib import Path

# Read the data
project_root = Path(__file__).parent.parent.parent
input_dir = project_root / "data" / "processed"
log2_fc_df = pd.read_excel(input_dir / "5B_log2fc_abundance.xlsx")

# Define timepoints of interest
experimental_timepoints = ['24h', '48h', 'disch']

# Calculate means for each timepoint
for timepoint in experimental_timepoints:
    # Adjust column selection to match your naming pattern
    timepoint_cols = [col for col in log2_fc_df.columns if f'_{timepoint}_base_ratio' in col]
    log2_fc_df[f'Mean_log2fc_{timepoint}'] = (
        log2_fc_df[timepoint_cols].mean(axis=1)
    )

# Function to calculate p-values (unchanged)
def calculate_p_value(row, cols):
    values = row[cols].dropna().astype(float)
    if len(values) < 2:
        return np.nan
    try:
        t_stat, p_val = ttest_1samp(values, 0)
        return p_val
    except:
        return np.nan

# Get all ratio columns
ratio_cols = [col for col in log2_fc_df.columns if 'base_ratio' in col]

# Ensure numeric data types
for col in ratio_cols:
    log2_fc_df[col] = pd.to_numeric(log2_fc_df[col], errors='coerce')

# Calculate p-values for each timepoint
for timepoint in experimental_timepoints:
    timepoint_cols = [col for col in ratio_cols if f'_{timepoint}_base_ratio' in col]
    log2_fc_df[f'P-value_{timepoint}'] = (
        log2_fc_df.apply(lambda row: calculate_p_value(row, timepoint_cols), axis=1)
    )

# Save results
log2_fc_df.to_excel(input_dir / "log2fc_protein.xlsx", index=False)

