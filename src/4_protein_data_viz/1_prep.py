# This script prepares and reshapes the data for heat-mapping

### [1] SET UP WORKSPACE ###
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
import numpy as np
import matplotlib.patches as patches
from pathlib import Path

# Set project root directory and paths
project_root = Path(__file__).parent.parent.parent
input_dir = project_root / "data" / "processed"
figures_dir = project_root / "data" / "figures" /"heatmaps"
results_dir = project_root / "data" / "processed" / "results"

# Load in data
abundance_df = pd.read_excel(input_dir / "5B_log2fc_abundance.xlsx")
poi_24h = pd.read_excel(results_dir/"percept_poi_24h.xlsx")
poi_48h = pd.read_excel(results_dir/"percept_poi_48h.xlsx")

### [2] PREPARE DATA FOR HEAT-MAPPING ###
# Get list of proteins significant at each timepoint
proteins_24h = set(poi_24h['Protein.Names'])
proteins_48h = set(poi_48h['Protein.Names'])

# Find proteins significant at both timepoints
significant_proteins = proteins_24h.union(proteins_48h)

print(f"Proteins significant at 24h: {len(proteins_24h)}")
print(f"Proteins significant at 48h: {len(proteins_48h)}")
print(f"Proteins significant at both timepoints: {len(significant_proteins)}")

# Filter abundance data for these proteins
poi_abundance = abundance_df[abundance_df['Protein.Names'].isin(significant_proteins)]

# Filter in relevant columns ("Genes", and containing "_48h")
selected_columns = ['Genes'] + [col for col in poi_abundance.columns if '_24h' in col or '_48h' in col]
poi_filtered = poi_abundance[selected_columns]

# Rename columns to remove '_48h' suffix
poi_wide_df = poi_filtered.rename(columns=lambda x: x.replace('_base_ratio', ''))

# Transpose
# First melt the dataframe to get it into long format
df_long = pd.melt(
    poi_wide_df, 
    id_vars=['Genes'], 
    var_name='Patient_ID', 
    value_name='Value'
)

# Create new column identifying timepoint
df_long['Timepoint'] = df_long['Patient_ID'].map(lambda x: '24h' if '24h' in x else '48h')

# Clean Patient_ID to remove timepoint reference
df_long['Patient_ID'] = df_long['Patient_ID'].str.replace('_24h', '').str.replace('_48h', '')

# Now pivot to get desired format
df_reshaped = df_long.pivot(
    index='Patient_ID',
    columns=['Genes', 'Timepoint'],
    values='Value'
)

# Flatten column names
df_reshaped.columns = [f'{timepoint}_{gene}' for gene, timepoint in df_reshaped.columns]

# Reset index to make Patient_ID a regular column
df_reshaped = df_reshaped.reset_index()

df_reshaped.to_excel(input_dir/"POI_log2fc_long.xlsx",index=False)

