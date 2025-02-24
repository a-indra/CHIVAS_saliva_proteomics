# This script applies pre-calculated PERCEPT thresholds to identify proteins of interest (POI)
# Note: PERCEPT threshold values are hardcoded for reproducibility but can be modified if needed

### [1] SET UP WORKSPACE ###
import pandas as pd
import numpy as np
from scipy.stats import ttest_1samp
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import seaborn as sns
import kaleido
from pathlib import Path

# Set project root directory and paths
project_root = Path(__file__).parent.parent.parent
input_dir = project_root / "data" / "processed"
figures_dir = project_root / "data" / "figures" / "percept"
results_dir = project_root / "data" / "processed" / "results"

# Load in data
log2_fc_df = pd.read_excel(input_dir / "log2fc_protein.xlsx")
ratio_cols = [col for col in log2_fc_df.columns if col.endswith('_base_ratio')]

### [2] DEFINE FUNCTIONS ###
# PERCEPT scaling function
def percept(m0, m1, F, p):
   return m0 + ((m0 - m1) / -(F**p))

# Function to apply PERCEPT scaling
def apply_percept(data, hypothetical_mean, penalty):
   tval, pval = ttest_1samp(data.dropna(), popmean=hypothetical_mean)
   sample_mean = np.mean(data.dropna())
   scaled_value = percept(m0=hypothetical_mean, m1=sample_mean, F=penalty, p=pval)
   return scaled_value, pval

# Function to create scatter plot
def create_percept_scatter(df, x_col, y_col, p_col, title, percept_lower, percept_upper):
    # Color based on PERCEPT thresholds
    df['color'] = 'gray'
    df.loc[df[y_col] > percept_upper, 'color'] = 'blue'   # "Upregulated"
    df.loc[df[y_col] < percept_lower, 'color'] = 'red'    # "Downregulated"

    fig = go.Figure()

    for color in ['gray', 'blue', 'red']:
       mask = df['color'] == color
       fig.add_trace(go.Scatter(
           x=df.loc[mask, x_col],
           y=df.loc[mask, y_col],
           mode='markers',
           marker=dict(
               color=color,
               size=5,
               opacity=0.7
           ),
           name=color,
           text=df.loc[mask, 'Protein.Names'],
           hovertemplate='<b>%{text}</b><br>' +
                         f'{x_col}: ' + '%{x:.3f}<br>' +
                         f'{y_col}: ' + '%{y:.3f}<br>' +
                         f'{p_col}: ' + '%{customdata:.3e}',
           customdata=df.loc[mask, p_col]
       ))

    fig.update_layout(
       title=title,
       xaxis_title=x_col,
       yaxis_title=y_col,
       hovermode="closest"
   )
     
    return fig

# Function to extract significant proteins
def extract_significant_proteins(df, percept_col, log2fc_col, pvalue_col, lower_threshold, upper_threshold):
   return df[(df[percept_col] < lower_threshold) | (df[percept_col] > upper_threshold)][[
       'Protein.Names', 'Protein.Ids', percept_col, log2fc_col, pvalue_col
   ]]

### [3] APPLY PERCEPT SCALING ###
# Set parameters
m0 = 0
n = 20
F = 10 * n

# Define PERCEPT thresholds
percept_thresholds = {
    '95': (-0.384025325, 0.356768907),
    '99': (-0.73164318, 0.682285458)
}

# Apply scaling for each timepoint
for timepoint in ['24h', '48h', 'disch']:
   cols = [col for col in ratio_cols if timepoint in col]
   results = log2_fc_df[cols].apply(lambda row: apply_percept(row, m0, F), axis=1)
   log2_fc_df[f'PERCEPT_scaled_{timepoint}'] = results.apply(lambda x: x[0])
   log2_fc_df[f'P-value_{timepoint}'] = results.apply(lambda x: x[1])

PERCEPT_df = log2_fc_df
PERCEPT_df.to_excel(results_dir / "PERCEPT_scaled_all_proteins.xlsx", index=False)

### [4] CREATE PLOTS AND EXTRACT PROTEINS OF INTEREST ###
for timepoint in ['24h', '48h', 'disch']:
   output_file = results_dir / f"percept_poi_{timepoint}.xlsx"
   
   with pd.ExcelWriter(output_file) as writer:
       for variance, (lower, upper) in percept_thresholds.items():
           # Create scatter plot
           scatter = create_percept_scatter(PERCEPT_df, 
                                         f'Mean_log2fc_{timepoint}', 
                                         f'PERCEPT_scaled_{timepoint}', 
                                         f'P-value_{timepoint}', 
                                         f'PERCEPT-scaled vs Mean Log2FC ({timepoint}) - {variance}% variance',
                                         lower,
                                         upper)
           
           # Save scatter plot
           scatter.write_html(figures_dir / f"percept_scatter_{timepoint}_{variance}percent.html")
           
           # Display the plot for 95% threshold
           if variance == '95':
              scatter.show()

           # Extract significant proteins
           significant_proteins = extract_significant_proteins(PERCEPT_df, 
                                                            f'PERCEPT_scaled_{timepoint}', 
                                                            f'Mean_log2fc_{timepoint}', 
                                                            f'P-value_{timepoint}',
                                                            lower,
                                                            upper)
           
           # Save to Excel sheet
           significant_proteins.to_excel(writer, sheet_name=f'{variance}varcutoff', index=False)

### [5] EXTRACT AND COMBINE SIGNIFICANT PROTEINS ACROSS TIME-POINTS, PER THRESHOLD ###
# Define input files and their timepoints
files = {
    "percept_poi_24h.xlsx": "24h",
    "percept_poi_48h.xlsx": "48h",
    "percept_poi_disch.xlsx": "Disch"
}

# Process both thresholds
for threshold in ["95varcutoff", "99varcutoff"]:
    # List to store DataFrames for each timepoint
    dfs = []
    
    # Process each file
    for filename, timepoint in files.items():
        filepath = results_dir / filename
        if filepath.exists():
            df = pd.read_excel(filepath, sheet_name=threshold)
            df_timepoint = pd.DataFrame({
                "Protein.Names": df["Protein.Names"],
                "Timepoint": timepoint
            })
            dfs.append(df_timepoint)
    
    # Combine all DataFrames
    if dfs:
        combined_df = pd.concat(dfs, ignore_index=True)
    else:
        combined_df = pd.DataFrame(columns=["Protein.Names", "Timepoint"])
    
    # Group by protein and aggregate timepoints
    grouped = combined_df.groupby("Protein.Names")["Timepoint"].agg(lambda x: ','.join(sorted(x)))
    result_df = grouped.reset_index()
    
    # Save results
    output_path = results_dir / f"significant_proteins_{threshold}.xlsx"
    result_df.to_excel(output_path, index=False)
    print(f"Results saved to: {output_path}")