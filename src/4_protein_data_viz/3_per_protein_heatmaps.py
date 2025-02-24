# This file generates heatmaps per protein of interest where
# each column is a time point and 
# each row is a patient (hierarchically clustered)

### [1] SET UP WORKSPACE ###
# Import modules
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from sklearn.impute import KNNImputer
from matplotlib.patches import Rectangle

# Set up paths
project_root = Path(__file__).parent.parent.parent
input_dir = project_root / "data" / "processed"
results_dir = project_root / "data" / "processed" / "results"
figures_dir = project_root / "data" / "figures" / "heatmaps"  

# Read input files
sig_proteins = pd.read_excel(results_dir / "significant_proteins_99varcutoff.xlsx")
expression_data = pd.read_excel(input_dir / "4_pg_chivas_filtered.xlsx")

### [2] FILTER ABUNDANCE DATA BY SIGNIFICANT PROTEINS ###
# Merge two dataframes
filtered_df = pd.merge(sig_proteins, expression_data, on="Protein.Names", how="inner")

# Keep only necessary columns (Protein info + sample columns with underscores)
cols_to_keep = ["Protein.Names", "Genes", "Timepoint"] + [col for col in filtered_df.columns if '_' in col]
filtered_df = filtered_df[cols_to_keep]

# Get sample columns (those containing underscores)
sample_cols = [col for col in filtered_df.columns if '_' in col]

# Group sample columns by patient ID
patient_groups = {}
for col in sample_cols:
    patient_id = col.split('_')[0]
    if patient_id not in patient_groups:
        patient_groups[patient_id] = []
    patient_groups[patient_id].append(col)

### [3] CALCULATE AND SAVE Z-SCORES ###
# Calculate z-scores for each patient separately
for patient_id, columns in patient_groups.items():
    # Extract the data for this patient's samples
    patient_data = filtered_df[columns].values.astype(float)
    # Calculate z-scores for each protein across this patient's timepoints
    z_scores = stats.zscore(patient_data, axis=1, nan_policy='omit')
    # Replace the original values with z-scores
    filtered_df[columns] = z_scores

# Save results
output_path = results_dir / "99_poi_zscores.xlsx"
filtered_df.to_excel(output_path, index=False)

### [4] SET UP HEAT MAP ELEMENTS ###
# Define constants
TIMEPOINT_ORDER = ['screen', 'prior', 'chal', '24h', '48h', 'disch', '1week', '1mth', '3mth']
TIMEPOINT_LABELS = [
    'Screening', 'Evening\nbefore\nchallenge', '8h post-\nchallenge',
    '24h post-\nchallenge', '48h post-\nchallenge', 'Day of\ndischarge',
    '1 week\npost-\ndischarge', '1 month\npost-\ndischarge', '3 months\npost-\ndischarge'
]

PHARYNGITIS_DATA = {
    '75': 'Intense', '9': 'Intense', '20': 'Intense', '25': 'Medium', '43': 'Medium',
    '64': 'Medium', '55': 'Medium', '57': 'Medium', '61': 'Mild', '32': 'Mild',
    '17': 'Mild', '26': 'Mild', '71': 'Mild', '68': 'Mild', '33': 'Mild',
    '59': 'Mild', '2': 'Mild', '10': 'Mild', '13': 'None', '30': 'None'
}

GRADE_COLORS = {
    'Intense': '#433684',
    'Medium': '#734AE5',
    'Mild': '#C7B7F5',
    'None': '#DDDDDD'
}

# Define plotting functions
def setup_figure():
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.5, 0.2, 0.02, 20])
    return fig, (fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), 
                fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[0, 3]))

def create_heatmap(protein_data, protein_name):
    # Reshape data for heatmap
    heatmap_data = protein_data.pivot(index='Patient_ID', columns='Timepoint', values='Z_score')
    
    # Only use timepoints that exist in the data
    available_timepoints = [t for t in TIMEPOINT_ORDER if t in heatmap_data.columns]
    heatmap_data = heatmap_data.reindex(columns=available_timepoints)
    
    # Cluster data
    cluster_data = pd.DataFrame(
        KNNImputer(n_neighbors=2).fit_transform(heatmap_data),
        index=heatmap_data.index,
        columns=heatmap_data.columns
    )
    linkage = hierarchy.ward(pdist(cluster_data))
    
    # Create figure and axes
    fig, (ax_dend, ax_grade, ax_spacer, ax_main) = setup_figure()
    
    # Plot dendrogram
    dn = hierarchy.dendrogram(linkage, ax=ax_dend, orientation='left', 
                             no_labels=True, color_threshold=None, 
                             above_threshold_color='black')
    ax_dend.axis('off')
    
    # Reorder data based on clustering
    order = dn['leaves'][::-1]
    heatmap_data = heatmap_data.iloc[order]
    heatmap_data = heatmap_data.apply(pd.to_numeric, errors='coerce')
    
    # Create heatmap
    mask = heatmap_data.isnull()
    sns.heatmap(heatmap_data, ax=ax_main, cmap='RdBu_r', center=0, 
                vmin=-2, vmax=2, mask=mask, linewidths=0, 
                linecolor='none', cbar=False)
    
    # Add pharyngitis grade colors
    for i, patient in enumerate(heatmap_data.index):
        grade = PHARYNGITIS_DATA.get(str(patient), 'Unknown')
        ax_grade.add_patch(Rectangle((0, i), 1, 1, 
                                   facecolor=GRADE_COLORS.get(grade, '#FFFFFF'), 
                                   edgecolor='none'))
    
    # Configure axes
    ax_grade.set(xlim=(0, 1), ylim=(0, len(heatmap_data.index)))
    ax_grade.set_xticks([])
    ax_grade.set_yticks([])
    ax_grade.invert_yaxis()
    ax_spacer.axis('off')
    
    # Add missing value annotations
    for i, patient in enumerate(heatmap_data.index):
        for j, timepoint in enumerate(heatmap_data.columns):
            if pd.isnull(heatmap_data.iloc[i, j]):
                is_na = f"{patient}_{timepoint}" not in protein_data['Patient_Timepoint'].values
                text = 'N/A' if is_na else 'N.D.'
                color = 'white' if is_na else 'black'
                weight = 'bold' if is_na else 'normal'
                
                ax_main.text(j + 0.5, i + 0.5, text, ha='center', va='center',
                            color=color, fontsize=12, fontweight=weight,
                            fontname='Arial', fontstyle='italic' if is_na else 'normal')
                ax_main.add_patch(Rectangle((j, i), 1, 1, fill=True,
                                          facecolor='#b3b3b3', edgecolor='none'))
    
    # Style figure
    fig.suptitle(protein_name, fontname='Arial', fontweight='bold', fontsize=21)
    ax_main.set_xlabel("Timepoint", fontname='Arial', fontweight='bold', fontsize=18)
    fig.text(0.01, 0.5, "Patient ID", va='center', rotation='vertical',
             fontname='Arial', fontweight='bold', fontsize=18)
    
    # Set tick labels
    # Only show labels for available timepoints
    timepoint_indices = [TIMEPOINT_ORDER.index(t) for t in available_timepoints]
    ax_main.set_xticks(np.arange(len(available_timepoints)) + 0.5)
    ax_main.set_xticklabels([TIMEPOINT_LABELS[i] for i in timepoint_indices],
                           rotation=0, ha='center', va='top',
                           fontname='Arial', fontsize=14)
    ax_main.set_yticks(np.arange(len(heatmap_data.index)) + 0.5)
    ax_main.set_yticklabels(heatmap_data.index, va='center',
                           fontname='Arial', rotation=0, fontsize=14)
    
    plt.rcParams['font.family'] = 'Arial'
    plt.tight_layout(rect=[0.09, 0, 1, 0.95])
    
    return fig

### [5] EXECUTE FUNCTION WITH ERROR HANDLING ###
df = pd.read_excel(results_dir / "99_poi_zscores.xlsx")

for i, protein_name in enumerate(df['Protein.Names'].unique(), 1):
    print(f"Processing protein {i} of {len(df['Protein.Names'].unique())}: {protein_name}")
    
    try:
        # Prepare data
        protein_data = df[df['Protein.Names'] == protein_name].melt(
            id_vars=['Protein.Names', 'Genes'],
            value_vars=[col for col in df.columns if '_' in col],
            var_name='Patient_Timepoint',
            value_name='Z_score'
        )
        protein_data[['Patient_ID', 'Timepoint']] = (
            protein_data['Patient_Timepoint'].str.split('_', n=1, expand=True)
        )
        
        # Create and save heatmap
        fig = create_heatmap(protein_data, protein_name)
        fig.savefig(figures_dir / f"{protein_name.replace('/', '_')}_heatmap.svg",
                    format='svg', dpi=1000, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error processing {protein_name}: {str(e)}")
        plt.close()
        continue

print(f"All heatmaps have been generated and saved in: {figures_dir}")