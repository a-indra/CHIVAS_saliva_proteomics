# This file generates a 24- and 48-hour clustered heatmap of proteins of interest
# where each column is a patient 
# and each row is a protein

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
figures_dir = project_root / "data" / "figures" /"heatmaps"
input_dir = project_root / "data" / "processed"
results_dir = project_root / "data" / "processed" / "results"

# Read the data
df = pd.read_excel(input_dir/"POI_log2fc_long.xlsx")
patient_info = pd.read_excel(project_root/"data"/"patient_pharyngitis_grade.xlsx")
patient_info = patient_info.set_index('Patient_ID')

### [2] DEFINE PREPATORY FUNCTIONS ###
def prepare_dataframe(df, time_point):
    cols = ['Patient_ID'] + [col for col in df.columns if col.startswith(f'{time_point}_')]
    new_df = df[cols].copy()
    new_df.columns = ['Patient_ID'] + [col.replace(f'{time_point}_', '') for col in new_df.columns if col != 'Patient_ID']
    return new_df

def process_dataframe(df):
    df = df.set_index('Patient_ID')
    df = df.dropna(how='all')
    df_t = df.T
    nan_threshold = 0.50  # Kept as in your preferred version
    df_t = df_t.dropna(thresh=int(df_t.shape[1] * (1 - nan_threshold)))
    df_t = df_t.fillna(0)
    return df_t

### [3] DEFINE FUNCTION FOR CLUSTERED HEATMAP ###
def create_clustered_heatmap(data, title, patient_info, filename):
    # Perform clustering using Ward linkage and Euclidean distance
    row_linkage = hierarchy.linkage(data, method='ward', metric='euclidean')
    col_linkage = hierarchy.linkage(data.T, method='ward', metric='euclidean')
    
    # Create the heatmap
    g = sns.clustermap(data, 
                       cmap='RdBu_r', 
                       vmin=-3, 
                       vmax=3, 
                       center=0,
                       row_linkage=row_linkage,
                       col_linkage=col_linkage,
                       xticklabels=True,
                       yticklabels=True,
                       figsize=(20, 15))
    
    # Customize the plot
    g.fig.suptitle(title, fontsize=16, y=1.02)
    
    # Rotate x-axis labels
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    
    # Add pharyngitis grade coloring
    patient_order = data.columns[g.dendrogram_col.reordered_ind]
    grades = patient_info.loc[patient_order, 'P_grade_num']
    
    # Define colors for each grade
    grade_colors = {1: '#fee5d9', 2: '#fcae91', 3: '#fb6a4a', 4: '#cb181d'}
    
    # Get the positions of the heatmap columns
    xpos = g.ax_heatmap.get_xticks()
    
    # Create a new axes for the grade colors
    ax_grade = g.fig.add_axes([g.ax_heatmap.get_position().x0, 
                               g.ax_heatmap.get_position().y1,
                               g.ax_heatmap.get_position().width, 0.02])
    
    # Plot colored rectangles for each patient
    for x, grade in zip(xpos, grades):
        rect = patches.Rectangle((x, 0), 1, 1, facecolor=grade_colors.get(grade, 'white'))
        ax_grade.add_patch(rect)
    
    ax_grade.set_xlim(xpos[0], xpos[-1] + 1)
    ax_grade.set_ylim(0, 1)
    ax_grade.axis('off')
    
    # Add a legend
    legend_elements = [patches.Patch(facecolor=color, edgecolor='black',
                                     label=f'Grade {grade}')
                       for grade, color in grade_colors.items()]
    g.fig.legend(handles=legend_elements, title='Pharyngitis Grade',
                 loc='upper right', bbox_to_anchor=(0.98, 0.98))
    
    # Adjust layout to prevent cutting off labels
    g.fig.tight_layout()
    
    # Save the plot as SVG
    g.savefig(filename, format='svg', dpi=300, bbox_inches='tight')
    print(f"Heatmap saved as {filename}")
    
    # Close the figure to free up memory
    plt.close()

### [4] EXECUTE MAIN FUNCTION ###
if __name__ == "__main__":
    # Prepare separate dataframes for 24h and 48h
    df_24h = prepare_dataframe(df, '24h')
    df_48h = prepare_dataframe(df, '48h')

    # Save the new dataframes to Excel files
    df_24h.to_excel(results_dir/'patient_POI_ratio_log2_24h.xlsx', index=False)
    df_48h.to_excel(results_dir/'patient_POI_ratio_log2_48h.xlsx', index=False)

    # Process dataframes
    df_24h_processed = process_dataframe(df_24h)
    df_48h_processed = process_dataframe(df_48h)

    # Print statistics
    for name, df in [('24h', df_24h_processed), ('48h', df_48h_processed)]:
        print(f"\n{name} data:")
        print(f"Number of proteins after filtering: {df.shape[0]}")
        print(f"Number of patients: {df.shape[1]}")

    # Create and save heatmaps
    create_clustered_heatmap(df_24h_processed, 'Clustered Protein Expression at 24 Hours\nWard Linkage, Euclidean Distance', patient_info, figures_dir/'heatmap_24h.svg')
    create_clustered_heatmap(df_48h_processed, 'Clustered Protein Expression at 48 Hours\nWard Linkage, Euclidean Distance', patient_info, figures_dir/'heatmap_48h.svg')