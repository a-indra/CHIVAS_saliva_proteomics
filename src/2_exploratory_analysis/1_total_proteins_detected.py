# This script calculates and visualises the total number of filtered
# and unfiltered proteins detected commonly across the cohort 

### [1] SET UP WORKSPACE ###

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set project root directory
project_root = Path(__file__).parent.parent.parent

# Define input/output paths
processed_dir = project_root / "data" / "processed"
figures_dir = project_root / "data" / "figures" / "exploratory"

# Load the data
unfiltered_df = pd.read_excel(processed_dir / "3_pg_chivas_no_keratin.xlsx")
filtered_df = pd.read_excel(processed_dir / "4_pg_chivas_filtered.xlsx")

### [2] WRITE FUNCTION TO CALCULATE CUMULATIVE/COMMON COUNTS ###
def visualize_protein_detection(df, title):
    # Identify the sample columns
    sample_columns = [col for col in df.columns if '_' in col]

    # Function to get patient ID from column name
    def get_patient_id(column_name):
        return int(column_name.split('_')[0])

    # Group columns by patient
    patient_groups = {}
    for col in sample_columns:
        patient_id = get_patient_id(col)
        if patient_id not in patient_groups:
            patient_groups[patient_id] = []
        patient_groups[patient_id].append(col)

    # Count in how many patients each protein is detected
    df['detection_count'] = df.apply(lambda row: sum(row[patient_groups[patient]].any() for patient in patient_groups), axis=1)

    # Calculate cumulative counts
    total_proteins = len(df)
    max_patients = len(patient_groups)
    cumulative_counts = [(i, (df['detection_count'] >= i).sum()) for i in range(1, max_patients + 1)]

    # Create a DataFrame for plotting
    plot_df = pd.DataFrame(cumulative_counts, columns=['Patients', 'Proteins'])

    # Create the column graph
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.barplot(x='Patients', y='Proteins', data=plot_df, ax=ax)

    ax.set_title(title)
    ax.set_xlabel('Number of Patients')
    ax.set_ylabel('Number of Proteins Detected')

    # Add value labels on top of each bar
    for i, row in plot_df.iterrows():
        ax.text(i, row['Proteins'], str(row['Proteins']), ha='center', va='bottom')

    plt.tight_layout()
    plt.show()

    # Print out the data for verification
    print(f"Data for {title}:")
    print(plot_df)
    
    return fig  

### [3] GENERATE AND SAVE FIGURES ###
# Apply the visualization to both dataframes
filtered_fig = visualize_protein_detection(filtered_df, "Protein Detection Across Patients (Filtered Data)")
unfiltered_fig = visualize_protein_detection(unfiltered_df, "Protein Detection Across Patients (Unfiltered Data)")

# Save figures as SVG with proper paths
filtered_fig.savefig(figures_dir / "protein_detection_filtered.svg", format='svg', dpi=300, bbox_inches='tight')
unfiltered_fig.savefig(figures_dir / "protein_detection_unfiltered.svg", format='svg', dpi=300, bbox_inches='tight')

