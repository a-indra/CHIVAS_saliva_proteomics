# This file generates individual lineplots for each protein of interest 

### [1] SET UP WORKSPACE ###
# Import modules
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# Set up paths
project_root = Path(__file__).parent.parent.parent
results_dir = project_root / "data" / "processed" / "results"
figures_dir = project_root / "data" / "figures" / "lineplots"
figures_dir.mkdir(parents=True, exist_ok=True)

### [2] DEFINE FUNCTION ###
# Define timepoint parameters
TIMEPOINT_ORDER = ['screen', 'prior', 'chal', '24h', '48h', 'disch', '1week', '1mth', '3mth']
TIMEPOINT_LABELS = [
    'Screening',
    'Evening\nbefore\nchallenge',
    '8h post-\nchallenge',
    '24h post-\nchallenge',
    '48h post-\nchallenge',
    'Day of\ndischarge',
    '1 week\npost-\ndischarge',
    '1 month\npost-\ndischarge',
    '3 months\npost-\ndischarge'
]
TIMEPOINT_MAP = dict(zip(TIMEPOINT_ORDER, TIMEPOINT_LABELS))

def plot_patient_data(ax, data, timepoints, label, marker, color):
    x_vals = []
    y_vals = []
    first_segment = True  # Track if this is the first segment
    
    for tp in TIMEPOINT_ORDER:
        if tp in timepoints:
            val = data.get(tp)
            if pd.notnull(val):
                x_vals.append(TIMEPOINT_MAP[tp])
                y_vals.append(val)
            else:
                if x_vals:  # Plot accumulated points
                    ax.plot(x_vals, y_vals, 
                           label=label if first_segment else "_nolegend_",
                           marker=marker, color=color, linestyle='-')
                    first_segment = False
                x_vals = []
                y_vals = []
                
    if x_vals:  # Plot any remaining points
        ax.plot(x_vals, y_vals, 
               label=label if first_segment else "_nolegend_",
               marker=marker, color=color, linestyle='-')

def create_lineplot(protein_data, protein_name):
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Process cohort data (excluding patients 13 and 30)
    cohort_data = {}
    cohort_std = {}
    for tp in TIMEPOINT_ORDER:
        cols = [col for col in protein_data.columns if '_' + tp in col 
                and not col.startswith('13_') and not col.startswith('30_')]
        if cols:
            values = protein_data[cols].values[0]
            values = values[~np.isnan(values)]
            if len(values) > 0:
                cohort_data[tp] = np.mean(values)
                cohort_std[tp] = np.std(values)
    
    # Plot cohort average
    if cohort_data:
        timepoints = list(cohort_data.keys())
        values = [cohort_data[tp] for tp in timepoints]
        std_values = [cohort_std[tp] for tp in timepoints]
        timepoint_labels = [TIMEPOINT_MAP[tp] for tp in timepoints]
        
        # Plot mean line
        ax.plot(timepoint_labels, values, color='blue', marker='o',
                label='Cohort average', zorder=2)
        
        # Add confidence interval
        ax.fill_between(timepoint_labels, 
                       np.array(values) - np.array(std_values),
                       np.array(values) + np.array(std_values),
                       color='blue', alpha=0.2, zorder=1)
    
    # Process and plot patient data
    for patient_id, color, marker in [('13', 'orange', 's'), ('30', 'green', '^')]:
        patient_data = {}
        for tp in TIMEPOINT_ORDER:
            cols = [col for col in protein_data.columns if col.startswith(f'{patient_id}_{tp}')]
            if cols:
                val = protein_data[cols].values[0][0]
                patient_data[tp] = val
        
        plot_patient_data(ax, patient_data, patient_data.keys(), 
                         f'Patient {patient_id}', marker, color)
    
    # Style the plot
    ax.set_title(f"{protein_name} Expression", fontsize=14, pad=20)
    ax.set_xlabel("Time-point", fontsize=12, labelpad=10)
    ax.set_ylabel("Protein Expression (Z-score)", fontsize=12, labelpad=10)
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Set x-axis
    plt.xticks(range(len(TIMEPOINT_LABELS)), TIMEPOINT_LABELS, rotation=45, ha='right')
    plt.tight_layout()
    
    return fig

### [3] EXECUTE LINEPLOT FUNCTION 
# Read in data
df = pd.read_excel(results_dir / "99_poi_zscores.xlsx")

# Execute lineplots
for i, protein_name in enumerate(df['Protein.Names'].unique(), 1):
    print(f"Processing protein {i} of {len(df['Protein.Names'].unique())}: {protein_name}")
    
    try:
        protein_data = df[df['Protein.Names'] == protein_name]
        fig = create_lineplot(protein_data, protein_name)
        
        output_path = figures_dir / f"{protein_name.replace('/', '_')}_lineplot.svg"
        fig.savefig(output_path, format='svg', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Successfully saved: {protein_name}")
        
    except Exception as e:
        print(f"Error processing {protein_name}: {str(e)}")
        plt.close()
        continue

print(f"All lineplots have been generated and saved in: {figures_dir}")