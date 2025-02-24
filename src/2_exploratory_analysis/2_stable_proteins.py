# This script identifies proteins whose expression is stable in baseline (uninfected)
# samples based on CV < 25%. It then considers and outputs how many patients
# a given protein is stable in. 

### [1] SET UP WORK SPACE ###
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set project root directory and paths
project_root = Path(__file__).parent.parent.parent
input_dir = project_root / "data" / "processed"
figures_dir = project_root / "data" / "figures" / "exploratory"
results_dir = project_root / "data" / "processed" / "results"

# Read data and identify baseline samples
data = pd.read_excel(input_dir /"4_pg_chivas_filtered.xlsx")
baseline_timepoints = ['screen', 'prior', '1mth', '3mth']
baseline_columns = [col for col in data.columns if any(tp in col for tp in baseline_timepoints)]
baseline_data = data[['Protein.Names'] + baseline_columns].copy()

### [2] DEFINE 'CALCULATE CV' FUNCTION ###
# Get patient IDs
patient_ids = sorted(list(set(int(col.split('_')[0]) for col in baseline_columns)))

def calculate_patient_cv(data, patient_id, baseline_timepoints):
    patient_cols = [col for col in data.columns 
                   if col.startswith(f"{patient_id}_") and 
                   any(tp in col for tp in baseline_timepoints)]
      
    patient_data = data[patient_cols].astype(float)
    
    cv_df = pd.DataFrame({
        'Protein.Names': data['Protein.Names'],
        'Mean': patient_data.mean(axis=1),
        'Std': patient_data.std(axis=1),
        'N_measurements': patient_data.notna().sum(axis=1)
    })
    
    cv_df['CV'] = (cv_df['Std'] / cv_df['Mean']) * 100
    return cv_df

### [3] APPLY FUNCTION AND PLOT ###
# Calculate CVs for each patient
patient_cvs = {}
for patient_id in patient_ids:
    cv_df = calculate_patient_cv(baseline_data, patient_id, baseline_timepoints)
    patient_cvs[patient_id] = cv_df

# Plot CV distributions
n_patients = len(patient_cvs)
n_cols = 4
n_rows = (n_patients + n_cols - 1) // n_cols

plt.figure(figsize=(20, 4*n_rows))

for i, (patient_id, cv_data) in enumerate(patient_cvs.items()):
    plt.subplot(n_rows, n_cols, i+1)
    sns.histplot(data=cv_data, x='CV', bins=50)
    
    for percentile, color in zip([25, 50, 75], ['r', 'g', 'b']):
        cutoff = cv_data['CV'].quantile(percentile/100)
        plt.axvline(x=cutoff, color=color, linestyle='--', 
                   label=f'{percentile}th percentile: {cutoff:.1f}')
    
    plt.title(f'Patient {patient_id}\nN={cv_data["N_measurements"].iloc[0]} measurements')
    plt.xlabel('Coefficient of Variation (%)')
    plt.ylabel('Count')
    plt.legend()
    plt.xlim(0, cv_data['CV'].quantile(0.95))

plt.tight_layout()
plt.savefig(figures_dir / 'participant_cv_distributions.svg', format='svg', dpi=300, bbox_inches='tight')
plt.show()

# Analyze protein stability
cv_threshold = 25
stable_proteins = {}
for patient_id, cv_data in patient_cvs.items():
    stable_in_patient = cv_data[cv_data['CV'] <= cv_threshold]['Protein.Names'].tolist()
    stable_proteins[patient_id] = set(stable_in_patient)

# Count stable proteins across patients
all_proteins = set(patient_cvs[patient_ids[0]]['Protein.Names'])
protein_stability_counts = {}
for protein in all_proteins:
    count = sum(1 for patient_proteins in stable_proteins.values() 
               if protein in patient_proteins)
    protein_stability_counts[protein] = count

# Plot stability distribution
cohort_sizes = range(1, len(patient_ids) + 1)
proteins_per_cohort = [sum(1 for count in protein_stability_counts.values() 
                         if count >= cohort_size) 
                      for cohort_size in cohort_sizes]

plt.figure(figsize=(12, 6))
plt.bar(cohort_sizes, proteins_per_cohort)
plt.xlabel('Number of Patients')
plt.ylabel('Number of Stable Proteins')
plt.title(f'Proteins with CV ≤ {cv_threshold}% Across Patients')
plt.grid(True, alpha=0.3)
plt.ylim(0, 1400)
plt.xticks(range(1, len(patient_ids) + 1))

for i, v in enumerate(proteins_per_cohort):
    plt.text(i + 1, v, str(v), ha='center', va='bottom')

plt.savefig(figures_dir / 'stability_distribution.svg', format='svg', dpi=300, bbox_inches='tight')
plt.show()

# Create and save stability summary
stability_table = pd.DataFrame({
    'Protein_Name': list(protein_stability_counts.keys()),
    'Number_of_Patients_Stable': list(protein_stability_counts.values())
})

stability_table = stability_table.sort_values(
    ['Number_of_Patients_Stable', 'Protein_Name'], 
    ascending=[False, True]
).reset_index(drop=True)

print("\nFirst 10 rows of protein stability table:")
print(stability_table.head(10))

stability_table.to_excel(results_dir / 'protein_stability_summary.xlsx', index=False)
