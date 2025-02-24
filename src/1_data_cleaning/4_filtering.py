# This file removes proteins whose detection does not meet defined criteria:
# 1. Must be detected in at least three patients, where
# 2. "Detection" is defined as quantitated in at least three time-points

### [1] SET-UP WORKSPACE ###
# Import modules
import pandas as pd
import numpy as np
import re
from pathlib import Path
import random

# Set project root directory
project_root = Path(__file__).parent.parent.parent

# Define input/output paths
processed_dir = project_root / "data" / "processed"

# Load in data
df = pd.read_excel(processed_dir / "3_pg_chivas_no_keratin.xlsx")

### [2] FILTER IN IF DETECTED REPEATEDLY ###
# Function to extract patient ID from column name
def get_patient_id(column_name):
    match = re.match(r'(\d+)_', column_name)
    return match.group(1) if match else None

# Function to check if a protein meets the detection criteria for a single patient
def meets_criteria(patient_data):
    return patient_data.notna().sum(axis=1) >= 3

# Group the dataframe by patient
patient_ids = [get_patient_id(col) for col in df.columns]
patient_groups = [group for _, group in df.groupby(patient_ids, axis=1)]

# Apply the criteria to each patient group
meets_criteria_by_patient = pd.DataFrame([meets_criteria(group) for group in patient_groups]).T

# Count in how many patients each protein meets the criteria
patients_meeting_criteria = meets_criteria_by_patient.sum(axis=1)

# Filter proteins that meet the criteria in at least 3 patients
filtered_df = df[patients_meeting_criteria >= 3]

### [3] INSPECT RESULTS ###
print(f"Original number of proteins: {len(df)}")
print(f"Number of proteins after filtering: {len(filtered_df)}")

# Function to get detection counts for a protein
def get_detection_counts(protein_data):
    patient_groups = df.groupby(patient_ids, axis=1)
    return {patient: sum(protein_data[cols].notna()) 
            for patient, cols in patient_groups.groups.items()}

# Function to print protein information
def print_protein_info(protein_index, protein_name, counts, status):
    print(f"\n{protein_index} - {protein_name} ({status}):")
    for patient, count in counts.items():
        print(f"  Patient {patient}: detected in {count} time points")

# Get proteins filtered out and kept
proteins_filtered_out = set(df.index) - set(filtered_df.index)
proteins_kept = set(filtered_df.index)

# Sample proteins to display
sample_size = min(5, len(proteins_filtered_out), len(proteins_kept))
sample_filtered_out = random.sample(list(proteins_filtered_out), sample_size)
sample_kept = random.sample(list(proteins_kept), sample_size)

print("\nSample of proteins filtered out:")
for protein in sample_filtered_out:
    protein_name = df.loc[protein, 'Protein.Names']
    counts = get_detection_counts(df.loc[protein])
    print_protein_info(protein, protein_name, counts, "Filtered Out")

print("\nSample of proteins kept:")
for protein in sample_kept:
    protein_name = df.loc[protein, 'Protein.Names']
    counts = get_detection_counts(df.loc[protein])
    print_protein_info(protein, protein_name, counts, "Kept")

# Save the filtered dataframe to a new Excel file
filtered_df.to_excel(processed_dir / '4_pg_chivas_filtered.xlsx', index=False)




