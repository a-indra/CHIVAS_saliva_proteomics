# This script pulls in the .tsv generated by Fragpipe and 
# cleans up the protein group file (pg_matrix) by renaming
# the columns and removing duplicates

### [1] SET-UP WORKSPACE ###
# Import modules
import pandas as pd
import numpy as np
import re
from pathlib import Path

# Set project root directory
project_root = Path(__file__).parent.parent.parent

# Define input/output paths
input_file = project_root / "data" / "report.pg_matrix.tsv"
processed_dir = project_root / "data" / "processed"

# Read the TSV file and save as Excel for inspection
pg_matrix = pd.read_csv(input_file, sep='\t')
pg_matrix.to_excel(processed_dir / "0_pg_matrix.xlsx", index=False)

# Read in Excel file
pg_matrix = pd.read_excel(processed_dir / "0_pg_matrix.xlsx")

### [2] REPLACE FILE NAMES ###
# Function to rename columns
def rename_column(col):
    # Split the string by '/' and take the last part
    new_name = col.split('/')[-1]
    # Remove '.mzML' from the end
    new_name = new_name.replace('.mzML', '')
    return new_name

# Get the list of columns that need to be renamed (excluding the first 5 metadata columns)
columns_to_rename = pg_matrix.columns[5:]

# Create a dictionary mapping old names to new names
rename_dict = {col: rename_column(col) for col in columns_to_rename}

# Rename the columns
pg_matrix.rename(columns=rename_dict, inplace=True)

# Print the new column names to verify
print(pg_matrix.columns)

# Change shortened file-name to sample ID
def rename_column(col):
    parts = col.split('_')
    
    # Category 1: "19AUG24_AI_CHIVAS_68-4"
    if col.startswith("19AUG24_AI_CHIVAS_"):
        return parts[-1]
    
    # For other categories
    if len(parts) > 6:  # This is a rerun
        sample_num = parts[2].replace('Sample', '')
        run_num = parts[5]
        return f"{sample_num}-{run_num}-rerun"
    elif len(parts) == 6:  # Regular run
        sample_num = parts[2].replace('Sample', '')
        run_num = parts[5]
        return f"{sample_num}-{run_num}"
    
    # If it doesn't match any category, return the original name
    return col

# Get the list of columns that need to be renamed (excluding the first 5 metadata columns)
columns_to_rename = pg_matrix.columns[5:]

# Create a dictionary mapping old names to new names
rename_dict = {col: rename_column(col) for col in columns_to_rename}

# Rename the columns
pg_matrix.rename(columns=rename_dict, inplace=True)

# Function to count valid IDs in a column (non-empty/non-NaN values)
def count_valid_ids(df, column):
    return df[column].notna().sum()

# Find all columns with '-rerun' suffix and their original counterparts
rerun_cols = [col for col in pg_matrix.columns if col.endswith('-rerun')]
original_cols = [col.replace('-rerun', '') for col in rerun_cols]

print("\nFound these rerun pairs:")
for rerun, orig in zip(rerun_cols, original_cols):
    if orig in pg_matrix.columns:
        rerun_count = count_valid_ids(pg_matrix, rerun)
        orig_count = count_valid_ids(pg_matrix, orig)
        print(f"{rerun}: {rerun_count} IDs, {orig}: {orig_count} IDs")

# Columns to drop (original runs with fewer IDs)
cols_to_drop = []
for rerun, orig in zip(rerun_cols, original_cols):
    if orig in pg_matrix.columns:
        rerun_count = count_valid_ids(pg_matrix, rerun)
        orig_count = count_valid_ids(pg_matrix, orig)
        if rerun_count > orig_count:
            cols_to_drop.append(orig)
        else:
            cols_to_drop.append(rerun)

# Drop the columns with fewer IDs
print("\nDropping these columns:", cols_to_drop)
pg_matrix = pg_matrix.drop(columns=cols_to_drop)

# Handle duplicate column names
duplicate_cols = pg_matrix.columns[pg_matrix.columns.duplicated()].unique()
if len(duplicate_cols) > 0:
    print("\nFound duplicate columns:", duplicate_cols)
    for col in duplicate_cols:
        # Get all columns with this name
        dup_data = pg_matrix.loc[:, pg_matrix.columns == col]
        # Keep the column with the most valid IDs
        valid_counts = dup_data.notna().sum()
        keep_col = valid_counts.idxmax()
        drop_cols = [i for i in range(len(pg_matrix.columns)) 
                    if pg_matrix.columns[i] == col and i != keep_col]
        pg_matrix = pg_matrix.drop(pg_matrix.columns[drop_cols], axis=1)

# Check for any remaining rerun suffixes
remaining_reruns = [col for col in pg_matrix.columns if col.endswith('-rerun')]
if remaining_reruns:
    print("\nThese columns still have -rerun suffix:", remaining_reruns)
    
    # Remove the -rerun suffix
    new_names = {col: col.replace('-rerun', '') for col in remaining_reruns}
    pg_matrix = pg_matrix.rename(columns=new_names)
    print("Removed -rerun suffix from these columns")

# Sort columns
def custom_sort_key(col_name):
    if col_name in ['Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description']:
        return (0, col_name)  # Keep these columns first
    
    parts = col_name.split('-')
    if len(parts) >= 2:
        sample_num = int(parts[0])
        run_num = int(parts[1].split('-')[0])  # In case there's a "-rerun"
        return (1, sample_num, run_num)
    
    return (2, col_name)  # For any columns that don't match the expected format

sorted_columns = sorted(pg_matrix.columns, key=custom_sort_key)

# Reorder the DataFrame columns
pg_matrix_sorted = pg_matrix.reindex(columns=sorted_columns)
print(pg_matrix_sorted.columns)

### [3] RENAME TIME-POINTS WITH MORE DESCRIPTIVE TITLES ###
# Create the time_point_mapping dictionary
time_point_mapping = {
    "2-1": "2_screen", "2-2": "2_prior", "2-3": "2_chal", "2-4": "2_24h", "2-5": "2_48h", "2-6": "2_disch", "2-7": "2_1week", "2-8": "2_1mth", "2-9": "2_3mth",
    "9-1": "9_screen", "9-2": "9_prior", "9-3": "9_chal", "9-4": "9_24h", "9-5": "9_disch", "9-6": "9_1week", "9-7": "9_1mth", "9-8": "9_3mth",
    "10-1": "10_screen", "10-2": "10_prior", "10-3": "10_chal", "10-4": "10_24h", "10-5": "10_48h", "10-6": "10_72h", "10-7": "10_disch", "10-8": "10_1week", "10-9": "10_1mth", "10-10": "10_3mth",
    "13-1": "13_screen", "13-2": "13_prior", "13-3": "13_chal", "13-4": "13_24h", "13-5": "13_48h", "13-6": "13_72h", "13-7": "13_96h", "13-8": "13_disch", "13-9": "13_1week", "13-10": "13_1mth",
    "17-1": "17_screen", "17-2": "17_prior", "17-3": "17_chal", "17-4": "17_24h", "17-5": "17_48h", "17-6": "17_disch", "17-7": "17_1week",
    "20-1": "20_screen", "20-2": "20_prior", "20-3": "20_chal", "20-4": "20_24h", "20-5": "20_48h", "20-6": "20_disch", "20-7": "20_1week", "20-8": "20_1mth", "20-9": "20_3mth",
    "25-1": "25_screen", "25-2": "25_prior", "25-3": "25_chal", "25-4": "25_24h", "25-5": "25_48h", "25-6": "25_disch", "25-7": "25_1week", "25-8": "25_1mth", "25-9": "25_3mth",
    "26-1": "26_screen", "26-2": "26_prior", "26-3": "26_chal", "26-4": "26_24h", "26-5": "26_48h", "26-6": "26_72h", "26-7": "26_96h", "26-8": "26_disch", "26-9": "26_1week", "26-10": "26_1mth", "26-11": "26_3mth",
    "30-1": "30_screen", "30-2": "30_prior", "30-3": "30_chal", "30-4": "30_24h", "30-5": "30_48h", "30-6": "30_72h", "30-7": "30_96h", "30-8": "30_disch", "30-9": "30_1week", "30-10": "30_1mth", "30-11": "30_3mth",
    "32-1": "32_screen", "32-2": "32_prior", "32-3": "32_chal", "32-4": "32_24h", "32-5": "32_48h", "32-6": "32_disch", "32-7": "32_1week", "32-8": "32_1mth", "32-9": "32_3mth",
    "33-1": "33_screen", "33-2": "33_prior", "33-3": "33_chal", "33-4": "33_24h", "33-5": "33_48h", "33-6": "33_disch", "33-7": "33_1week", "33-8": "33_1mth", "33-9": "33_3mth",
    "43-1": "43_screen", "43-2": "43_prior", "43-3": "43_chal", "43-4": "43_24h", "43-5": "43_disch", "43-6": "43_1week", "43-7": "43_1mth", "43-8": "43_3mth",
    "55-1": "55_screen", "55-2": "55_prior", "55-3": "55_chal", "55-4": "55_24h", "55-5": "55_48h", "55-6": "55_72h", "55-7": "55_disch", "55-8": "55_1week", "55-9": "55_1mth", "55-10": "55_3mth",
    "57-1": "57_screen", "57-2": "57_prior", "57-3": "57_chal", "57-4": "57_24h", "57-5": "57_48h", "57-6": "57_72h", "57-7": "57_96h", "57-8": "57_disch", "57-9": "57_1week", "57-10": "57_1mth", "57-11": "57_3mth",
    "59-1": "59_screen", "59-2": "59_prior", "59-3": "59_chal", "59-4": "59_24h", "59-5": "59_48h", "59-6": "59_disch", "59-7": "59_1week", "59-8": "59_1mth", "59-9": "59_3mth",
    "61-1": "61_screen", "61-2": "61_prior", "61-3": "61_chal", "61-4": "61_24h", "61-5": "61_48h", "61-6": "61_disch", "61-7": "61_1week", "61-8": "61_1mth", "61-9": "61_3mth",
    "64-1": "64_screen", "64-2": "64_prior", "64-3": "64_chal", "64-4": "64_24h", "64-5": "64_48h", "64-6": "64_disch", "64-7": "64_1week", "64-8": "64_1mth", "64-9": "64_3mth",
    "68-1": "68_screen", "68-2": "68_prior", "68-3": "68_chal", "68-4": "68_24h", "68-5": "68_disch", "68-6": "68_1week", "68-7": "68_1mth",
    "71-1": "71_screen", "71-2": "71_prior", "71-3": "71_chal", "71-4": "71_24h", "71-5": "71_48h", "71-6": "71_disch", "71-7": "71_1week", "71-8": "71_1mth",
    "75-1": "75_screen", "75-2": "75_prior", "75-3": "75_chal", "75-4": "75_24h", "75-5": "75_48h", "75-6": "75_disch", "75-7": "75_1week"
}

def rename_column(col_name):
    if col_name in ['Protein.Group', 'Protein.Ids', 'Protein.Names', 'Genes', 'First.Protein.Description']:
        return col_name  # Keep these columns as they are
    if col_name in time_point_mapping:
        return time_point_mapping[col_name]
    return col_name  # If no mapping found, keep the original name

# Create a list of tuples containing old and new column names
old_and_new_columns = [(col, rename_column(col)) for col in pg_matrix_sorted.columns]

# Rename the columns
pg_matrix_sorted.columns = [new_col for _, new_col in old_and_new_columns]

# Print the old and new column names side by side
print("Old Name".ljust(30) + "New Name")
print("-" * 60)
for old_col, new_col in old_and_new_columns:
    print(f"{old_col.ljust(30)}{new_col}")

# Print summary of changes
unchanged = sum(1 for old, new in old_and_new_columns if old == new)
changed = len(old_and_new_columns) - unchanged
print(f"\nSummary:")
print(f"Total columns: {len(old_and_new_columns)}")
print(f"Changed columns: {changed}")
print(f"Unchanged columns: {unchanged}")

pg_ordered_matrix = pg_matrix_sorted.drop('9-9', axis=1)
print(pg_ordered_matrix.columns)

pg_ordered_matrix.to_excel(processed_dir / "1_pg_chivas_renamed_cols.xlsx", index=False)