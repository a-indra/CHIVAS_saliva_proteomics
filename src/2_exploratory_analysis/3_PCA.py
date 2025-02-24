# This script transforms data into a format amenable to PCA grouped by time-point.
# Then, PCA is executed and plotted 

### [1] SET UP WORK SPACE ###
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set project root directory and paths
project_root = Path(__file__).parent.parent.parent
input_dir = project_root / "data" / "processed" 
figures_dir = project_root / "data" / "figures" / "exploratory"
results_dir = project_root /"data" / "processed" / "results"

# Read the data
original_df = pd.read_excel(input_dir/"4_pg_chivas_filtered.xlsx")

### [2] RESHAPE DATA ###
# Get protein columns and sample columns
protein_metadata = ['First.Protein.Description', 'Genes', 'Protein.Group', 'Protein.Ids', 'Protein.Names','Protein.Abbrev']
sample_columns = [col for col in original_df.columns if col not in protein_metadata]

# Transform data
rows = []
for column in sample_columns:
    patient_id, time_point = column.split('_', 1)
    
    # Create row with Patient_ID and Time_Point
    row_dict = {
        'Patient_ID': int(patient_id),
        'Time_Point': time_point
    }
    
    # Add protein values
    protein_values = original_df[column].values
    for protein_name, value in zip(original_df['Protein.Names'], protein_values):
        row_dict[protein_name] = value
    
    rows.append(row_dict)

# Convert to DataFrame and sort
transformed_df = pd.DataFrame(rows)
transformed_df = transformed_df.sort_values(['Patient_ID', 'Time_Point'])

# Exclude specific timepoints 
df = transformed_df[~transformed_df['Time_Point'].isin(['72h', '96h', '3mth'])]

### [3] DEAL WITH MISSING VALUES ### 
# Temporarily drop 'metadata' to allow specific handling of protein abundance values
metadata = df[['Patient_ID', 'Time_Point']]
protein_data = df.drop(['Patient_ID', 'Time_Point'], axis=1)

# Remove proteins with too many missing values (above 50%)
missing_percentage = protein_data.isnull().mean()
proteins_to_keep = missing_percentage[missing_percentage < 0.5].index
protein_data = protein_data[proteins_to_keep]
print(f"\nNumber of proteins after removing those with >50% missing values: {protein_data.shape[1]}")

# Replace remaining missing values with 0
protein_data_filled = protein_data.fillna(0)

# Aggregate data by time point (using mean)
aggregated_data = protein_data_filled.groupby(metadata['Time_Point']).mean()

### [4] STANDARDISE DATA AND EXECUTE PCA ###
# Standardise the aggregated data
scaler = StandardScaler()
aggregated_data_scaled = pd.DataFrame(
    scaler.fit_transform(aggregated_data), 
    columns=aggregated_data.columns, 
    index=aggregated_data.index
)

# Perform PCA
pca = PCA()
pca_result = pca.fit_transform(aggregated_data_scaled)
explained_variance_ratio = pca.explained_variance_ratio_

# Create DataFrame with PCA results
pca_df = pd.DataFrame(
    data=pca_result[:, :2], 
    columns=['PC1', 'PC2'], 
    index=aggregated_data.index
)

### [5] CREATE AND SAVE PLOTS ###
# Create PCA plot
plt.figure(figsize=(12, 8))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', s=100)

# Add labels for each point
for i, txt in enumerate(pca_df.index):
    plt.annotate(txt, (pca_df.iloc[i, 0], pca_df.iloc[i, 1]), 
                xytext=(5, 5), textcoords='offset points')

# Add plot labels
plt.title('PCA of Proteomics Data by Timepoint')
plt.xlabel(f'PC1 ({explained_variance_ratio[0]:.2%} explained variance)')
plt.ylabel(f'PC2 ({explained_variance_ratio[1]:.2%} explained variance)')
plt.tight_layout()

# Save and display plot
plt.savefig(figures_dir / 'PCA_plot_by_timepoint.svg', format='svg', dpi=300, bbox_inches='tight')
plt.show()

# Print results
print(f"\nTotal explained variance by PC1 and PC2: {sum(explained_variance_ratio[:2]):.2%}")
print("\nPCA coordinates for each time point:")
print(pca_df)