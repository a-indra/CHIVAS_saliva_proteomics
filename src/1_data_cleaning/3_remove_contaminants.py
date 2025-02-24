### [1] SET-UP WORKSPACE ###
import pandas as pd
from pathlib import Path

# Set project root directory
project_root = Path(__file__).parent.parent.parent

# Define input/output paths
processed_dir = project_root / "data" / "processed"

# Load in data
chivas_protein = pd.read_excel(processed_dir / "2_pg_chivas_with_protein_names.xlsx")

### [2] REMOVE KERATIN  ###
# Print the initial number of rows
print(f"Initial number of rows: {len(chivas_protein)}")

# Remove rows containing 'keratin' (case-insensitive) in the Protein.Names column
chivas_protein_filtered = chivas_protein[~chivas_protein['Protein.Names'].str.contains('keratin', case=False, na=False)]

# Print the number of rows after filtering
print(f"Number of rows after removing keratin: {len(chivas_protein_filtered)}")

# Print the number of rows removed
print(f"Number of rows removed: {len(chivas_protein) - len(chivas_protein_filtered)}")

# Export to Excel
chivas_protein_filtered.to_excel(processed_dir / "3_pg_chivas_no_keratin.xlsx", index=False)