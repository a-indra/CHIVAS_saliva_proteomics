# This script connects to UniProt API to add protein names to each entry

### [1] SET-UP WORKSPACE ###
# Import modules
import requests
import pandas as pd
from tqdm import tqdm
import time
import re
from pathlib import Path

# Set project root directory
project_root = Path(__file__).parent.parent.parent

# Define input/output paths
processed_dir = project_root / "data" / "processed"
input_file = processed_dir / "1_pg_chivas_renamed_cols.xlsx"

# Load in data
df = pd.read_excel(input_file)

### [2] DEFINE FUNCTION FOR RETRIEVING NAME ###
def get_protein_info(protein_id):
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        full_name = data['proteinDescription']['recommendedName']['fullName']['value'] if 'recommendedName' in data['proteinDescription'] else data['proteinDescription']['submissionNames'][0]['fullName']['value']
        short_name = data['uniProtkbId']
        return full_name, short_name
    else:
        return "Not found", "Not found"

# Create new columns
df['Protein.Names'] = ''
df['Protein.Abbrev'] = ''

# Process each unique protein ID
unique_proteins = df['Protein.Group'].unique()

for protein_id in tqdm(unique_proteins, desc="Processing proteins"):
    try:
        full_name, short_name = get_protein_info(protein_id)
        mask = df['Protein.Group'] == protein_id
        df.loc[mask, 'Protein.Names'] = full_name
        df.loc[mask, 'Protein.Abbrev'] = short_name
    except Exception as e:
        print(f"Error processing {protein_id}: {str(e)}")
    time.sleep(0.1)  # To avoid overloading the API

# Export to excel
df.to_excel(processed_dir / "2_pg_chivas_with_protein_names.xlsx", index=False)