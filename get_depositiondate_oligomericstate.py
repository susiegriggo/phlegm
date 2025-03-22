# imports 
import pandas as pd
import re
from collections import Counter
import requests
import time
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# get in the pharokka dataframe
pharokka = pd.read_csv('/home/grig0076/scratch/phlegm/all_pdb/pharokka/pharokka_proteins_full_merged_output.tsv', sep = '\t')
pharokka_phrog = pharokka[pharokka['phrog'] != 'No_PHROGs_HMM']
pharokka_phrog['entry'] = [re.split('_',i)[0] for i in pharokka_phrog['ID']] 

# drop the duplicates rows 
pharokka_phrog = pharokka_phrog[pharokka_phrog.index.isin(pharokka_phrog.drop('ID', axis=1).drop_duplicates().index)] 

# get the counts
entry_chain_counts = pd.DataFrame.from_dict(Counter(pharokka_phrog['entry']), orient = 'index').sort_values(0,axis=0)

# determine which entries are heteromeric complexes 
heteromers = entry_chain_counts[entry_chain_counts[0] > 1].index.to_list()

# just get the homomeric complexes
pharokka_phrogs_homomers = pharokka_phrog[~pharokka_phrog['entry'].isin(heteromers)]

# List of PDB IDs
pdb_ids = pharokka_phrogs_homomers['entry'].to_list()

def fetch_pdb_info(pdb_id):
    """Fetch deposition date and oligomeric state from RCSB PDB API with error handling."""
    base_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    assembly_base_url = f"https://data.rcsb.org/rest/v1/core/assembly/{pdb_id}"

    try:
        logging.info(f"Fetching data for PDB ID: {pdb_id}")

        # Fetch main entry data
        response = requests.get(base_url)
        response.raise_for_status()
        data = response.json()

        # Get deposition date
        deposition_date = data.get('rcsb_accession_info', {}).get("deposit_date", "Unknown")

        # Get available assemblies
        available_assemblies = data.get("rcsb_entry_container_identifiers", {}).get("assembly_ids", [])
        
        # Set the defaults
        oligomeric_state = "Unknown"
        symmetry = "No symmetry data available"

        # Loop through assemblies to find oligomeric state
        for assembly_id in available_assemblies:
            assembly_url = f"{assembly_base_url}/{assembly_id}"
            assembly_response = requests.get(assembly_url)
            
            if assembly_response.status_code == 200:
                assembly_data = assembly_response.json()
                
                # Extract symmetry info
                if "rcsb_struct_symmetry" in assembly_data:
                    symmetry_info = assembly_data["rcsb_struct_symmetry"]
                    if isinstance(symmetry_info, list) and symmetry_info:
                        symmetry = ", ".join(
                            [f"{sym.get('type', 'Unknown')} ({sym.get('symbol', 'N/A')})"
                             for sym in symmetry_info if "type" in sym and "symbol" in sym]
                        )

                # Extract oligomeric state
                oligomeric_state = assembly_data.get("pdbx_struct_assembly", {}).get("oligomeric_details", "Unknown")
                break  # Stop after first valid assembly

        logging.info(f"Fetched data for PDB ID: {pdb_id} - Deposition Date: {deposition_date}, Oligomeric State: {oligomeric_state}")
        return pdb_id, deposition_date, oligomeric_state

    except requests.exceptions.RequestException as e:
        logging.error(f"Request error for PDB ID: {pdb_id} - {e}")
        return pdb_id, "Error", "Error"
    except Exception as e:
        logging.error(f"Unexpected error for PDB ID: {pdb_id} - {e}")
        return pdb_id, "Error", "Error"

# Fetch data for all PDB IDs
results = []
batch_size = 1000  # Show progress every 1000 entries
output_file = "/home/grig0076/scratch/phlegm/all_pdb/pharokka/pdb_oligomeric_states.tsv"

# Initialize the output file
df = pd.DataFrame(columns=["PDB_ID", "Deposition_Date", "Oligomeric_State"])
df.to_csv(output_file, index=False, sep='\t')

# loop through the pdb ids 
for i, pdb_id in enumerate(pdb_ids, start=1):
    results.append(fetch_pdb_info(pdb_id))
    
    # Save progress every 1000 sequences
    if i % batch_size == 0 or i == len(pdb_ids):
        logging.info(f"[{i}/{len(pdb_ids)}] Processed {i} sequences...")
        
        # Convert results to DataFrame and append to the file
        batch_df = pd.DataFrame(results, columns=["PDB_ID", "Deposition_Date", "Oligomeric_State"])
        batch_df.to_csv(output_file, mode='a', header=False, index=False, sep='\t')
        
        # Clear results list
        results = []

    time.sleep(0.5)  # Small delay to avoid overloading the server

logging.info(f"\nProcessing complete. Results saved to {output_file}")
