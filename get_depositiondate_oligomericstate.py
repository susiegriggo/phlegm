import pandas as pd
import requests
import time
import logging
import pickle
import re
from collections import Counter

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Load and process pharokka data
#pharokka = pd.read_csv('/home/grig0076/scratch/phlegm/all_pdb/pharokka/pharokka_proteins_full_merged_output.tsv', sep='\t')
#pharokka_phrog = pharokka[pharokka['phrog'] != 'No_PHROGs_HMM']
#pharokka_phrog['entry'] = pharokka_phrog['ID'].apply(lambda x: re.split('_', x)[0])
#pharokka_phrog = pharokka_phrog[pharokka_phrog.index.isin(pharokka_phrog.drop('ID', axis=1).drop_duplicates().index)]

# Identify homomeric complexes
#entry_chain_counts = Counter(pharokka_phrog['entry'])
#heteromers = {entry for entry, count in entry_chain_counts.items() if count > 1}
#homomeric_entries = set(entry_chain_counts.keys()) - heteromers

#mmseqs_phrogs = pd.read_csv('/home/grig0076/scratch/phlegm/all_pdb/mmseqs/results.tsv', sep='\t', header=None)
#mmseqs_phrogs.columns = ['phrog', 'pdb_id','alnScore', 'seqIdentity', 'eVal', 'qStart','qEnd','qLen','tStart','tEnd','tLen'] 
#mmseqs_phrogs['pdb'] = [re.split('_', p)[0] for p  in mmseqs_phrogs['pdb_id']]  

# Load the m8 file
df = pd.read_csv("/home/grig0076/scratch/phlegm/PHROGs/phrog_represenative_pdb_seqres_minseqid0.3_c0.4.m8", sep="\t", header=None)
df.columns = [
    "query", "target", "pident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bitscore"
]
df['pdb'] = [re.split('_', p)[0] for p  in df['target']]  


# List of PDB IDs
#pdb_ids = list(set(mmseqs_phrogs['pdb']))
pdb_ids = list(set(df['pdb']))
print('Number of pdb ids to get: ' + str(len(pdb_ids)))

def fetch_pdb_info(pdb_id):
    """Fetch deposition date and oligomeric states for all assemblies from RCSB PDB API."""
    base_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    assembly_base_url = f"https://data.rcsb.org/rest/v1/core/assembly/{pdb_id}"
    
    try:
        logging.info(f"Fetching data for PDB ID: {pdb_id}")
        
        # Fetch main entry data
        response = requests.get(base_url)
        response.raise_for_status()
        data = response.json()
        
        deposition_date = data.get('rcsb_accession_info', {}).get("deposit_date", "Unknown")
        available_assemblies = data.get("rcsb_entry_container_identifiers", {}).get("assembly_ids", [])

        assemblies = {}

        for assembly_id in available_assemblies:
            assembly_url = f"{assembly_base_url}/{assembly_id}"
            assembly_response = requests.get(assembly_url)
            if assembly_response.status_code == 200:
                assembly_data = assembly_response.json()

                # Extract oligomeric state
                oligomeric_states = []
                if "rcsb_struct_symmetry" in assembly_data:
                    oligomeric_states = [sym.get('oligomeric_state', "Unknown") for sym in assembly_data["rcsb_struct_symmetry"]]
                
                # Fallback to `pdbx_struct_assembly.oligomeric_details` if no symmetry info is found
                if not oligomeric_states:
                    oligomeric_states = [assembly_data.get("pdbx_struct_assembly", {}).get("oligomeric_details", "Unknown")]

                assemblies[assembly_id] = oligomeric_states
            
        logging.info(f"Fetched data for PDB ID: {pdb_id} - Deposition Date: {deposition_date}, Assemblies: {assemblies}")
        return {"PDB_ID": pdb_id, "Deposition_Date": deposition_date, "Assemblies": assemblies}
    
    except requests.exceptions.RequestException as e:
        logging.error(f"Request error for PDB ID: {pdb_id} - {e}")
        return {"PDB_ID": pdb_id, "Deposition_Date": "Error", "Assemblies": {}}
    except Exception as e:
        logging.error(f"Unexpected error for PDB ID: {pdb_id} - {e}")
        return {"PDB_ID": pdb_id, "Deposition_Date": "Error", "Assemblies": {}}

# Fetch data for all PDB IDs
results = {}
batch_size = 1000
output_file = "/home/grig0076/scratch/phlegm/all_pdb/pharokka/pdb_oligomeric_states_mmseqs2_2025-04-08.pkl"

for i, pdb_id in enumerate(pdb_ids, start=1):
    results[pdb_id] = fetch_pdb_info(pdb_id)
    
    if i % batch_size == 0 or i == len(pdb_ids):
        logging.info(f"[{i}/{len(pdb_ids)}] Processed {i} sequences...")
        
        # Save results to pickle file
        with open(output_file, 'wb') as f:
            pickle.dump(results, f)
        
        # Comment out the line that clears the results dictionary
        # results.clear()
    
    time.sleep(1.0)

logging.info(f"Processing complete. Results saved to {output_file}")
