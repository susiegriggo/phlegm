# imports 
import pandas as pd
from Bio import SeqIO
import re
import requests
from collections import Counter
from tqdm import tqdm

df = pd.read_csv("/home/grig0076/GitHubs/phlegm/files/phrog_represenative_pdb_seqres_minseqid0.3_c0.7.m8", sep="\t", header=None)
df.columns = [
    "query", "target", "pident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bitscore"
]

# Load query and target lengths from FASTA files
def load_lengths(fasta_file):
    return {record.id: len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

query_lengths = load_lengths("/home/grig0076/GitHubs/phlegm/files/nonsingleton_representative_sequences.fasta") # PHROG representative sequences
target_lengths = load_lengths("/home/grig0076/GitHubs/phlegm/files/pdb_seqres.cleaned.txt") # Sequences in the protein databank 

# Map lengths to dataframe
df["qlen"] = df["query"].map(query_lengths)
df["tlen"] = df["target"].map(target_lengths)

# Compute coverage
df["query_coverage"] = (df["qend"] - df["qstart"] + 1) / df["qlen"]
df["target_coverage"] = (df["tend"] - df["tstart"] + 1) / df["tlen"]
df['pdb'] = [re.split('_',i)[0] for i in df['target']]

# list of pdb ids to extract info for 
pdb_ids = list(set(df['pdb'].to_list()))

def get_experimental_method(pdb_id):
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        method = data.get("exptl", [{}])[0].get("method", "Not available")
        return method
    else:
        return f"Error {response.status_code}"

# Fetch methods
results = []
for pdb_id in tqdm(pdb_ids, desc="Fetching experimental methods"):
    method = get_experimental_method(pdb_id)
    results.append({'PDB_ID': pdb_id.upper(), 'Experimental_Method': method})

# Convert to DataFrame
method_df = pd.DataFrame(results)
print(method_df)

# Optional: Save to CSV
method_df.to_csv("pdb_experimental_methods.csv", index=False)
