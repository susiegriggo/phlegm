# phlegm
**P**hage **H**omomer **L**ikelihood **E**stimator and **G**enerator **M**ethod 
<p align="center">
  <img src="https://github.com/susiegriggo/phlegm/blob/main/phlegm.png" width="600" title="phlegm logo" alt="phlegm logo">
</p> 

## How-to

### Workflow Overview
The PHLEGM method predicts oligomeric states from protein sequences using AlphaFold structure prediction and ipTM scoring:

1. **Prepare input sequences**
   ```bash
   # Split multi-FASTA file into individual FASTA files (one per protein)
   bash src/split_fasta.sh test_data/test.faa test_data_fasta
   ```

2. **Generate multiple sequence alignments (MSAs)**
   ```bash
   # Use ColabSearch to build MSAs for each sequence
   mkdir test_data_search
   cd test_data_search
   qsub -v DIR=../test_data_fasta ../src/collab_search_dir.sh
   ```

3. **Generate alignments for different subunit stoichiometries**
   ```bash
   # Create MSAs for 2-mer, 3-mer, 4-mer, etc. configurations
   bash src/generate_alignments.sh test_data_search/ test_data_alignments
   ```

4. **Run AlphaFold structure prediction**
   ```bash
   # Use ColabBatch to predict structures for different subunit counts
   # Configure for 5 models and 3 recycles per subunit configuration
   # This generates structures for each protein at different oligomeric states
   # Output: protein.2_unrelaxed_rank_001.pdb, protein.3_unrelaxed_rank_001.pdb, etc.
   ```

5. **Extract comprehensive metrics** (optional - for detailed analysis)
   ```bash
   python extract_metrics.py -i path/to/alphafold/output -o output_metrics
   ```

6. **Predict oligomeric state**
   ```bash
   # Run prediction on directory containing all AlphaFold outputs
   python predict_oligomeric_state.py -i path/to/alphafold/output -o predictions.pkl
   ```

### Output Structure Expected
After ColabBatch, your output directory should contain files like:
```
protein_name.2_unrelaxed_rank_001.pdb
protein_name.2_scores_rank_001.json
protein_name.3_unrelaxed_rank_001.pdb  
protein_name.3_scores_rank_001.json
protein_name.4_unrelaxed_rank_001.pdb
protein_name.4_scores_rank_001.json
...
```

The `predict_oligomeric_state.py` script will:
- Auto-detect the protein name prefix
- Extract ipTM scores from the JSON files
- Calculate mean ipTM for each subunit count
- Apply the 0.65 threshold to determine oligomeric state
- Output binned predictions (Monomer, 2-mer, 3-mer, 4-mer, 5-6-mer, 7+-mer)

### Scripts

#### extract_metrics.py
Extracts comprehensive metrics from AlphaFold Colab output including ipTM, pTM, pDockQ2, and interface contacts.

Usage:
```bash
python extract_metrics.py -i /path/to/alphafold/output -o metrics_output
python extract_metrics.py -i results.tar.gz -o metrics_output -p protein_name
```

#### predict_oligomeric_state.py
Predicts oligomeric state based on ipTM scores using the approach from the predict_phrogs.ipynb notebook. Applies the 0.65 ipTM threshold and bins results into oligomeric state categories (Monomer, 2-mer, 3-mer, 4-mer, 5-6-mer, 7+-mer).

Usage:
```bash
# Basic usage with auto-detected prefix
python predict_oligomeric_state.py -i /path/to/alphafold/output

# With tar.gz input and custom output
python predict_oligomeric_state.py -i results.tar.gz -o predictions.pkl

# With custom prefix and threshold
python predict_oligomeric_state.py -i /path/to/output -p my_protein -t 0.7 -o results.pkl
```

Options:
- `-i, --input`: Path to directory or tar.gz file containing AlphaFold output
- `-o, --output`: Optional output pickle file to save results
- `-p, --prefix`: Prefix used in colab output (auto-detected if not provided)
- `-t, --threshold`: ipTM threshold for confident prediction (default: 0.65)
- `-m, --min-structures`: Minimum number of structures needed for prediction (default: 2)

Output contains:
- `predicted_state_bin`: Binned oligomeric state (e.g., "2-mer", "5-6-mer")
- `max_iptm`: Highest mean ipTM score found
- `max_iptm_subunits`: Subunit count with highest mean ipTM
- `per_subunits_mean_iptm`: Mean ipTM scores for each subunit count tested

### Dependencies

The scripts require the following Python packages:

#### For extract_metrics.py:
- `click` - Command line interface
- `numpy` - Numerical computations
- `pandas` - Data manipulation
- `biopython` - PDB file parsing
- Standard library: `json`, `pickle`, `glob`, `re`, `tarfile`, `os`, `tempfile`, `shutil`, `collections`

#### For predict_oligomeric_state.py:
- `click` - Command line interface  
- `numpy` - Numerical computations
- Standard library: `json`, `pickle`, `glob`, `re`, `tarfile`, `os`, `tempfile`, `shutil`, `collections`

#### For notebooks (predict_phrogs.ipynb):
- `pandas` - Data manipulation
- `numpy` - Numerical computations
- `matplotlib` - Plotting
- `seaborn` - Statistical plotting
- `scikit-learn` - Machine learning (TSNE)
- `umap-learn` - UMAP dimensionality reduction
- `tqdm` - Progress bars
- `wesanderson` - Color palettes
- Standard library: `pickle`, `glob`, `re`, `collections`

Install dependencies:
```bash
pip install click numpy pandas biopython matplotlib seaborn scikit-learn umap-learn tqdm wesanderson
```

## Citation

If you use PHLEGM in your research, please cite this repository:

```
Grigson, S., Bouras, G., Dutilh, B.E., Edwards, R.A. PHLEGM: Phage Homomer Level Estimator and Generator Method. 
https://github.com/susiegriggo/phlegm
```

### Acknowledgments

Thankyou to Laura Inglis for designing the PHELGM logo! 

This project includes code adapted from [FoldDock](https://gitlab.com/ElofssonLab/FoldDock), originally licensed under the Apache License, Version 2.0.

Copyright 2021 Patrick Bryant, Gabriele Pozzati and Arne Elofsson

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at:

   http://www.apache.org/licenses/LICENSE-2.0
