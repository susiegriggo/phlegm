[![DOI](https://zenodo.org/badge/548652990.svg)](https://doi.org/10.5281/zenodo.8128898)
# phlegm

**P**hage **H**omomer **L**evel **E**stimate and **G**eneration **M**ethod 
<p align="center">
  <img src="https://github.com/susiegriggo/phlegm/blob/main/phlegm.png" width="600" title="phlegm logo" alt="phlegm logo">
</p> 

## Quick start

### 1. Install (one command)

Requires [conda/mamba](https://github.com/conda-forge/miniforge).

```bash
git clone https://github.com/susiegriggo/phlegm.git
cd phlegm
bash install.sh
conda activate phlegm
```

`install.sh` creates a `phlegm` conda environment with ColabFold v1.5.5, MMseqs2 15.6f452 and all Python dependencies (the versions used in the manuscript). If it finds an NVIDIA GPU, it also installs GPU-enabled JAX (CUDA 12 via pip, which needs NVIDIA driver ≥ 525). It then downloads the AlphaFold-Multimer v3 weights. Options: `-n <env name>`, `--cpu`, `--no-weights`. To create the environment without the extras: `conda env create -f environment.yml`.

### 2. Run the whole pipeline from a FASTA file

```bash
python phlegm.py -i proteins.fasta -o phlegm_out
```

This runs every step below automatically for every sequence in the FASTA file. It builds one MSA per protein and reuses it for each stoichiometry, predicts the 2-mer to 10-mer complexes with AlphaFold-Multimer, and calls the oligomeric state from the mean ipTM. Results are written to `phlegm_out/phlegm_predictions.tsv`:

| column | description |
|---|---|
| `protein_id` | first word of the FASTA header (file-name safe) |
| `predicted_subunits` / `predicted_state` | predicted number of subunits and binned state (Monomer, 2-mer, 3-mer, 4-mer, 5-6-mer, 7+-mer) |
| `max_mean_iptm` | highest mean ipTM across stoichiometries |
| `top_multimer_subunits` | stoichiometry with the highest mean ipTM (before applying the threshold) |
| `mean_iptm_<n>mer` | mean ipTM across models for each stoichiometry |
| `status` | `ok`, `incomplete`, or the reason a protein failed |

Useful options (`python phlegm.py --help` lists all):

| option | purpose |
|---|---|
| `--db /path/to/colabfold_db` | build MSAs with `colabfold_search` against a local ColabFold database instead of the public MMseqs2 server (recommended for large datasets) |
| `--msa-dir DIR` | use precomputed monomer MSAs named `<id>.a3m` |
| `--min-oligomer` / `--max-oligomer` | range of stoichiometries to model (default 2–10) |
| `-t/--threshold` | mean ipTM threshold for calling a multimer (default 0.65) |
| `--colabfold-args "..."` | pass extra arguments to `colabfold_batch` |
| `--extended-metrics` | also compute pDockQ2, interface contacts and PAE with `extract_metrics.py` |
| `--stop-after msa` | stop after one step, e.g. build MSAs on a CPU node and continue on a GPU node |
| `--keep-going` | keep going with the remaining proteins if one fails |
| `--dry-run` | print the commands without running them |

Runs can be resumed: re-running the same command skips finished MSAs and structures.

Output layout:
```
phlegm_out/
├── phlegm_predictions.tsv     # one row per protein
├── phlegm_predictions.pkl     # full per-protein results (all ipTM scores)
├── phlegm.log
├── msas/<id>.a3m              # monomer MSA (generated once per protein)
├── oligomer_msas/<id>/<id>.<n>.a3m
└── structures/<id>/           # colabfold_batch outputs for every stoichiometry
```

## How-to (individual steps)

The steps below are what `phlegm.py` runs. They are documented for users who want to run or customize each step by hand, for example as separate jobs on an HPC scheduler.

### Workflow Overview
The PHLEGM method predicts oligomeric states from protein sequences using [ColabFold](https://github.com/sokrypton/ColabFold) for AlphaFold structure prediction and ipTM scoring. If you use this method, please also cite ColabFold and AlphaFold Multimer:

```
Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S, Steinegger M. ColabFold: Making protein folding accessible to all.
Nature Methods (2022) doi: 10.1038/s41592-022-01488-1

Evans R, O'Neill M, Pritzel A, et al. Protein complex prediction with AlphaFold-Multimer.
bioRxiv (2022) doi: 10.1101/2021.10.04.463034
```

The workflow consists of four main steps:

1. **Generate MSAs using colabfold_search**
   ```bash
   # For a single protein sequence
   colabfold_search protein.fasta /path/to/database protein_msas
   ```

2. **Generate alignments for different subunit stoichiometries**
   ```bash
   # Create MSAs for 2-mer through 10-mer configurations (9 alignments per protein)
   # Takes a directory of MSAs and generates modified MSAs for multimeric prediction
   bash src/generate_alignments.sh msas_output/ multimer_alignments
   ```
   
   This step creates one alignment file per stoichiometry for each protein, named like:
   ```
   protein1.2.a3m  # For 2-mer prediction
   protein1.3.a3m  # For 3-mer prediction
   ...
   protein1.10.a3m # For 10-mer prediction
   ```

3. **Run AlphaFold multimer prediction**
   ```bash
   # Create output directory
   mkdir -p alphafold_output
   
   # Run AlphaFold on each multimer configuration separately
   colabfold_batch --model-type alphafold2_multimer_v3 \
                   --num-recycle 3 \
                   --num-models 5 \
                   multimer_alignments/protein1.2.a3m \
                   alphafold_output/protein1.2
   
   # Repeat for each subunit configuration (3-mer, 4-mer, etc.)
   colabfold_batch --model-type alphafold2_multimer_v3 \
                   --num-recycle 3 \
                   --num-models 5 \
                   multimer_alignments/protein1.3.a3m \
                   alphafold_output/protein1.3
   
   # For a full directory of alignments (processing each protein and subunit count)
   for protein in multimer_alignments/*.a3m; do
     base=$(basename $protein .a3m)
     
     # Create a dedicated output directory
     mkdir -p alphafold_output/$base
     
     colabfold_batch --model-type alphafold2_multimer_v3 \
                     --num-recycle 3 \
                     --num-models 5 \
                     $protein \
                     alphafold_output/$base
   done
   
   # After all predictions are complete, collect all outputs into one directory for analysis
   mkdir -p alphafold_all_outputs
   find alphafold_output -name "*.pdb" -o -name "*.json" | xargs -I{} cp {} alphafold_all_outputs/
   ```

   ### Expected Output Structure
   After running AlphaFold, your combined output directory should contain files like:
   ```
   protein_name.2_unrelaxed_rank_001.pdb  # Structure for 2-mer
   protein_name.2_scores_rank_001.json    # Scores for 2-mer
   protein_name.3_unrelaxed_rank_001.pdb  # Structure for 3-mer
   protein_name.3_scores_rank_001.json    # Scores for 3-mer
   ...
   protein_name.10_unrelaxed_rank_001.pdb # Structure for 10-mer
   protein_name.10_scores_rank_001.json   # Scores for 10-mer
   ```

4. **Analyze results and predict oligomeric state**

   You can choose between two analysis approaches depending on your needs:
   
   **Option A: Fast prediction with `predict_oligomeric_state.py`**
   
   Use this when you need quick oligomeric state predictions:
   - Extracts ipTM scores only
   - Applies the validated 0.65 threshold
   - Provides simple binned classifications (Monomer, 2-mer, 3-mer, 4-mer, 5-6-mer, 7+-mer)
   - Ideal for high-throughput screening or automated pipelines
   
   ```bash
   # Basic usage (run on the directory with all protein outputs combined)
   python predict_oligomeric_state.py -i alphafold_all_outputs -o predictions.pkl
   
   # With custom threshold
   python predict_oligomeric_state.py -i alphafold_all_outputs -t 0.7 -o predictions.pkl
   ```
   
   **Option B: Comprehensive analysis with `extract_metrics.py`**
   
   Use this for detailed structural analysis:
   - Extracts ipTM, pTM, pLDDT scores
   - Calculates pDockQ2 scores per chain
   - Measures interface contacts at multiple distance thresholds
   - Computes average interface PAE values
   - Ideal for research, benchmarking, or generating ML datasets
   
   ```bash
   python extract_metrics.py -i alphafold_all_outputs -o metrics_output
   ```

### Dependencies

All dependencies for the pipeline are installed by `bash install.sh` (see [Quick start](#quick-start)). The scripts require the following Python packages:

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
https://doi.org/10.5281/zenodo.8128898
```

### Acknowledgments

Thankyou to Laura Inglis for designing the PHELGM logo! 

This project includes code adapted from [FoldDock](https://gitlab.com/ElofssonLab/FoldDock), originally licensed under the Apache License, Version 2.0.

Copyright 2021 Patrick Bryant, Gabriele Pozzati and Arne Elofsson

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at:

   http://www.apache.org/licenses/LICENSE-2.0
