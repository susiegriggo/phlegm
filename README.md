[![DOI](https://zenodo.org/badge/548652990.svg)](https://doi.org/10.5281/zenodo.8128898)
# phlegm

**P**hage **H**omomer **L**evel **E**stimate and **G**eneration **M**ethod 

📄 Preprint: [Computational prediction resolves thousands of homooligomeric phage protein structures](https://doi.org/10.64898/2026.05.24.727406) (bioRxiv, 2026)
<p align="center">
  <img src="https://github.com/susiegriggo/phlegm/blob/main/phlegm.png" width="600" title="phlegm logo" alt="phlegm logo">
</p> 

## Quick start

### 1. Install

Requires [conda/mamba](https://github.com/conda-forge/miniforge).

```bash
git clone https://github.com/susiegriggo/phlegm.git
cd phlegm
mamba env create -f environment.yml    # use -n NAME or -p /path to choose where it goes
conda activate phlegm
bash install.sh
```

Use `mamba` (included with [Miniforge](https://github.com/conda-forge/miniforge)) rather than `conda` to create the environment. Older conda versions (before 23.10) use a slow solver that can take 30+ minutes on this environment, while mamba takes a few minutes. With conda ≥ 23.10, `conda env create -f environment.yml` is equally fast.

`environment.yml` contains ColabFold v1.5.5, MMseqs2 15.6f452 and all Python dependencies (the versions used in the manuscript). You can create the environment however you like, or add these dependencies to an existing environment with `mamba env update -f environment.yml`.

`install.sh` then sets up PHLEGM inside whichever conda environment is active. It installs the `phlegm` command and, on Linux, GPU-enabled JAX (CUDA 12 via pip; needs NVIDIA driver ≥ 525 to use a GPU). GPU JAX falls back to the CPU when there's no GPU, so on a cluster you can install from the login node and use the same environment on CPU and GPU nodes. Use `bash install.sh --cpu` to skip it (saves ~3 GB if you'll never run `fold`). On macOS, JAX is CPU-only.

Nothing large is downloaded at install time:
- **AlphaFold-Multimer weights (~4 GB)** are downloaded automatically the first time you run `phlegm fold` (or `phlegm run`). On HPC, GPU nodes often have no internet access, so run `phlegm download-weights` once on the login node. The job scripts in `examples/` do this for you. Weights are stored in ColabFold's cache (`~/.cache/colabfold` on Linux); use `--weights-dir` or `$PHLEGM_WEIGHTS_DIR` to keep one shared copy elsewhere, e.g. for a whole lab.
- **ColabFold sequence databases** are not downloaded at all (see [Using a local database](#using-a-local-colabfold-database)).

### Test your installation

`test_data/phrog_807.faa` contains one protein: PHROG 807, a major tail protein, which the manuscript uses as a worked example. Its closest PDB homolog (8rk8) is a hexamer.

```bash
phlegm run -i test_data/phrog_807.faa -o phlegm_test
```

This models 2- to 10-mers with 5 models each (the 10-mer is ~2,600 residues), so expect it to take roughly 6–12 hours on a single A40/A100 (mostly the 8- to 10-mers) and to need a GPU with plenty of memory (≥ 40 GB is safest for the largest assemblies). It uses the public ColabFold server for the MSA, so it needs internet access. `phlegm_test/phlegm_predictions.tsv` should report `predicted_subunits` = 6 (`predicted_state` = `5-6-mer`), with a max mean ipTM of around 0.75.

For a quicker check that everything runs (not a full prediction), limit the stoichiometries and models:

```bash
phlegm run -i test_data/phrog_807.faa -o phlegm_test --max-oligomer 6 --num-models 2
```

### 2. Look up a phage genome in the precomputed resource (no GPU)

The PHLEGM resource holds predictions for the representative sequences of 22,554 PHROGs (`data/phlegm_phrog_predictions.csv.gz`, also the manuscript's supplementary table). If your proteins already have PHROG assignments, for example from [pharokka](https://github.com/gbouras13/pharokka) or [phold](https://github.com/gbouras13/phold), `phlegm lookup` gives each protein the prediction of its PHROG in seconds. It needs no GPU, MSAs or AlphaFold weights:

```bash
phlegm lookup -i phold_out/phold.gbk -o lookup_out        # or pharokka_out/pharokka.gbk
```

The GenBank output of phold or pharokka is the simplest input: it holds each CDS's PHROG (`/phrog`) and protein sequence (`/translation`). `-i` also accepts phold's `*_per_cds_predictions.tsv`, pharokka's `*_cds_final_merged_output.tsv`, or any TSV/CSV whose first two columns are a protein ID and a PHROG (`phrog_123` or `123`). With table input, add the protein sequences with `--proteins` (optional), e.g. `--proteins pharokka_out/prodigal-gv.faa`. Given the sequences, PHLEGM reports how similar each protein is to its PHROG's representative and writes `to_fold.faa`, the proteins without a lookup result, ready for step 3. Use unmasked sequences (the GenBank translations or pharokka's `prodigal-gv.faa`, not phold's `phold_aa.fasta`): residues masked as `X` are ignored when computing identity.

Results are written to `lookup_out/phlegm_lookup.tsv`:

| column | description |
|---|---|
| `protein_id`, `phrog` | protein and its PHROG assignment |
| `status` | `predicted`, `not_in_resource` (the PHROG representative was > 350 aa or its monomer model had pLDDT ≤ 70, so it was not modelled) or `no_phrog` |
| `predicted_state` / `predicted_subunits` / `max_mean_iptm` | the PHLEGM prediction for the PHROG representative |
| `phrog_annotation`, `phrog_category` | PHROG functional annotation |
| `protein_length_aa`, `representative_length_aa` | lengths of the protein and of the PHROG representative |
| `identity_to_representative`, `coverage_of_protein` | local alignment identity (%) and the percentage of the protein covered (needs sequences: GenBank input or `--proteins`) |
| `structure_url` | the PHLEGM structure of the representative: the predicted oligomer, or the monomer for monomer calls |
| `note` | why a protein has no prediction, or a warning if it is only distantly related to the representative |

A protein inherits the prediction made for its PHROG's representative, so the call is most reliable when the protein is similar to the representative. PHROGs assigned by structure (phold) can be remote homologs with low sequence identity, so treat those calls with more caution, or fold the protein itself with `phlegm run`. Proteins not in the resource can be predicted directly:

```bash
phlegm run -i lookup_out/to_fold.faa -o phlegm_out
```

### 3. Run the whole pipeline from a FASTA file

```bash
phlegm run -i proteins.fasta -o phlegm_out
```

This runs the whole pipeline for every sequence in the FASTA file. It builds one MSA per protein and reuses it for each stoichiometry, predicts the 2-mer to 10-mer complexes with AlphaFold-Multimer, and calls the oligomeric state from the mean ipTM. By default MSAs come from the public ColabFold MMseqs2 server, so no database is needed. This suits a few to tens of proteins on a single GPU machine.

Results are written to `phlegm_out/phlegm_predictions.tsv`:

| column | description |
|---|---|
| `protein_id` | first word of the FASTA header (file-name safe) |
| `predicted_subunits` / `predicted_state` | predicted number of subunits and binned state (Monomer, 2-mer, 3-mer, 4-mer, 5-6-mer, 7+-mer) |
| `max_mean_iptm` | highest mean ipTM across stoichiometries |
| `top_multimer_subunits` | stoichiometry with the highest mean ipTM (before applying the threshold) |
| `mean_iptm_<n>mer` | mean ipTM across models for each stoichiometry |
| `status` | `ok`, `incomplete`, `not folded`, or the reason a protein failed |

### 4. Running the steps separately (HPC)

The pipeline has three steps with different hardware needs, and each is its own command:

| command | what it does | runs on | needs |
|---|---|---|---|
| `phlegm msa -i proteins.fasta -o out [--db DB]` | one MSA per protein | CPU node (high memory with `--db`) | ColabFold database, or internet for the server |
| `phlegm fold -o out [--shard I/N]` | predicts every stoichiometry with AlphaFold-Multimer | GPU node | AlphaFold weights only |
| `phlegm predict -o out` | calls oligomeric states, writes the TSV | anywhere | — |

`phlegm run` simply runs all three. All commands share the output directory, and any command can be re-run to resume. Before starting, `fold` checks three things and stops with a clear message if any fails:
- every MSA exists (it never builds MSAs itself, so it won't quietly fall back to the public server);
- JAX can see a GPU, rather than running AlphaFold on the CPU for days (override with `--allow-cpu`);
- the AlphaFold weights are present. Without `--shard`, it downloads them if they're missing. With `--shard`, it never downloads, so parallel jobs can't all download into the same folder at once; run `phlegm download-weights` first.

`--shard I/N` makes `fold` process only shard `I` of `N`, so the proteins can be split across a GPU job array. Ready-made submission scripts chain the three steps with job dependencies:

```bash
# SLURM: MSA job -> GPU job array (20 shards) -> predict job
bash examples/submit_slurm.sh proteins.fasta phlegm_out /path/to/colabfold_db 20

# PBS Pro (e.g. NCI Gadi, which has no job arrays): MSA job -> 20 GPU jobs -> predict job
bash examples/submit_pbs.sh proteins.fasta phlegm_out /path/to/colabfold_db 20
```

Edit the resource settings at the top of each script (account/project, queues, memory, walltime) for your cluster.

### Using a local ColabFold database

For more than a few dozen proteins, use a local database instead of the public server. PHLEGM does not install it: it is large (~1–2 TB with indices) and is often shared on clusters already. To build one, follow ColabFold's [`setup_databases.sh`](https://github.com/sokrypton/ColabFold/blob/v1.5.5/setup_databases.sh). Then pass `--db /path/to/colabfold_db` to `phlegm msa` or `phlegm run`.

### Options

`phlegm <command> --help` lists every option for a command. The main ones:

| option | command(s) | purpose |
|---|---|---|
| `--db DIR` | msa, run | local ColabFold database for `colabfold_search` |
| `--msa-dir DIR` | msa, run | use precomputed monomer MSAs named `<id>.a3m` |
| `--threads N` | msa, run | threads for `colabfold_search` |
| `--min-oligomer` / `--max-oligomer` | fold, predict, run | stoichiometries to model (default 2–10) |
| `--num-models` / `--num-recycle` | fold, run | AlphaFold-Multimer settings (default 5 models, 3 recycles, as in the manuscript) |
| `--colabfold-args "..."` | fold, run | extra arguments for `colabfold_batch` |
| `--shard I/N` | fold | process one shard of the proteins (job arrays) |
| `--weights-dir DIR` | fold, run, download-weights | AlphaFold weights location (default: ColabFold's cache, or `$PHLEGM_WEIGHTS_DIR`) |
| `--keep-going` | fold, run | continue with other proteins if one fails |
| `-t/--threshold` | predict, run | mean ipTM threshold for calling a multimer (default 0.65) |
| `--extended-metrics` | predict, run | also compute pDockQ2, interface contacts and PAE with `extract_metrics.py` |
| `--dry-run` | msa, fold, run | print commands without running them |
| `--proteins FASTA` / `--resource TABLE` | lookup | protein sequences for table input (GenBank input already has them); an alternative per-PHROG table |

Output layout:
```
phlegm_out/
├── phlegm_predictions.tsv     # one row per protein
├── phlegm_predictions.pkl     # full per-protein results (all ipTM scores)
├── input/                     # validated sequences used by every step
├── logs/                      # one log per command (and per fold shard)
├── msas/<id>.a3m              # monomer MSA (generated once per protein)
├── oligomer_msas/<id>/<id>.<n>.a3m
└── structures/<id>/           # colabfold_batch outputs for every stoichiometry
```

## How-to (individual steps)

The steps below are what the `phlegm` command runs. They are documented for users who want to run or customize each step by hand, for example as separate jobs on an HPC scheduler.

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

If you use PHLEGM in your research, please cite:

```
Grigson, S.R., Geliashvili, N., Schubert, T., Bouras, G., Mallawaarachchi, V., Bogacz, M.,
Hellmich, U.A., Edwards, R.A., Dutilh, B.E. (2026) Computational prediction resolves
thousands of homooligomeric phage protein structures. bioRxiv.
https://doi.org/10.64898/2026.05.24.727406
```

### Acknowledgments

Thankyou to Laura Inglis for designing the PHELGM logo! 

This project includes code adapted from [FoldDock](https://gitlab.com/ElofssonLab/FoldDock), originally licensed under the Apache License, Version 2.0.

Copyright 2021 Patrick Bryant, Gabriele Pozzati and Arne Elofsson

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at:

   http://www.apache.org/licenses/LICENSE-2.0
