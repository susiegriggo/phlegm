# phlegm
**P**hage **H**omomer **L**evel **E**stimator and **G**enerator **M**ethod 
<p align="center">
  <img src="https://github.com/susiegriggo/phlegm/blob/main/phlegm.png" width="600" title="phlegm logo" alt="phlegm logo">
</p> 

## How-to

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
   
   This step creates 9 different alignment files for each protein, named like:
   ```
   protein1_2.a3m  # For 2-mer prediction
   protein1_3.a3m  # For 3-mer prediction
   ...
   protein1_10.a3m # For 10-mer prediction
   ```

3. **Run AlphaFold multimer prediction**
   ```bash
   # Create output directory
   mkdir -p alphafold_output
   
   # Run AlphaFold on each multimer configuration separately
   colabfold_batch --model-type alphafold2_multimer_v3 \
                   --num-recycle 3 \
                   --num-models 5 \
                   --templates false \
                   --amber false \
                   multimer_alignments/protein1_2.a3m \
                   alphafold_output/protein1_2
   
   # Repeat for each subunit configuration (3-mer, 4-mer, etc.)
   colabfold_batch --model-type alphafold2_multimer_v3 \
                   --num-recycle 3 \
                   --num-models 5 \
                   --templates false \
                   --amber false \
                   multimer_alignments/protein1_3.a3m \
                   alphafold_output/protein1_3
   
   # For a full directory of alignments (processing each protein and subunit count)
   for protein in multimer_alignments/*.a3m; do
     base=$(basename $protein .a3m)
     
     # Create a dedicated output directory
     mkdir -p alphafold_output/$base
     
     colabfold_batch --model-type alphafold2_multimer_v3 \
                     --num-recycle 3 \
                     --num-models 5 \
                     --templates false \
                     --amber false \
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
