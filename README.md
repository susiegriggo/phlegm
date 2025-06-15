# phlegm
**P**hage **H**omomer **L**ikelihood **E**stimator and **G**enerator **M**ethod 
<p align="center">
  <img src="https://github.com/susiegriggo/phlegm/blob/main/phlegm.png" width="600" title="phlegm logo" alt="phlegm logo">
</p> 

## How-to
1. split fasta file into separate fasta files
bash src/split_fasta.sh test_data/test.faa test_data_fasta

2. Colabsearch to build MSA for each sequence
mkdir test_data_search
cd test_data_search
qsub -v DIR=../test_data_fasta ../src/collab_search_dir.sh

3. Generate an alignment for each number of subunits to test
bash src/generate_alignments.sh test_data_search/  test_data_alignments

4. Colabbatch

5. Run metrics

6. Compare metrics to make state prediction for each sequence 

### Acknowledgments

Thankyou to Laura Inglis for designing the PHELGM logo! 

This project includes code adapted from [FoldDock](https://gitlab.com/ElofssonLab/FoldDock), originally licensed under the Apache License, Version 2.0.

Copyright 2021 Patrick Bryant, Gabriele Pozzati and Arne Elofsson

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at:

   http://www.apache.org/licenses/LICENSE-2.0
