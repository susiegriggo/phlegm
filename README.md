# phlegm
Phage Homomer Likelihood Estimator and Generator Method 


1. split fasta file into separate fasta files
bash src/split_fasta.sh test_data/test.faa test_data_fasta

2. Colabsearch to build MSA for each sequence
mkdir test_data_search
cd test_data_search
qsub -v DIR=../test_data_fasta ../src/collab_search_dir.sh

3. Generate an alignment for each number of subunits to test


4. Colabbatch

5. Run metrics and state prediction
