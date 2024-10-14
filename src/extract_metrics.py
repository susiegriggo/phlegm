#!/usr/bin/env python 

"""
Extract metrics for a directory of directories of collabfold output 
"""

from sysconfig import get_makefile_filename
from webbrowser import get
from src import get_metrics
import glob
import json
import numpy as np
import re
import click 
import pickle


@click.command()
@click.option('--input', '-i', type=click.Path(exists=True), help='Path to the input file or directory')
@click.option('--output', '-o', type=click.Path(), help='Path to the output file or directory')


def main(input, output):

    # create a dictionary to store information 
    colabfold_metrics_per_sequence = dict() 

    # path of the files
    test_complex_paths = glob.glob(input + '/*')
    #test_complex_paths = [i for i in test_complex_paths if '.2' in i]
    
    # loop through each of the structures 
    for t in test_complex_paths: 

        # need to handle the three different rank models separateley 
        max_rank = 3 

        # files for the current complex to investigate 
        files = glob.glob(t + '/*')
        
        # create dictionary to store models 
        sequence_metrics = dict() 

        # loop through each of the models for this sequence
        for i in range(1,max_rank + 1): 
            
            # create info to store this particular replicate 
            complex_metrics = dict() 

            # get the pdb file and json corresponding to this entry 
            s_p = 'unrelaxed_rank_00' + str(i) 
            s_j = '_rank_00' + str(i) 
            p = [i for i in files if s_p in i]

            if len(p) > 0: 
                p = p[0] 
                j = [j for j in [i for i in files if s_j in i] if 'json' in j][0]

                # compute the metrics 
                complex_metrics = dict() 
                f = json.load(open(j, 'r'))

                # store the relevant metrics that we would like to anlayse 
                complex_metrics['iptm'] = f.get('iptm')
                complex_metrics['ptm'] = f.get('ptm')
                complex_metrics['max_pae'] = f.get('max_pae')
                complex_metrics['mean_plddt'] = np.mean(f.get('plddt'))
                complex_metrics['median_plddt'] = np.median(f.get('plddt'))

                # see how i can extract which scores are at the interface so I can get pae/plddt at the interface 

                # get the pdockq2 score - one score per chain - read why 
                pdockq2, avgif_pae, plddt_lst = get_metrics.calc_pdockq2(j,p)
                complex_metrics['pdockq2'] = pdockq2.to_list()
                complex_metrics['avgif_pae'] = avgif_pae
                complex_metrics['plddt_lst'] = plddt_lst

                # compute the number of contacts 
                chain_coords, chain_plddt = get_metrics.read_pdb(p)
                for a in range(1,10): # loop over different thresholds for computing the number of contacts 
                    complex_metrics['n_contacts_' + str(a) + 'A'] = get_metrics.calc_interface_contacts(chain_coords, t = a)
            
                # update the sequence metrics dictionary 
                sequence_metrics['rep_' + str(i)] = complex_metrics
                print(sequence_metrics)

        # update the overall dictionary 
        colabfold_metrics_per_sequence[re.split('/', t)[-1]] = sequence_metrics
            
    # save the dictionary to a pickle file 
    print('saving to: ' + output + '.pkl', flush=True)
    pickle.dump(colabfold_metrics_per_sequence, open(output + '.pkl', 'wb'))

if __name__ == '__main__':
    main()
