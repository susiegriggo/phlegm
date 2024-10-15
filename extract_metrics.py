""" 
Compute metrics for a directory of structures with different numbers of subunits for the same oligomeric state 
""" 

# imports 
from src import metrics 
import click 
import glob 
import pickle 
import json
import numpy as np 
import re 

@click.command() 
@click.option('--input', '-i', type=click.Path(exists=True), help='Path to the input file or directory')
@click.option('--output', '-o', type=click.Path(), help='Name of output pickle file containing metrics')
@click.option('--max_subunits', '-m', type=int, help='Maximum number of subunits in a complex')
@click.option('--num_structures', '-n', type=int, help='Number of structures to consider for deach subunit size', default=3)

def main(input, output, max_subunits, num_structures): 

    # read the files in the input
    files = glob.glob(input + '/*')

    # get the prefix from the directory name 
    prefix = input + re.split('/', input)[-2]

    # dictionary to store the metrics 
    sequence_metrics = dict() 

    for s in range(2, max_subunits+1):   

        subunit_metrics = dict() 

        for i in range(1, num_structures+1): 

            # creat dicionary to sore this replicate 
            complex_metrics = dict() 

            # get the pdb files corresponding to this entry 
            s_p = prefix + '.' + str(s) + '_unrelaxed_rank_00' + str(i)
            s_j = prefix + '.' + str(s) + '_scores_rank_00' + str(i)
        
            # get the json file for this entry
            j = [i for i in files if s_j in i and 'json' in i ][0] 
            f = json.load(open(j, 'r'))
            p = [i for i in files if s_p in i and 'pdb' in i][0]
        

            # extract the metrics from the json file 
            complex_metrics['iptm'] = f.get('iptm')
            complex_metrics['ptm'] = f.get('ptm')
            complex_metrics['max_pae'] = f.get('max_pae')
            complex_metrics['mean_plddt'] = np.mean(f.get('plddt'))
            complex_metrics['median_plddt'] = np.median(f.get('plddt'))

            # get the pdockq2 score - one score per chain - read why 
            pdockq2, avgif_pae, plddt_lst = metrics.calc_pdockq2(j,p)
            complex_metrics['pdockq2'] = pdockq2.to_list()
            complex_metrics['avgif_pae'] = avgif_pae
            complex_metrics['plddt_lst'] = plddt_lst

            # compute the number of contacts 
            chain_coords, chain_plddt = metrics.read_pdb(p)
            for a in range(1,10): # loop over different thresholds in angstrongs 
                complex_metrics['n_contacts_' + str(a) + 'A'] = metrics.calc_interface_contacts(chain_coords, t = a)

            # update the sequence metrics dictionary 
            subunit_metrics['rep_' + str(i)] = complex_metrics

        # update with the number of subunits 
        sequence_metrics[str(s) + '_subunits'] = subunit_metrics 
              
    # save this to a pickle file 
    pickle.dump(sequence_metrics, open(output + '.pkl', 'wb'))
    print(sequence_metrics)

if __name__ == '__main__': 
    main()