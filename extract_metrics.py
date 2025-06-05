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
import tarfile
import os
import tempfile
import shutil

def extract_tar_gz(input_path):
    if tarfile.is_tarfile(input_path):
        with tarfile.open(input_path, "r:gz") as tar:
            temp_dir = tempfile.mkdtemp()  # Create a temporary directory
            tar.extractall(temp_dir)  # Extract the contents
            return temp_dir
    else:
        raise ValueError(f"{input_path} is not a valid .tar.gz file")

def detect_prefix(files):
    """
    Automatically detect the prefix used in the files.
    
    Returns:
        tuple: (detected_prefix, is_consistent)
    """
    # Regular expression to match the pattern: something.{digits}_scores_rank or something.{digits}_unrelaxed_rank
    pattern = r"(.+?)\.(\d+)_(?:unrelaxed|scores)_rank_\d+"
    
    prefixes = {}
    
    for file in files:
        basename = os.path.basename(file)
        match = re.search(pattern, basename)
        
        if match:
            prefix = match.group(1)
            if prefix in prefixes:
                prefixes[prefix] += 1
            else:
                prefixes[prefix] = 1
    
    # If we found any prefixes
    if prefixes:
        # Find the most common prefix
        most_common_prefix = max(prefixes.items(), key=lambda x: x[1])[0]
        
        # Check if all files use the same prefix
        is_consistent = len(prefixes) == 1
        
        return most_common_prefix, is_consistent, prefixes
    
    return None, False, {}

def determine_subunit_info(files, prefix):
    """
    Determine the maximum number of subunits and number of structures
    based on the files in the input directory.
    
    Returns:
        tuple: (max_subunits, num_structures_dict, missing_subunits)
    """
    subunit_structures = {}
    
    for file in files:
        basename = os.path.basename(file)
        
        # Match patterns like: prefix.{subunit}_unrelaxed_rank_00{number} or prefix.{subunit}_scores_rank_00{number}
        pattern = rf"{re.escape(prefix)}\.(\d+)_(?:unrelaxed|scores)_rank_00(\d+)"
        match = re.search(pattern, basename)
        
        if match:
            try:
                subunit_count = int(match.group(1))
                structure_num = int(match.group(2))
                
                if subunit_count not in subunit_structures:
                    subunit_structures[subunit_count] = set()
                
                subunit_structures[subunit_count].add(structure_num)
            except (ValueError, IndexError):
                continue
    
    # Check for missing subunit counts
    if not subunit_structures:
        return 0, {}, []
    
    max_subunits = max(subunit_structures.keys())
    expected_subunits = set(range(2, max_subunits + 1))  # Starting from 2
    actual_subunits = set(subunit_structures.keys())
    missing_subunits = sorted(expected_subunits - actual_subunits)
    
    # Determine max structures per subunit count
    num_structures_dict = {s: max(nums) + 1 for s, nums in subunit_structures.items()}
    
    return max_subunits, num_structures_dict, missing_subunits

@click.command() 
@click.option('--input', '-i', type=click.Path(exists=True), help='Path to the input file or directory')
@click.option('--output', '-o', type=click.Path(), help='Name of output pickle file containing metrics')
@click.option('--prefix', '-p', type=str, help='prefix used in colab output (auto-detected if not provided)')
def main(input, output, prefix): 

    # Check if input is a .tar.gz file and extract it
    temp_dir = None
    if input.endswith('.tar.gz'):
        try:
            temp_dir = extract_tar_gz(input)  # Extract the tar.gz file
            input = temp_dir
        except ValueError as e:
            print(f"Error: {str(e)}")
            return
    
    try:
        # read the files in the input
        files = glob.glob(input + '/*')
        
        # Auto-detect prefix if not provided
        if not prefix:
            detected_prefix, is_consistent, prefix_counts = detect_prefix(files)
            if not detected_prefix:
                raise ValueError("Could not auto-detect prefix from files. Please provide a prefix using the -p option.")
            
            if not is_consistent:
                prefixes_info = ', '.join([f"'{p}' ({c} files)" for p, c in prefix_counts.items()])
                print(f"Warning: Multiple prefixes detected: {prefixes_info}")
                print(f"Using most common prefix: '{detected_prefix}'")
            else:
                print(f"Detected prefix: '{detected_prefix}'")
            
            prefix = detected_prefix
        else:
            # Verify if the provided prefix matches files
            _, _, prefix_counts = detect_prefix(files)
            if prefix not in prefix_counts:
                raise ValueError(f"Provided prefix '{prefix}' does not match any files in the directory.")
        
        # Determine max_subunits, num_structures, and missing subunits dynamically
        max_subunits, num_structures_dict, missing_subunits = determine_subunit_info(files, prefix)
        
        if max_subunits == 0:
            raise ValueError(f"No valid structure files found matching prefix '{prefix}'. Please check prefix and input directory.")
        
        print(f"Detected maximum subunits: {max_subunits}")
        print(f"Detected structures per subunit count: {num_structures_dict}")
        
        if missing_subunits:
            print(f"Warning: Missing subunit counts: {missing_subunits}. These will be skipped.")
        
        # dictionary to store the metrics 
        sequence_metrics = dict() 

        for s in range(2, max_subunits+1):   
            print(s) 
            
            # Skip if no structures for this subunit count
            if s not in num_structures_dict:
                print(f"No structures found for subunit count {s}, skipping...")
                continue
                
            subunit_metrics = dict() 
            
            # Get number of structures for this subunit count
            num_structures = num_structures_dict[s]

            for i in range(1, num_structures): 
                print(i) 
                # creat dicionary to sore this replicate 
                complex_metrics = dict() 

                # get the pdb files corresponding to this entry 
                s_p = prefix + '.' + str(s) + '_unrelaxed_rank_00' + str(i)
                s_j = prefix + '.' + str(s) + '_scores_rank_00' + str(i)
            
                # get the json file for this entry
                j = [i for i in files if s_j in i and 'json' in i ]
                if not j:
                    print(f"Warning: No JSON file found for {s_j}, skipping...")
                    continue
                j = j[0]
                
                f = json.load(open(j, 'r'))
                
                p = [i for i in files if s_p in i and 'pdb' in i]
                if not p:
                    print(f"Warning: No PDB file found for {s_p}, skipping...")
                    continue
                p = p[0]
            
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
        print(f"Metrics saved successfully to {output}.pkl")
    
    except ValueError as e:
        print(f"Error: {str(e)}")
    finally:
        # Clean up temporary directory if created
        if temp_dir and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

if __name__ == '__main__': 
    main()
