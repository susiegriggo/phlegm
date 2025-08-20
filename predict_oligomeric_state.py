"""
Predict oligomeric state from AlphaFold output based on ipTM scores.
Uses the approach from predict_phrogs.ipynb with ipTM > 0.65 threshold.
"""

import click
import glob
import json
import numpy as np
import re
import tarfile
import os
import tempfile
import shutil
from collections import defaultdict
import pickle

def extract_tar_gz(input_path):
    """Extract tar.gz file to temporary directory."""
    if tarfile.is_tarfile(input_path):
        with tarfile.open(input_path, "r:gz") as tar:
            temp_dir = tempfile.mkdtemp()
            tar.extractall(temp_dir)
            return temp_dir
    else:
        raise ValueError(f"{input_path} is not a valid .tar.gz file")

def detect_prefix(files):
    """Automatically detect the prefix used in the files."""
    pattern = r"(.+?)\.(\d+)_(?:unrelaxed|scores)_rank_\d+"
    prefixes = {}
    
    for file in files:
        basename = os.path.basename(file)
        match = re.search(pattern, basename)
        
        if match:
            prefix = match.group(1)
            prefixes[prefix] = prefixes.get(prefix, 0) + 1
    
    if prefixes:
        most_common_prefix = max(prefixes.items(), key=lambda x: x[1])[0]
        is_consistent = len(prefixes) == 1
        return most_common_prefix, is_consistent, prefixes
    
    return None, False, {}

def determine_subunit_info(files, prefix):
    """Determine available subunit counts and structures."""
    subunit_structures = {}
    
    for file in files:
        basename = os.path.basename(file)
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
    
    if not subunit_structures:
        return 0, {}
    
    max_subunits = max(subunit_structures.keys())
    num_structures_dict = {s: max(nums) + 1 for s, nums in subunit_structures.items()}
    
    return max_subunits, num_structures_dict

def extract_iptm_scores(files, prefix, max_subunits, num_structures_dict):
    """Extract ipTM scores for all available structures."""
    iptm_data = defaultdict(list)
    
    for s in range(2, max_subunits + 1):
        if s not in num_structures_dict:
            continue
            
        num_structures = num_structures_dict[s]
        
        for i in range(1, num_structures):
            s_j = prefix + '.' + str(s) + '_scores_rank_00' + str(i)
            
            # Find JSON file
            json_files = [f for f in files if s_j in f and 'json' in f]
            if not json_files:
                continue
                
            try:
                with open(json_files[0], 'r') as f:
                    data = json.load(f)
                    iptm = data.get('iptm')
                    if iptm is not None:
                        iptm_data[s].append(iptm)
            except (json.JSONDecodeError, IOError) as e:
                print(f"Warning: Could not read {json_files[0]}: {e}")
                continue
    
    return dict(iptm_data)

def predict_oligomeric_state(iptm_data, threshold=0.65, min_structures=2):
    """
    Predict oligomeric state based on ipTM scores using the approach from predict_phrogs.ipynb.
    
    Args:
        iptm_data: Dict mapping subunit count to list of ipTM scores
        threshold: Minimum mean ipTM threshold for confident prediction
        min_structures: Minimum number of structures needed for prediction
    
    Returns:
        Dict with prediction results
    """
    # Dictionary for mapping to a homo-oligomeric state (from notebook)
    subunits_state_dict = {
        1: 'Monomer', 
        2: '2-mer',  
        3: '3-mer', 
        4: '4-mer', 
        5: '5-6-mer', 
        6: '5-6-mer',
        7: '7+-mer',
        8: '7+-mer',
        9: '7+-mer',
        10: '7+-mer',
    }
    
    results = {
        'predicted_state': None,
        'predicted_state_bin': None,
        'max_iptm': 0,
        'max_iptm_subunits': 1,
        'all_scores': iptm_data
    }
    
    # Calculate mean ipTM for each subunit count (like in notebook)
    per_subunits_iptm = {}
    for subunit_count, scores in iptm_data.items():
        if len(scores) >= min_structures:
            per_subunits_iptm[subunit_count] = np.mean(scores)
    
    if not per_subunits_iptm:
        # No valid predictions possible
        results['predicted_state_bin'] = 'Monomer'
        return results
    
    # Find the maximum mean ipTM and corresponding subunit count
    max_iptm = max(per_subunits_iptm.values())
    max_iptm_subunits = max(per_subunits_iptm.keys(), key=per_subunits_iptm.get)
    
    # Apply threshold logic from notebook
    if max_iptm < threshold:
        max_iptm_subunits = 1
    
    results['max_iptm'] = max_iptm
    results['max_iptm_subunits'] = max_iptm_subunits
    results['predicted_state'] = max_iptm_subunits
    results['predicted_state_bin'] = subunits_state_dict.get(max_iptm_subunits, 'Monomer')
    results['per_subunits_mean_iptm'] = per_subunits_iptm
    
    return results

def format_output(results, prefix):
    """Format prediction results for display."""
    print(f"\n=== Oligomeric State Prediction for {prefix} ===")
    print(f"Predicted State: {results['predicted_state_bin']}")
    print(f"Max mean ipTM: {results['max_iptm']:.3f}")
    print(f"Best subunit count: {results['max_iptm_subunits']}")
    
    if results['per_subunits_mean_iptm']:
        print(f"\nMean ipTM by subunit count:")
        for subunit_count, mean_iptm in sorted(results['per_subunits_mean_iptm'].items()):
            print(f"  {subunit_count}-mer: {mean_iptm:.3f}")

@click.command()
@click.option('--input', '-i', type=click.Path(exists=True), 
              help='Path to the input file or directory containing AlphaFold output')
@click.option('--output', '-o', type=click.Path(), 
              help='Optional output pickle file to save results')
@click.option('--prefix', '-p', type=str, 
              help='Prefix used in colab output (auto-detected if not provided)')
@click.option('--threshold', '-t', type=float, default=0.65,
              help='ipTM threshold for confident prediction (default: 0.65)')
@click.option('--min-structures', '-m', type=int, default=2,
              help='Minimum number of structures needed for prediction (default: 2)')
def main(input, output, prefix, threshold, min_structures):
    """Predict oligomeric state from AlphaFold output based on ipTM scores."""
    
    temp_dir = None
    
    try:
        # Handle tar.gz input
        if input.endswith('.tar.gz'):
            temp_dir = extract_tar_gz(input)
            input = temp_dir
        
        # Get files
        files = glob.glob(os.path.join(input, '*'))
        
        # Auto-detect or verify prefix
        if not prefix:
            detected_prefix, is_consistent, prefix_counts = detect_prefix(files)
            if not detected_prefix:
                raise ValueError("Could not auto-detect prefix. Please provide using -p option.")
            
            if not is_consistent:
                prefixes_info = ', '.join([f"'{p}' ({c} files)" for p, c in prefix_counts.items()])
                print(f"Warning: Multiple prefixes detected: {prefixes_info}")
                print(f"Using most common prefix: '{detected_prefix}'")
            
            prefix = detected_prefix
        else:
            _, _, prefix_counts = detect_prefix(files)
            if prefix not in prefix_counts:
                raise ValueError(f"Provided prefix '{prefix}' does not match any files.")
        
        print(f"Using prefix: '{prefix}'")
        
        # Determine available subunit counts and structures
        max_subunits, num_structures_dict = determine_subunit_info(files, prefix)
        
        if max_subunits == 0:
            raise ValueError(f"No valid structure files found matching prefix '{prefix}'.")
        
        print(f"Found structures for subunit counts: {sorted(num_structures_dict.keys())}")
        
        # Extract ipTM scores
        iptm_data = extract_iptm_scores(files, prefix, max_subunits, num_structures_dict)
        
        if not iptm_data:
            raise ValueError("No ipTM scores could be extracted from the files.")
        
        # Predict oligomeric state
        results = predict_oligomeric_state(iptm_data, threshold, min_structures)
        
        # Display results
        format_output(results, prefix)
        
        # Save results if requested
        if output:
            with open(output, 'wb') as f:
                pickle.dump(results, f)
            print(f"\nResults saved to {output}")
    
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")
    finally:
        # Clean up temporary directory
        if temp_dir and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

if __name__ == '__main__':
    main()
main()
