import pandas as pd
import pickle
import json
import numpy as np
import argparse
from pathlib import Path
from step3_visualization import analyze_tree_distribution
from step3_analyze import *
from step3_optimize import *
import os

def process_bootstrap_data(
    patient: str,
    bootstrap_list: list[int],
    bootstrap_parent_dir: Path,
    output_dir: Path,
    method: str = 'phylowgs'
) -> None:
    """
    Process bootstrap data for a given patient from a specific directory containing bootstrap replicates.
    
    Args:
        patient (str): Patient ID (primarily for naming output files).
        bootstrap_list (list[int]): List of bootstrap numbers to process.
        bootstrap_parent_dir (Path): The direct path to the directory containing bootstrapN subdirectories 
                                     (e.g., /path/to/data/CRUK0044/initial/bootstraps).
        output_dir (Path): The directory where aggregation results will be saved.
        method (str): Method used (default: phylowgs).
    """
    # Use the provided output directory instead of hardcoding it
    aggregation_output_dir = output_dir
    
    print(f"\nProcessing bootstrap data for patient: {patient}")
    print(f"Bootstrap parent directory: {bootstrap_parent_dir}")
    print(f"Looking for {len(bootstrap_list)} bootstrap replicates.")
    print(f"Aggregation results will be saved to: {aggregation_output_dir}")
    
    # Create aggregation output directory if it doesn't exist
    aggregation_output_dir.mkdir(exist_ok=True, parents=True)
    
    # Initialize storage dictionaries
    tree_distribution = {
        'cp_tree': [], 'node_dict': [], 'node_dict_name': [], 
        'node_dict_re': [], 'tree_structure': [], 'freq': [], 
        'clonal_freq': [], 'vaf_frac': []
    }
    tree_aggregation = {
        'cp_tree': [], 'node_dict': [], 'node_dict_name': [], 
        'node_dict_re': [], 'tree_structure': [], 'freq': [], 
        'clonal_freq': [], 'vaf_frac': []
    }
    
    processed_bootstraps = 0
    for bootstrap_idx in bootstrap_list:
        print(f"\nProcessing bootstrap {bootstrap_idx}")
        
        # Construct paths directly within the bootstrap_parent_dir
        current_bootstrap_replicate_dir = bootstrap_parent_dir / f"bootstrap{bootstrap_idx}"
        summ_file = current_bootstrap_replicate_dir / "result.summ.json.gz"
        muts_file = current_bootstrap_replicate_dir / "result.muts.json.gz"
        mutass_file = current_bootstrap_replicate_dir / "result.mutass.zip"
        
        if not current_bootstrap_replicate_dir.exists():
            print(f"Directory not found: {current_bootstrap_replicate_dir}")
            continue

        required_files = [summ_file, muts_file, mutass_file]
        missing_files = [f for f in required_files if not f.exists()]
        if missing_files:
            print(f"Skipping bootstrap {bootstrap_idx} in {current_bootstrap_replicate_dir} - missing files:")
            for f in missing_files:
                print(f"  - {f}")
            continue
        
        processed_bootstraps += 1
        tree_structure, node_dict, node_dict_name, node_dict_re, final_tree_cp, prev_mat, clonal_freq, vaf_frac = process_phylowgs_output(
            summ_file, muts_file, mutass_file
        )
        
        tree_distribution = combine_tree(
            node_dict, node_dict_name, node_dict_re, tree_structure, 
            final_tree_cp, clonal_freq, vaf_frac, method, tree_distribution
        )
        
        tree_aggregation['tree_structure'].append(tree_structure)
        tree_aggregation['cp_tree'].append(final_tree_cp)
        tree_aggregation['node_dict'].append(node_dict)
        tree_aggregation['node_dict_re'].append(node_dict_re)
        tree_aggregation['node_dict_name'].append(node_dict_name)
        tree_aggregation['freq'].append(1)
        tree_aggregation['clonal_freq'].append(clonal_freq)
        tree_aggregation['vaf_frac'].append(vaf_frac)
    
    print(f"\nSuccessfully processed {processed_bootstraps} out of {len(bootstrap_list)} bootstrap replicates.")
    
    if processed_bootstraps == 0:
        print(f"No bootstrap data successfully processed for patient {patient} in {bootstrap_parent_dir}")
        print("This could be because required bootstrap subdirectories or their result files were missing.")
        return

    # Analyze and save results to the aggregation_output_dir
    analysis_type_for_output_name = 'initial'
    analyze_tree_distribution(tree_distribution, aggregation_output_dir, patient, analysis_type_for_output_name, fig=True)
    
    best_bootstrap_idx = np.argmax(tree_distribution['freq'])
    results_dict = {
        'node_dict_name': tree_distribution['node_dict'][best_bootstrap_idx],
        'tree_structure': tree_distribution['tree_structure'][best_bootstrap_idx]
    }
    
    with open(aggregation_output_dir / f"{patient}_results_bootstrap_{analysis_type_for_output_name}_best.json", 'w') as f:
        json.dump(results_dict, f)
    
    with open(aggregation_output_dir / f'{method}_bootstrap_summary.pkl', 'wb') as g:
        pickle.dump(tree_distribution, g)
    
    with open(aggregation_output_dir / f'{method}_bootstrap_aggregation.pkl', 'wb') as g:
        pickle.dump(tree_aggregation, g)

    print(f"Aggregation results saved in {aggregation_output_dir}")

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Process and aggregate bootstrap data for phylogenetic analysis.')
    
    parser.add_argument('patient', type=str,
                      help='Patient ID (used for naming output files).')
    
    parser.add_argument('--bootstrap-list', type=int, nargs='+',
                      help='List of bootstrap numbers to process (e.g., 1 2 3 ... 100).')
    
    parser.add_argument('--bootstrap-parent-dir', type=str, required=True,
                      help='The direct path to the directory containing bootstrapN subdirectories (e.g., /path/to/data/CRUK0044/initial/bootstraps).')

    parser.add_argument('--output-dir', type=str, required=True,
                      help='Output directory for aggregation results.')

    parser.add_argument('--method', type=str, default='phylowgs',
                      help='Method used (default: phylowgs).')
        
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    
    bootstrap_parent_dir = Path(args.bootstrap_parent_dir)
    output_dir = Path(args.output_dir)
    
    if not args.bootstrap_list:
        print("Error: --bootstrap-list must be provided.")
        exit(1)

    process_bootstrap_data(
        patient=args.patient,
        bootstrap_list=args.bootstrap_list,
        bootstrap_parent_dir=bootstrap_parent_dir,
        output_dir=output_dir,
        method=args.method
    )

"""
Example usage with explicit output directory:
python 3-aggregation/aggregate.py CRUK0044 \
    --bootstrap-list $(seq 1 100) \
    --bootstrap-parent-dir /path/to/data/CRUK0044/initial/bootstraps \
    --output-dir /path/to/data/CRUK0044/initial/aggregation_results \
    --method phylowgs
"""