#!/usr/bin/env python3
"""
Multi-sample compatible marker selection script.

This script works with the new multi-sample DataFrame format while maintaining
compatibility with the original marker selection algorithms and tree distribution data.
"""

from step4_optimize import *
from step4_optimize_fraction import *
from step4_convert_ssm import convert_ssm_to_dataframe_multi
import pandas as pd
import pickle
import argparse
import matplotlib.pyplot as plt
import os
import sys


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Run marker selection analysis with multi-sample SSM input.')
    
    parser.add_argument('patient', type=str,
                      help='Patient ID')
    
    parser.add_argument('-r', '--read-depth', type=int, default=1500,
                      help='Read depth for analysis (default: 1500)')
    
    parser.add_argument('-a', '--aggregation-dir', type=str, required=True,
                        help='Path to the directory containing phylowgs_bootstrap_aggregation.pkl')
    
    parser.add_argument('-s', '--ssm-file', type=str, required=True,
                        help='Path to the ssm.txt file for this patient')
    
    parser.add_argument('-o', '--output-dir', type=str,
                        help='Path to output directory for marker selection results')
    
    # VAF filtering options
    parser.add_argument('-f', '--filter-strategy', type=str, 
                       choices=['any_high', 'all_high', 'majority_high', 'specific_samples'],
                       default='any_high',
                       help='VAF filtering strategy (default: any_high)')
    
    parser.add_argument('-t', '--filter-threshold', type=float, default=0.9,
                       help='VAF threshold for filtering (default: 0.9)')
    
    parser.add_argument('--filter-samples', type=int, nargs='+',
                       help='Sample indices for specific_samples filtering strategy')
    
    return parser.parse_args()


def create_backward_compatible_dataframe(multi_sample_df):
    """
    Create a backward-compatible DataFrame that mimics the old cf/st format.
    
    This function takes a multi-sample DataFrame and creates the two-column
    VAF format expected by the original marker selection code.
    
    Args:
        multi_sample_df: DataFrame with VAF_sample_0, VAF_sample_1, etc. columns
        
    Returns:
        DataFrame with Variant_Frequencies_cf and Variant_Frequencies_st columns
    """
    # Get VAF columns
    vaf_columns = [col for col in multi_sample_df.columns if col.startswith('VAF_sample_')]
    
    if len(vaf_columns) < 2:
        print("Warning: Less than 2 samples available. Using available samples for cf/st.")
        
    # Create backward-compatible DataFrame
    compat_df = multi_sample_df.copy()
    
    # Use first two samples as cf and st equivalents
    if len(vaf_columns) >= 1:
        compat_df['Variant_Frequencies_cf'] = multi_sample_df[vaf_columns[0]]
    else:
        compat_df['Variant_Frequencies_cf'] = 0.0
        
    if len(vaf_columns) >= 2:
        compat_df['Variant_Frequencies_st'] = multi_sample_df[vaf_columns[1]]
    else:
        compat_df['Variant_Frequencies_st'] = 0.0
    
    # Remove the VAF_sample_* columns to match old format
    for col in vaf_columns:
        compat_df = compat_df.drop(columns=[col])
    
    return compat_df


def validate_tree_compatibility(gene_list, tree_distribution):
    """
    Validate that the gene list is compatible with the tree distribution data.
    
    Args:
        gene_list: List of gene IDs from filtered DataFrame
        tree_distribution: Tree distribution data from aggregation
        
    Returns:
        bool: True if compatible, False otherwise
    """
    try:
        # Get the node dictionaries which contain the mutation assignments
        node_list = tree_distribution['node_dict']
        if not node_list:
            print("Warning: No node dictionary data found in tree distribution")
            return False
        
        # Extract all mutations from the first tree's node dictionary
        # (all trees should have the same mutations, just different assignments)
        first_tree_node_dict = node_list[0]
        tree_mutations = set()
        
        for node, mutations in first_tree_node_dict.items():
            tree_mutations.update(mutations)
        
        expected_gene_count = len(tree_mutations)
        actual_gene_count = len(gene_list)
        
        # Also check that the gene IDs match
        gene_set = set(gene_list)
        
        if expected_gene_count != actual_gene_count:
            print(f"ERROR: Gene count mismatch!")
            print(f"  Tree data expects: {expected_gene_count} genes")
            print(f"  Filtered data has: {actual_gene_count} genes")
            print(f"  This indicates the filtering produced a different mutation set than expected.")
            print(f"  The tree distribution was computed on a different set of mutations.")
            return False
        
        # Check if the gene sets match
        if tree_mutations != gene_set:
            missing_in_tree = gene_set - tree_mutations
            missing_in_genes = tree_mutations - gene_set
            print(f"ERROR: Gene set mismatch!")
            if missing_in_tree:
                print(f"  Genes in filtered data but not in tree: {sorted(missing_in_tree)}")
            if missing_in_genes:
                print(f"  Genes in tree but not in filtered data: {sorted(missing_in_genes)}")
            return False
            
        print(f"Validation successful: {expected_gene_count} genes match between tree and filtered data")
        return True
        
    except Exception as e:
        print(f"Error validating tree compatibility: {e}")
        return False


def main():
    args = parse_args()
    patient = args.patient
    read_depth = args.read_depth

    # Set up paths
    aggregation_dir = args.aggregation_dir
    ssm_file_path = args.ssm_file
    
    # Define output directory
    output_dir = args.output_dir if args.output_dir else os.path.join(args.aggregation_dir, 'marker_selection_output')
    os.makedirs(output_dir, exist_ok=True)

    # Define file paths
    tree_distribution_file = os.path.join(aggregation_dir, 'phylowgs_bootstrap_aggregation.pkl')

    # Verify files exist
    if not os.path.exists(tree_distribution_file):
        print(f"Error: Tree distribution file not found at {tree_distribution_file}")
        sys.exit(1)
    if not os.path.exists(ssm_file_path):
        print(f"Error: SSM file not found at {ssm_file_path}")
        sys.exit(1)

    # Load tree distribution from aggregation directory
    print("Loading tree distribution data...")
    with open(tree_distribution_file, 'rb') as f:
        tree_distribution = pickle.load(f)

    # Convert SSM file to multi-sample DataFrame format
    print(f"Converting SSM file with filtering strategy: {args.filter_strategy}")
    multi_sample_df = convert_ssm_to_dataframe_multi(
        ssm_file_path, 
        args.filter_strategy, 
        args.filter_threshold,
        args.filter_samples
    )
    
    print(f"Multi-sample DataFrame shape: {multi_sample_df.shape}")
    print("Multi-sample DataFrame columns:", list(multi_sample_df.columns))

    # Create backward-compatible DataFrame for the marker selection algorithms
    print("Creating backward-compatible DataFrame...")
    inter = create_backward_compatible_dataframe(multi_sample_df)
    
    print(f"Backward-compatible DataFrame shape: {inter.shape}")
    print("Backward-compatible DataFrame columns:", list(inter.columns))

    # Apply additional filtering if needed (this should be minimal since filtering was done in conversion)
    print("Applying final VAF filtering check...")
    original_count = len(inter)
    
    # Apply the same filtering logic as the old code (should be redundant but ensures compatibility)
    inter = inter[inter["Variant_Frequencies_cf"] < args.filter_threshold]
    inter = inter[inter["Variant_Frequencies_st"] < args.filter_threshold]
    
    filtered_count = len(inter)
    print(f"Final filtering check: {original_count} → {filtered_count} mutations")
    
    if filtered_count == 0:
        print("Error: All mutations were filtered out. Check VAF thresholds and filtering strategy.")
        sys.exit(1)
    
    calls = inter

    # Create gene indexing exactly as in old code
    gene2idx = {'s' + str(i): i for i in range(len(inter))}
    gene_list = list(gene2idx.keys())
    
    print(f"Created gene indexing: {len(gene_list)} genes (s0 to s{len(gene_list)-1})")

    # Validate compatibility with tree distribution
    if not validate_tree_compatibility(gene_list, tree_distribution):
        print("\nERROR: The filtered mutation set is incompatible with the tree distribution data.")
        print("This likely means the tree distribution was computed on a different set of mutations.")
        print("Possible solutions:")
        print("1. Use a different filtering strategy")
        print("2. Regenerate the tree distribution with the current mutation set")
        print("3. Check that the SSM file matches the one used for tree computation")
        sys.exit(1)

    # Create gene names exactly as in old code
    gene_name_list = []
    gene_count = {}

    for i in range(inter.shape[0]):
        gene = calls.iloc[i]["Hugo_Symbol"]
        ref = calls.iloc[i]["Reference_Allele"]
        alt = calls.iloc[i]["Allele"]
        
        if pd.isna(gene) or not isinstance(gene, str):
            # If no gene name, create a label with chromosome, position, and mutation
            chrom = str(calls.iloc[i]["Chromosome"])
            pos = str(calls.iloc[i]["Start_Position"])
            gene = f"Chr{chrom}:{pos}({ref}>{alt})"
        else:
            # If gene name exists, add mutation info and handle duplicates
            mutation = f"({ref}>{alt})"
            gene_with_mut = f"{gene}{mutation}"
            if gene_with_mut in gene_name_list:
                gene_count[gene_with_mut] = gene_count.get(gene_with_mut, 1) + 1
                gene = f"{gene_with_mut}_{gene_count[gene_with_mut]}"
            else:
                gene = gene_with_mut
        gene_name_list.append(gene)

    print(f"Created gene names: {len(gene_name_list)} entries")
    print("Sample gene names:", gene_name_list[:5])

    # Extract tree distribution components
    tree_list, node_list, clonal_freq_list, tree_freq_list = (
        tree_distribution['tree_structure'], 
        tree_distribution['node_dict'],
        tree_distribution['vaf_frac'],
        tree_distribution['freq']
    )

    # Scrub node_list (same as old code)
    node_list_scrub = []
    for node_dict in node_list:
        temp = {}
        for key, values in node_dict.items():
            temp.setdefault(int(key), values)
        node_list_scrub.append(temp)

    clonal_freq_list_scrub = []
    for clonal_freq_dict in clonal_freq_list:
        temp = {}
        for key, values in clonal_freq_dict.items():
            temp.setdefault(int(key), values[0])
        clonal_freq_list_scrub.append(temp)

    print(f"Tree distribution loaded: {len(tree_list)} trees, {len(node_list)} node sets")

    # Save marker selection results to a text file
    results_file = os.path.join(output_dir, f'{patient}_marker_selection_results.txt')
    with open(results_file, 'w') as f:
        f.write(f"Marker Selection Results for Patient {patient}\n")
        f.write("=" * 50 + "\n")
        f.write(f"Input: {ssm_file_path}\n")
        f.write(f"Filter strategy: {args.filter_strategy}\n")
        f.write(f"Filter threshold: {args.filter_threshold}\n")
        if args.filter_samples:
            f.write(f"Filter samples: {args.filter_samples}\n")
        f.write(f"Mutations after filtering: {len(gene_list)}\n")
        f.write(f"Read depth: {read_depth}\n\n")

    # Method 1: Tracing fractions
    print("Running Method 1: Tracing fractions...")
    print(f"Will iterate through {len(gene_name_list)} marker counts (1 to {len(gene_name_list)})")
    selected_markers1_genename_ordered = []
    obj1_ordered = []

    for n_markers in range(1, len(gene_name_list) + 1):
        selected_markers1, obj = select_markers_fractions_weighted_overall(
            gene_list, n_markers, tree_list, node_list_scrub, 
            clonal_freq_list_scrub, gene2idx, tree_freq_list)
        
        # Handle case where optimization failed and returned empty results
        if not selected_markers1 or any(pd.isna([obj])):
            print(f"Warning: Optimization failed for n_markers={n_markers}. Skipping this iteration.")
            print(f"Selected markers: {selected_markers1}, Objective: {obj}")
            break
            
        selected_markers1_genename = [gene_name_list[int(i[1:])] for i in selected_markers1]
        obj1_ordered.append(obj)
        
        if len(selected_markers1_genename) == 1:
            selected_markers1_genename_ordered.append(selected_markers1_genename[0])
        else:
            diff_set = set(selected_markers1_genename).difference(set(selected_markers1_genename_ordered))
            if diff_set:  # Check if diff_set is not empty
                selected_markers1_genename_ordered.append(list(diff_set)[0])
            else:
                print(f"Warning: No new markers found for n_markers={n_markers}. This may indicate optimization issues.")
                # Use the first marker from selected_markers1_genename as fallback
                if selected_markers1_genename:
                    selected_markers1_genename_ordered.append(selected_markers1_genename[0])
                else:
                    print(f"Error: No markers selected for n_markers={n_markers}. Breaking loop.")
                    break
    
    # Save Method 1 results
    print(f"Method 1 completed with {len(selected_markers1_genename_ordered)} successful iterations out of {len(gene_name_list)} attempted")
    
    with open(results_file, 'a') as f:
        f.write("Method 1 (Tracing Fractions) Results:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Completed {len(selected_markers1_genename_ordered)} iterations out of {len(gene_name_list)} attempted\n")
        for i, (marker, obj) in enumerate(zip(selected_markers1_genename_ordered, obj1_ordered), 1):
            # Get the index of this marker in gene_name_list
            marker_idx = gene_name_list.index(marker)
            # Get position info
            chrom = str(calls.iloc[marker_idx]["Chromosome"])
            pos = str(calls.iloc[marker_idx]["Start_Position"])
            f.write(f"{i}. {marker} [Chr{chrom}:{pos}]: {obj}\n")
        f.write("\n")

    # Plot Method 1 results (only if we have results)
    if selected_markers1_genename_ordered and obj1_ordered:
        position1 = list(range(len(obj1_ordered)))
        plt.figure(figsize=(8, 5))
        plt.plot(position1, obj1_ordered, 'o-', label='tracing-fractions')
        plt.xticks(position1, selected_markers1_genename_ordered, rotation=30)
        plt.legend()
        plt.title(f'Patient {patient} - Tracing Fractions ({args.filter_strategy})')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{patient}_tracing_subclones.png'), format='png', dpi=300, bbox_inches='tight')
        plt.close()
        print("Method 1 plot saved successfully")
    else:
        print("Warning: No Method 1 results to plot")

    # Method 2: Tree-based selection with different parameters
    for lam1, lam2 in [(1, 0), (0, 1)]:
        print(f"Running Method 2: Tree-based selection (lam1={lam1}, lam2={lam2})...")
        selected_markers2_genename_ordered = []
        obj2_ordered = []
        
        for n_markers in range(1, len(gene_name_list) + 1):
            selected_markers2, obj_frac, obj_struct = select_markers_tree_gp(
                gene_list, n_markers, tree_list, node_list_scrub, clonal_freq_list_scrub, 
                gene2idx, tree_freq_list, read_depth=read_depth, lam1=lam1, lam2=lam2
            )
            
            # Handle case where optimization failed and returned empty results
            if not selected_markers2 or any(pd.isna([obj_frac, obj_struct])):
                print(f"Warning: Tree optimization failed for n_markers={n_markers} (lam1={lam1}, lam2={lam2}). Skipping this iteration.")
                print(f"Selected markers: {selected_markers2}, Objectives: frac={obj_frac}, struct={obj_struct}")
                break
                
            selected_markers2_genename = [gene_name_list[int(i[1:])] for i in selected_markers2]
            obj2_ordered.append((obj_frac, obj_struct))
            
            if len(selected_markers2_genename) == 1:
                selected_markers2_genename_ordered.append(selected_markers2_genename[0])
            else:
                diff_set = set(selected_markers2_genename).difference(set(selected_markers2_genename_ordered))
                if diff_set:  # Check if diff_set is not empty
                    selected_markers2_genename_ordered.append(list(diff_set)[0])
                else:
                    print(f"Warning: No new markers found for n_markers={n_markers} (lam1={lam1}, lam2={lam2}). This may indicate optimization issues.")
                    # Use the first marker from selected_markers2_genename as fallback
                    if selected_markers2_genename:
                        selected_markers2_genename_ordered.append(selected_markers2_genename[0])
                    else:
                        print(f"Error: No markers selected for n_markers={n_markers} (lam1={lam1}, lam2={lam2}). Breaking loop.")
                        break

        # Save Method 2 results
        with open(results_file, 'a') as f:
            f.write(f"\nMethod 2 Results (lam1={lam1}, lam2={lam2}):\n")
            f.write("-" * 40 + "\n")
            for i, (marker, (obj_frac, obj_struct)) in enumerate(zip(selected_markers2_genename_ordered, obj2_ordered), 1):
                # Get the index of this marker in gene_name_list
                marker_idx = gene_name_list.index(marker)
                # Get position info
                chrom = str(calls.iloc[marker_idx]["Chromosome"])
                pos = str(calls.iloc[marker_idx]["Start_Position"])
                f.write(f"{i}. {marker} [Chr{chrom}:{pos}]: fraction={obj_frac}, structure={obj_struct}\n")
            f.write("\n")

        obj2_frac_ordered = [obj2_ordered[i][0] for i in range(len(obj2_ordered))]
        obj2_struct_ordered = [obj2_ordered[i][1] for i in range(len(obj2_ordered))]
        position2 = list(range(len(obj2_ordered)))

        # Plot fractions
        plt.figure(figsize=(8, 5))
        plt.plot(position2, obj2_frac_ordered, 'o-', color='tab:orange', label='trees-fractions')
        plt.xticks(position2, selected_markers2_genename_ordered, rotation=30)
        plt.legend()
        plt.title(f'Patient {patient} - Tree Fractions (λ1={lam1}, λ2={lam2}, {args.filter_strategy})')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{patient}_trees_fractions_{lam1}_{lam2}_{read_depth}.png'), format='png', dpi=300, bbox_inches='tight')
        plt.close()

        # Plot structures
        plt.figure(figsize=(8, 5))
        plt.plot(position2, obj2_struct_ordered, 'o-', color='tab:green', label='trees-structure')
        plt.xticks(position2, selected_markers2_genename_ordered, rotation=30)
        plt.legend()
        plt.title(f'Patient {patient} - Tree Structures (λ1={lam1}, λ2={lam2}, {args.filter_strategy})')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{patient}_trees_structures_{lam1}_{lam2}_{read_depth}.png'), format='png', dpi=300, bbox_inches='tight')
        plt.close()

    print(f"\nMarker selection completed successfully!")
    print(f"Results saved to: {results_file}")
    print(f"Plots saved to: {output_dir}")


if __name__ == "__main__":
    main() 