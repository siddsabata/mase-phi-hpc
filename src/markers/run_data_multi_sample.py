#!/usr/bin/env python3
"""
Multi-sample compatible marker selection script.

This script works with the new multi-sample DataFrame format while maintaining
compatibility with the original marker selection algorithms and tree distribution data.
"""

from optimize import *
from optimize_fraction import select_markers_tree_gp
# Note: VAF filtering now handled in bootstrap stage, no need for convert_ssm
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
    
    # Note: VAF filtering is now handled in the bootstrap stage
    # The SSM file passed to this script should already be pre-filtered
    
    return parser.parse_args()


# Obsolete functions removed - VAF filtering now handled in bootstrap stage
# create_backward_compatible_dataframe() - No longer needed
# validate_tree_compatibility() - Tree compatibility guaranteed by design


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

    # Read pre-filtered SSM file directly (filtering done in bootstrap stage)
    print(f"Reading pre-filtered SSM file: {ssm_file_path}")
    try:
        ssm_df = pd.read_csv(ssm_file_path, sep='\t')
    except Exception as e:
        print(f"Error reading SSM file: {e}")
        sys.exit(1)
    
    print(f"Pre-filtered SSM DataFrame shape: {ssm_df.shape}")
    
    if ssm_df.empty:
        print("Error: SSM file is empty.")
        sys.exit(1)
    
    # Parse gene information from SSM file for backward compatibility
    def parse_gene_info(gene_string):
        """Parse gene string (expected format: SYMBOL_CHR_POS_REF>ALT or just SYMBOL)"""
        if pd.isna(gene_string) or not isinstance(gene_string, str):
            return {'Symbol': 'Unknown', 'Chromosome': 'N/A', 'Start_Position': 'N/A', 'Ref': 'N', 'Alt': 'N'}
        
        parts = gene_string.split('_')
        if len(parts) >= 4:
            symbol = parts[0]
            chromosome = parts[1] 
            position = parts[2]
            mutation_part = parts[3]
            if '>' in mutation_part:
                ref, alt = mutation_part.split('>', 1)
            else:
                ref, alt = 'N', 'N'
            return {'Symbol': symbol, 'Chromosome': chromosome, 'Start_Position': position, 'Ref': ref, 'Alt': alt}
        else:
            return {'Symbol': gene_string, 'Chromosome': 'N/A', 'Start_Position': 'N/A', 'Ref': 'N', 'Alt': 'N'}
    
    # Create backward-compatible DataFrame structure
    gene_info_list = ssm_df['gene'].apply(parse_gene_info)
    
    # Calculate VAFs from 'a' and 'd' columns for first two samples (backward compatibility)
    def calculate_sample_vafs(row):
        try:
            a_counts = [int(x) for x in str(row['a']).split(',') if x.strip()]
            d_counts = [int(x) for x in str(row['d']).split(',') if x.strip()]
            
            vafs = []
            for ref_count, total_depth in zip(a_counts, d_counts):
                if total_depth > 0:
                    vaf = (total_depth - ref_count) / total_depth
                    vafs.append(vaf)
                else:
                    vafs.append(0.0)
            return vafs
        except:
            return [0.0, 0.0]
    
    vaf_lists = ssm_df.apply(calculate_sample_vafs, axis=1)
    
    # Create DataFrame with marker selection expected format
    inter_data = []
    for i, (idx, row) in enumerate(ssm_df.iterrows()):
        gene_info = gene_info_list.iloc[i]
        vafs = vaf_lists.iloc[i]
        
        inter_data.append({
            'Hugo_Symbol': gene_info['Symbol'],
            'Reference_Allele': gene_info['Ref'],
            'Allele': gene_info['Alt'],
            'Chromosome': gene_info['Chromosome'],
            'Start_Position': gene_info['Start_Position'],
            'Variant_Frequencies_cf': vafs[0] if len(vafs) > 0 else 0.0,
            'Variant_Frequencies_st': vafs[1] if len(vafs) > 1 else 0.0
        })
    
    inter = pd.DataFrame(inter_data)
    print(f"Created backward-compatible DataFrame with {len(inter)} mutations")
    
    calls = inter

    # Create gene indexing exactly as in old code
    gene2idx = {'s' + str(i): i for i in range(len(inter))}
    gene_list = list(gene2idx.keys())
    
    print(f"Created gene indexing: {len(gene_list)} genes (s0 to s{len(gene_list)-1})")

    # Tree compatibility is guaranteed since bootstrap, PhyloWGS, and marker selection
    # all use the same pre-filtered mutation set from the bootstrap stage

    # Use original SSM gene identifiers directly (already informative)
    gene_name_list = ssm_df['gene'].tolist()

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
        f.write(f"Input (pre-filtered): {ssm_file_path}\n")
        f.write(f"Mutations in analysis: {len(gene_list)}\n")
        f.write(f"Read depth: {read_depth}\n\n")

    # Tree-based marker selection using two optimization strategies

    # Tree-based selection with two optimization strategies
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

        # Save Method 2 results with descriptive headers
        with open(results_file, 'a') as f:
            if lam1 == 1 and lam2 == 0:
                f.write(f"\nλ1=1, λ2=0 (Pure Fraction Optimization):\n")
            elif lam1 == 0 and lam2 == 1:
                f.write(f"\nλ1=0, λ2=1 (Pure Structure Optimization):\n")
            else:
                f.write(f"\nλ1={lam1}, λ2={lam2} (Mixed Optimization):\n")
            f.write("-" * 40 + "\n")
            for i, (marker, (obj_frac, obj_struct)) in enumerate(zip(selected_markers2_genename_ordered, obj2_ordered), 1):
                f.write(f"{i}. {marker}: fraction={obj_frac}, structure={obj_struct}\n")
            f.write("\n")

        obj2_frac_ordered = [obj2_ordered[i][0] for i in range(len(obj2_ordered))]
        obj2_struct_ordered = [obj2_ordered[i][1] for i in range(len(obj2_ordered))]
        position2 = list(range(len(obj2_ordered)))

        # Plot fractions
        plt.figure(figsize=(8, 5))
        plt.plot(position2, obj2_frac_ordered, 'o-', color='tab:orange', label='trees-fractions')
        plt.xticks(position2, selected_markers2_genename_ordered, rotation=30)
        plt.legend()
        plt.title(f'Patient {patient} - Tree Fractions (λ1={lam1}, λ2={lam2}, VAF Pre-filtered)')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{patient}_trees_fractions_{lam1}_{lam2}_{read_depth}.png'), format='png', dpi=300, bbox_inches='tight')
        plt.close()

        # Plot structures
        plt.figure(figsize=(8, 5))
        plt.plot(position2, obj2_struct_ordered, 'o-', color='tab:green', label='trees-structure')
        plt.xticks(position2, selected_markers2_genename_ordered, rotation=30)
        plt.legend()
        plt.title(f'Patient {patient} - Tree Structures (λ1={lam1}, λ2={lam2}, VAF Pre-filtered)')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{patient}_trees_structures_{lam1}_{lam2}_{read_depth}.png'), format='png', dpi=300, bbox_inches='tight')
        plt.close()

    print(f"\nMarker selection completed successfully!")
    print(f"Results saved to: {results_file}")
    print(f"Plots saved to: {output_dir}")


if __name__ == "__main__":
    main() 