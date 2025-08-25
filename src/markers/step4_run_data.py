from step4_optimize import *
from step4_optimize_fraction import *
import pandas as pd
from zipfile import ZipFile
import json
import gzip
import pickle
from pathlib import Path
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

def parse_gene_info(gene_string):
    """
    Parses a gene string (expected format: SYMBOL_CHR_POS_REF>ALT or just SYMBOL)
    into its components.
    Returns a dictionary with 'Symbol', 'Chromosome', and 'Start_Position'.
    Returns 'N/A' for fields that cannot be parsed.
    """
    parts = gene_string.split('_')
    if len(parts) >= 3:  # SYMBOL_CHR_POS...
        # Check if CHR and POS look like valid numbers, otherwise treat as part of symbol
        if parts[1].isdigit() and parts[2].isdigit():
            symbol = parts[0]
            chromosome = parts[1]
            position = parts[2]
            # Potentially more parts for REF>ALT etc.
            return {'Symbol': symbol, 'Chromosome': chromosome, 'Start_Position': position}
        else: # CHR or POS is not a simple number, assume it's part of a complex symbol
            return {'Symbol': gene_string, 'Chromosome': 'N/A', 'Start_Position': 'N/A'}
    elif len(parts) == 1:  # Just SYMBOL
        symbol = parts[0]
        return {'Symbol': symbol, 'Chromosome': 'N/A', 'Start_Position': 'N/A'}
    else:  # Other unexpected formats or too few parts after split (e.g. "SYMBOL_")
        return {'Symbol': gene_string, 'Chromosome': 'N/A', 'Start_Position': 'N/A'}

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Run marker selection analysis.')
    
    parser.add_argument('patient', type=str,
                      help='Patient ID')
    
    parser.add_argument('--read-depth', type=int, default=1500,
                      help='Read depth for analysis (default: 1500)')
    
    parser.add_argument('--aggregation-dir', type=str, required=True,
                        help='Path to the directory containing phylowgs_bootstrap_aggregation.pkl')
    
    parser.add_argument('--ssm-file', type=str, required=True,
                        help='Path to the ssm.txt file for this patient')
    
    parser.add_argument('--output-dir', type=str,
                        help='Path to output directory for marker selection results (default: {aggregation-dir}/marker_selection_output)')
    
    return parser.parse_args()

def main():
    args = parse_args()
    patient = args.patient
    read_depth = args.read_depth

    # Set up paths
    aggregation_dir = args.aggregation_dir
    ssm_file_path = args.ssm_file
    
    # Define output directory - either use provided or default
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
    with open(tree_distribution_file, 'rb') as f:
        tree_distribution = pickle.load(f)

    gene_list = []
    gene2idx = {}

    # Read from ssm.txt file
    ssm_df = pd.read_csv(ssm_file_path, sep='\t')

    # Parse gene information to create Chromosome and Start_Position columns
    parsed_info_list = ssm_df['gene'].apply(parse_gene_info)
    # ssm_df['Symbol_parsed'] = [info['Symbol'] for info in parsed_info_list] # If parsed symbol is needed later
    ssm_df['Chromosome'] = [info['Chromosome'] for info in parsed_info_list]
    ssm_df['Start_Position'] = [info['Start_Position'] for info in parsed_info_list]

    # --- START VAF Calculation and Filtering ---
    def get_vaf_list_for_filtering(row):
        """Helper function to calculate VAFs for a given mutation row from 'a' and 'd' columns."""
        try:
            # Ensure 'a' and 'd' are treated as strings for splitting
            a_counts_str = str(row['a']).split(',')
            d_counts_str = str(row['d']).split(',')

            # Handle cases where columns might be empty or just whitespace after split
            a_counts = [int(x) for x in a_counts_str if x.strip()]
            d_counts = [int(x) for x in d_counts_str if x.strip()]

            vafs = []
            if len(a_counts) != len(d_counts):
                # This case should ideally not happen with well-formed ssm.txt
                # print(f"Warning: Mismatch in a/d counts for row: {row.to_dict()}") # Optional warning
                return [] # Return empty list, will lead to filtering out this mutation by default VAFs
            
            for ref_r, tot_d in zip(a_counts, d_counts):
                if tot_d > 0:
                    vaf = (tot_d - ref_r) / tot_d
                    vafs.append(vaf)
                else:
                    vafs.append(0.0) # Or handle as per desired logic, e.g., np.nan then fillna
            return vafs
        except ValueError:
            # print(f"Warning: ValueError during VAF calculation for row: {row.to_dict()}") # Optional warning
            return [] # Error in parsing counts, treat as if no VAFs calculable
        except Exception as e:
            # print(f"Warning: Unexpected error {e} during VAF calculation for row: {row.to_dict()}") # Optional warning
            return []

    ssm_df['vaf_list_for_filter'] = ssm_df.apply(get_vaf_list_for_filtering, axis=1)

    # Remove the specific VAF_filter_s1 and VAF_filter_s2 columns
    # ssm_df['VAF_filter_s1'] = ssm_df['vaf_list_for_filter'].apply(lambda x: x[0] if len(x) > 0 else 1.0)
    # ssm_df['VAF_filter_s2'] = ssm_df['vaf_list_for_filter'].apply(lambda x: x[1] if len(x) > 1 else 1.0)

    original_mutation_count = len(ssm_df)
    
    # Apply VAF-based pre-filtering: if ANY sample VAF is >= 0.9, exclude the mutation.
    # A mutation is kept if ALL its VAFs are < 0.9 (or if it has no VAFs, in which case it's kept by default of an empty list). 
    # The get_vaf_list_for_filtering handles errors by returning empty list, which will pass this filter.
    # If a mutation must have VAFs to be kept, the condition could be `ssm_df['vaf_list_for_filter'].apply(lambda x: len(x) > 0 and all(vaf < 0.9 for vaf in x))`
    # For now, we keep it simple: if any VAF is high, filter out. If no VAFs or all VAFs low, keep.
    
    # Define a condition to filter: True if mutation should be KEPT (all VAFs < 0.9)
    # False if mutation should be DISCARDED (at least one VAF >= 0.9 OR no VAFs could be calculated and list is empty)
    # To correctly discard if ANY vaf >= 0.9:
    # We want to keep rows where NOT (ANY vaf in list is >= 0.9)
    # This is equivalent to: ALL vafs in list are < 0.9
    # If vaf_list_for_filter is empty (e.g. due to parsing error), all() on empty list is True, so it would be kept.
    # This might be okay, or we might want to explicitly filter out rows with empty vaf_list_for_filter if they signify problematic data.
    # Let's ensure that if vaf_list_for_filter is empty, it's treated as if it fails the < 0.9 condition for safety, effectively filtering it.
    # A simple way: keep if list is not empty AND all vafs are < 0.9.
    
    def should_keep_mutation(vaf_list):
        if not vaf_list: # If list is empty (e.g., parsing error, or no samples where VAFs could be computed)
            return True   # Keep the mutation (conservative approach, don't filter if VAFs are unknown/unparseable)
        return all(vaf < 0.9 for vaf in vaf_list) # Filter if any VAF is >= 0.9

    ssm_df = ssm_df[ssm_df['vaf_list_for_filter'].apply(should_keep_mutation)]
    
    # Debugging: Print out the filtered ssm_df to check formatting and contents
    print("\n[DEBUG] Filtered ssm_df (first 10 rows):")
    print(ssm_df.head(10))
    print("\n[DEBUG] ssm_df columns:", list(ssm_df.columns))
    print(f"[DEBUG] Number of mutations after filtering: {len(ssm_df)}\n")
    
    # Old filter lines:
    # ssm_df = ssm_df[ssm_df["VAF_filter_s1"] < 0.9]
    # ssm_df = ssm_df[ssm_df["VAF_filter_s2"] < 0.9]
    filtered_mutation_count = len(ssm_df)
    
    print(f"VAF PRE-FILTERING: Started with {original_mutation_count} mutations, {filtered_mutation_count} remaining after VAF < 0.9 filter across all samples.")
    if original_mutation_count > 0 and filtered_mutation_count == 0:
        print("Warning: VAF pre-filtering removed all mutations. Check VAF values and filtering thresholds.")
    # --- END VAF Calculation and Filtering ---

    calls = ssm_df # Use filtered ssm_df as the base for mutation data

    # Derive gene_list and gene2idx from the potentially filtered ssm_df ('id' column)
    gene_list = ssm_df['id'].tolist()
    gene2idx = {gene_id: i for i, gene_id in enumerate(gene_list)}

    gene_name_list = []
    gene_count = {}

    # Construct gene_name_list from ssm_df ('gene' column), ensuring uniqueness
    # The 'gene' column in ssm.txt is expected to be in SYMBOL_CHR_POS_REF>ALT format
    # or just SYMBOL. The existing logic handles creating unique names if needed.
    for i in range(calls.shape[0]):
        # The 'gene' column from ssm.txt is used directly here.
        # It might be just a SYMBOL or a more complex string like SYMBOL_CHR_POS_REF>ALT.
        # The original script had logic to parse Symbol, Ref, Alt from different columns
        # to form a base name. Here, we assume `calls.iloc[i]["gene"]` is the intended base name.
        gene_base_name = calls.iloc[i]["gene"] 

        if pd.isna(gene_base_name) or not isinstance(gene_base_name, str):
            # Fallback if 'gene' column is problematic or missing, though ssm.txt should have it.
            # This attempts to create a fallback similar to the original script's handling of missing Hugo_Symbol.
            # Requires parsing from the 'gene' column if it holds complex ID, or using 'id' if simple.
            # For ssm.txt, the 'id' column is s0, s1 etc. and 'gene' is the descriptive one.
            # If gene_base_name is NaN, we might need a robust way to get Chrom/Pos/Ref/Alt for the name.
            # However, the plan expects `ssm_df['gene']` to be the source for `gene_name_list`.
            # Let's assume `calls.iloc[i]["gene"]` is the primary source and is well-formed.
            # If it were truly NaN, we'd need a different strategy or error handling.
            # For now, let's stick to the plan of using the 'gene' column primarily.
            # The original script had a specific way to construct name from chr/pos/ref/alt if Hugo_Symbol was NaN.
            # We will adapt this if `gene_base_name` doesn't fit the expected pattern for unique naming. 
            # For now, if gene_base_name itself is an issue, this part might need refinement based on actual ssm.txt content.
            # Let's assume ssm_df['gene'] is the rich name like 'GPR65_14_88477948_G>T'
            # If it's just 'GPR65', the original logic for adding (ref>alt) is not directly applicable
            # without parsing ref/alt from somewhere else (e.g. other columns or the rich name itself if available).
            # Given the ssm.txt example, 'gene' column IS the rich name.
            
            # If gene name is NaN or not a string, construct a unique identifier from 'id'
            # This part of the original code is less likely to be hit if ssm.txt is well-formed,
            # but kept for robustness. The original used Chrom/Pos/Ref/Alt for this.
            # Since we don't have those parsed yet, we'll use the 'id' as a fallback name for uniqueness.
            gene_unique_name_candidate = str(calls.iloc[i]["id"]) # Fallback to ssm_df 'id'
        else:
            gene_unique_name_candidate = gene_base_name # This is typically SYMBOL_CHR_POS_REF>ALT

        # Original uniqueness logic: if a name (potentially with mutation info) is repeated,
        # append a counter. This should apply to the names from `ssm_df['gene']`.
        if gene_unique_name_candidate in gene_name_list:
            gene_count[gene_unique_name_candidate] = gene_count.get(gene_unique_name_candidate, 1) + 1
            final_gene_name = f"{gene_unique_name_candidate}_{gene_count[gene_unique_name_candidate]}"
        else:
            final_gene_name = gene_unique_name_candidate
        gene_name_list.append(final_gene_name)

    tree_list, node_list, clonal_freq_list, tree_freq_list = tree_distribution['tree_structure'], tree_distribution['node_dict'],tree_distribution['vaf_frac'],tree_distribution['freq']

    #scrub node_list
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

    # Run marker selection with different methods and parameters
    # Save marker selection results to a text file
    results_file = os.path.join(output_dir, f'{patient}_marker_selection_results.txt')
    with open(results_file, 'w') as f:
        f.write(f"Marker Selection Results for Patient {patient}\n")
        f.write("=" * 50 + "\n\n")

    # Method 1: Tracing fractions
    selected_markers1_genename_ordered = []
    obj1_ordered = []

    for n_markers in range(1, len(gene_name_list) + 1):
        selected_markers1, obj = select_markers_fractions_weighted_overall(gene_list, n_markers, tree_list, node_list_scrub, clonal_freq_list_scrub, gene2idx, tree_freq_list)
        selected_markers1_genename = [gene_name_list[int(i[1:])] for i in selected_markers1]
        obj1_ordered.append(obj)
        if len(selected_markers1_genename) == 1:
            selected_markers1_genename_ordered.append(selected_markers1_genename[0])
        else:
            diff_set = set(selected_markers1_genename).difference(set(selected_markers1_genename_ordered))
            selected_markers1_genename_ordered.append(list(diff_set)[0])
    
    # Save Method 1 results
    with open(results_file, 'a') as f:
        f.write("Method 1 (Tracing Fractions) Results:\n")
        f.write("-" * 40 + "\n")
        for i, (marker, obj) in enumerate(zip(selected_markers1_genename_ordered, obj1_ordered), 1):
            # Get the index of this marker in gene_name_list
            marker_idx = gene_name_list.index(marker)
            # Get position info
            chrom = str(calls.iloc[marker_idx]["Chromosome"])
            pos = str(calls.iloc[marker_idx]["Start_Position"])
            f.write(f"{i}. {marker} [Chr{chrom}:{pos}]: {obj}\n")
        f.write("\n")

    position1 = list(range(len(obj1_ordered)))
    plt.figure(figsize=(8, 5))
    plt.plot(position1, obj1_ordered, 'o-', label='tracing-fractions')
    plt.xticks(position1, selected_markers1_genename_ordered, rotation=30)
    plt.legend()
    plt.savefig(os.path.join(output_dir, f'{patient}_tracing_subclones.png'), format='png', dpi=300, bbox_inches='tight')
    plt.close()

    # Method 2: Tree-based selection with different parameters
    for lam1, lam2 in [(1, 0), (0, 1)]:
        selected_markers2_genename_ordered = []
        obj2_ordered = []
        
        for n_markers in range(1, len(gene_name_list) + 1):
            selected_markers2, obj_frac, obj_struct = select_markers_tree_gp(
                gene_list, n_markers, tree_list, node_list_scrub, clonal_freq_list_scrub, 
                gene2idx, tree_freq_list, read_depth=read_depth, lam1=lam1, lam2=lam2
            )
            selected_markers2_genename = [gene_name_list[int(i[1:])] for i in selected_markers2]
            obj2_ordered.append((obj_frac, obj_struct))
            if len(selected_markers2_genename) == 1:
                selected_markers2_genename_ordered.append(selected_markers2_genename[0])
            else:
                selected_markers2_genename_ordered.append(
                    list(set(selected_markers2_genename).difference(set(selected_markers2_genename_ordered)))[0])

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
        plt.savefig(os.path.join(output_dir, f'{patient}_trees_fractions_{lam1}_{lam2}_{read_depth}.png'), format='png', dpi=300, bbox_inches='tight')
        plt.close()

        # Plot structures
        plt.figure(figsize=(8, 5))
        plt.plot(position2, obj2_struct_ordered, 'o-', color='tab:green', label='trees-structure')
        plt.xticks(position2, selected_markers2_genename_ordered, rotation=30)
        plt.legend()
        plt.savefig(os.path.join(output_dir, f'{patient}_trees_structures_{lam1}_{lam2}_{read_depth}.png'), format='png', dpi=300, bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    main()
