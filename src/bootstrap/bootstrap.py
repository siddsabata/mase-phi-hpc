import pandas as pd
import numpy as np
import os
import sys
import argparse
from pathlib import Path

"""
The purpose of this script is to perform bootstrapping on mutation data from an SSM file.

This script takes an SSM (Simple Somatic Mutation) file, typically used as input
for tools like PhyloWGS. It performs bootstrapping on the read counts for each
mutation across multiple samples (if present) to quantify uncertainty.

The bootstrapping process for each mutation's samples:
1. Resamples read depths for each sample while maintaining the total coverage across samples for that mutation.
2. Resamples variant frequencies for each sample using its newly bootstrapped depth.
3. Repeats this process N times to create N bootstrapped SSM files.

Usage: 
python bootstrap.py -i <input_ssm_file> -o <output_directory> -n <number_of_bootstraps>
"""

def bootstrap_va_dt(AF_list, Depth_list, bootstrap_num):
    """
    Advanced bootstrapping of both depths and allele frequencies for a single mutation across its samples.
    
    Args:
        AF_list (list): List of variant allele frequencies for each sample of a mutation.
        Depth_list (list): List of read depths for each sample of a mutation.
        bootstrap_num (int): Number of bootstrap samples to generate.
    
    Returns:
        tuple: (bootstrapped_VAFs, bootstrapped_depths)
               bootstrapped_VAFs is a NumPy array (num_samples, bootstrap_num)
               bootstrapped_depths is a NumPy array (num_samples, bootstrap_num)
    """
    AF_array = np.array(AF_list)
    Depth_array = np.array(Depth_list)
    
    # Ensure no zero depths in pvals for multinomial if total_depth > 0
    # If a depth is 0, its pval should be 0. If all depths are 0, pvals sum to 0.
    # np.random.multinomial handles pvals summing to < 1 by distributing the remainder.
    # If sum(pvals) is 0 (all depths are 0), multinomial will fail.
    total_depth_sum = np.sum(Depth_array)

    # pvals for multinomial distribution of reads
    if total_depth_sum == 0:
        # If total depth is zero, all new depths will be zero.
        # VAFs will be NaN or 0 depending on binomial sampling with 0 depth.
        pvals = np.zeros_like(Depth_array, dtype=float)
    else:
        pvals = Depth_array / total_depth_sum

    # Counter for attempts to get non-zero depths
    count = 0
    while True:
        count += 1
        # new_Depth_list_transposed shape: (bootstrap_num, num_samples)
        if total_depth_sum == 0 :
             new_Depth_list_transposed = np.zeros((bootstrap_num, len(Depth_array)), dtype=int)
        else:
            new_Depth_list_transposed = np.random.multinomial(n=int(total_depth_sum), 
                                                              pvals=pvals, 
                                                              size=bootstrap_num)
        
        # Break if no zero depths or if all original depths were zero
        if not np.any(new_Depth_list_transposed == 0) or total_depth_sum == 0:
            break
            
        # If attempts exceed threshold, replace zeros with 1s (if original depth > 0)
        if count >= 10:
            # Only replace where original depth was non-zero
            mask_original_nonzero = Depth_array > 0
            mask_new_zero = new_Depth_list_transposed == 0
            # Apply replacement only where original was >0 and new is 0
            replacement_mask = mask_new_zero & mask_original_nonzero[np.newaxis, :]
            new_Depth_list_transposed[replacement_mask] = 1
            break
    
    # AF_list_update shape: (num_samples, bootstrap_num)
    AF_list_update = np.zeros((len(AF_list), bootstrap_num))
    
    # Transpose new_Depth_list for easier iteration: (num_samples, bootstrap_num)
    new_Depth_list = new_Depth_list_transposed.T 

    for i in range(len(AF_list)): # Iterate over samples
        for j in range(bootstrap_num): # Iterate over bootstrap iterations
            current_depth = new_Depth_list[i, j]
            original_vaf = AF_list[i]
            
            if current_depth == 0:
                AF_list_update[i, j] = 0 # Or np.nan, if preferred
            else:
                # Binomial sample for variant reads
                variant_reads_sample = np.random.binomial(n=current_depth, 
                                                          p=original_vaf, 
                                                          size=1)[0]
                AF_list_update[i, j] = variant_reads_sample / current_depth
    
    return AF_list_update, new_Depth_list

def write_bootstrapped_ssm_file(mutations_for_this_bootstrap_iter, bootstrap_iteration_num, output_dir):
    """
    Writes a single bootstrapped SSM file for a given bootstrap iteration.
    
    Args:
        mutations_for_this_bootstrap_iter (list): List of dicts, each representing a mutation
                                                  with its bootstrapped 'a' and 'd' values.
        bootstrap_iteration_num (int): The 1-based number of this bootstrap iteration.
        output_dir (str): The base directory to write bootstrap subdirectories.
    """
    bootstrap_sub_dir = Path(output_dir) / f'bootstrap{bootstrap_iteration_num}'
    os.makedirs(bootstrap_sub_dir, exist_ok=True)
    
    ssm_file_path = bootstrap_sub_dir / 'ssm.txt'
    
    df_bootstrapped_ssm = pd.DataFrame(mutations_for_this_bootstrap_iter)
    # Ensure correct SSM column order
    if not df_bootstrapped_ssm.empty:
        df_bootstrapped_ssm = df_bootstrapped_ssm[['id', 'gene', 'a', 'd', 'mu_r', 'mu_v']]
    else: # Handle case of no valid mutations processed
        df_bootstrapped_ssm = pd.DataFrame(columns=['id', 'gene', 'a', 'd', 'mu_r', 'mu_v'])

    df_bootstrapped_ssm.to_csv(ssm_file_path, sep='\t', index=False)
    
    # Create empty CNV file (required by PhyloWGS)
    cnv_file_path = bootstrap_sub_dir / 'cnv.txt'
    cnv_file_path.touch()

def apply_vaf_prefiltering(input_ssm_df, threshold=0.9):
    """
    Apply VAF pre-filtering using "any_high" strategy.
    Removes mutations where ANY sample has VAF >= threshold.
    
    Args:
        input_ssm_df (pd.DataFrame): Input SSM DataFrame
        threshold (float): VAF threshold for filtering (default: 0.9)
        
    Returns:
        pd.DataFrame: Filtered SSM DataFrame
    """
    print(f"Applying VAF pre-filtering (threshold >= {threshold})...")
    original_count = len(input_ssm_df)
    
    filtered_mutations = []
    
    for index, row in input_ssm_df.iterrows():
        mutation_id_val = row['id']
        
        try:
            # Parse reference and total counts
            if isinstance(row['a'], str):
                ref_counts_str = row['a'].split(',')
                depth_counts_str = row['d'].split(',')
            else:
                ref_counts_str = [str(row['a'])]
                depth_counts_str = [str(row['d'])]
                
            if len(ref_counts_str) != len(depth_counts_str):
                print(f"Warning: Sample count mismatch for {mutation_id_val}, skipping")
                continue
                
            ref_counts = [int(c) for c in ref_counts_str]
            depth_counts = [int(d) for d in depth_counts_str]
            
        except (ValueError, AttributeError, TypeError):
            print(f"Warning: Could not parse counts for {mutation_id_val}, skipping")
            continue
        
        # Calculate VAFs for all samples
        vafs = []
        valid_sample_found = False
        
        for ref_count, total_depth in zip(ref_counts, depth_counts):
            if total_depth < 0 or ref_count < 0 or ref_count > total_depth:
                continue  # Skip invalid samples
                
            valid_sample_found = True
            if total_depth == 0:
                vafs.append(0.0)
            else:
                vaf = (total_depth - ref_count) / total_depth
                vafs.append(vaf)
        
        if not valid_sample_found:
            continue  # Skip mutations with no valid samples
        
        # Apply "any_high" filtering: keep if ALL VAFs < threshold
        if all(vaf < threshold for vaf in vafs):
            filtered_mutations.append(row)
        else:
            # Mutation filtered out due to high VAF
            pass
    
    filtered_df = pd.DataFrame(filtered_mutations)
    filtered_count = len(filtered_df)
    
    print(f"VAF pre-filtering: {original_count} → {filtered_count} mutations")
    print(f"Removed {original_count - filtered_count} mutations with VAF >= {threshold}")
    
    return filtered_df


def process_and_bootstrap_ssm(input_ssm_df, num_bootstraps, output_dir):
    """
    Processes an input SSM DataFrame, performs bootstrapping, and writes output SSM files.
    """
    all_bootstrapped_iterations_data = [[] for _ in range(num_bootstraps)]

    for index, row in input_ssm_df.iterrows():
        mutation_id_val = row['id']
        gene_val = row['gene']
        mu_r_val = row['mu_r']
        mu_v_val = row['mu_v']

        try:
            # Handle both single-sample (integers) and multi-sample (comma-separated strings) formats
            if isinstance(row['a'], str):
                # Multi-sample format: comma-separated string
                ref_counts_str = row['a'].split(',')
                depth_counts_str = row['d'].split(',')
            else:
                # Single-sample format: integer values
                ref_counts_str = [str(row['a'])]
                depth_counts_str = [str(row['d'])]

            if len(ref_counts_str) != len(depth_counts_str):
                print(f"Warning: Mismatch in number of samples for 'a' and 'd' in mutation {mutation_id_val}. Skipping.")
                continue

            ref_counts = [int(c) for c in ref_counts_str]
            depth_counts = [int(d_val) for d_val in depth_counts_str] # Renamed d to d_val to avoid conflict
        except ValueError:
            print(f"Warning: Could not parse 'a' or 'd' columns for mutation {mutation_id_val}. Skipping.")
            continue
        except (AttributeError, TypeError): # Handles case where 'a' or 'd' might be unexpected types
             print(f"Warning: 'a' or 'd' columns for mutation {mutation_id_val} have unexpected format. Skipping.")
             continue


        vaf_list_for_mutation = []
        depth_list_for_mutation = []
        valid_sample_found = False
        for r, d_sample in zip(ref_counts, depth_counts): # Renamed d to d_sample
            if d_sample < 0 or r < 0 or r > d_sample : # Basic sanity check
                 print(f"Warning: Invalid read counts (d={d_sample}, a={r}) for a sample in mutation {mutation_id_val}. Skipping this sample for this mutation.")
                 # We still need to append something to keep lists aligned for bootstrap_va_dt, or handle this inside.
                 # For now, let's skip mutations with any invalid sample. A more robust way would be to skip only the sample.
                 # However, bootstrap_va_dt expects aligned lists.
                 # The current structure implies all samples for a mutation are processed together.
                 # Let's append 0/0 for problematics samples but they might skew p_vals in multinomial
                 # A better approach for now is to skip the entire mutation if any sample within it is bad.
                 # To implement skipping only a sample, bootstrap_va_dt would need more complex input handling.
                 continue # This skips the current sample, will be handled by checking valid_sample_found later

            valid_sample_found = True # Mark that at least one sample is processable
            if d_sample == 0:
                vaf_list_for_mutation.append(0.0) # VAF is 0 if depth is 0
            else:
                vaf_list_for_mutation.append((d_sample - r) / d_sample)
            depth_list_for_mutation.append(d_sample)
        
        if not valid_sample_found : # or if len(vaf_list_for_mutation) == 0:
            print(f"Warning: No valid samples found for mutation {mutation_id_val} after parsing/validation. Skipping this mutation.")
            continue


        # boot_vaf_array: shape (num_samples, num_bootstraps)
        # boot_depth_array: shape (num_samples, num_bootstraps)
        # These arrays are for the current single mutation, across its samples and all bootstrap iterations
        boot_vaf_array, boot_depth_array = bootstrap_va_dt(vaf_list_for_mutation, depth_list_for_mutation, num_bootstraps)

        for k_bootstrap_iter in range(num_bootstraps):
            # Get VAFs and Depths for all samples of the current mutation for the k-th bootstrap
            current_iter_vafs = boot_vaf_array[:, k_bootstrap_iter]
            current_iter_depths = boot_depth_array[:, k_bootstrap_iter]

            # Calculate new variant and reference counts
            # np.round is important here as counts must be integers
            new_variant_counts = np.round(current_iter_vafs * current_iter_depths).astype(int)
            new_ref_counts = current_iter_depths - new_variant_counts
            
            # Ensure ref_counts are not negative (depths from bootstrap_va_dt should be >=0 after fix)
            # and also not greater than depth
            new_ref_counts = np.maximum(new_ref_counts, 0)
            new_ref_counts = np.minimum(new_ref_counts, current_iter_depths)


            bootstrapped_a_str = ",".join(map(str, new_ref_counts))
            bootstrapped_d_str = ",".join(map(str, current_iter_depths))

            all_bootstrapped_iterations_data[k_bootstrap_iter].append({
                'id': mutation_id_val,
                'gene': gene_val,
                'a': bootstrapped_a_str,
                'd': bootstrapped_d_str,
                'mu_r': mu_r_val,
                'mu_v': mu_v_val
            })

    # Now write out each bootstrapped SSM file
    print(f"Writing {num_bootstraps} bootstrapped SSM files to {output_dir}...")
    for i in range(num_bootstraps):
        write_bootstrapped_ssm_file(all_bootstrapped_iterations_data[i], i + 1, output_dir)
    print("Bootstrap processing complete.")


def main():
    parser = argparse.ArgumentParser(description='Bootstrap mutation data from an SSM file.')
    parser.add_argument('-i', '--input', required=True,
                       help='Input SSM file (tab-separated). E.g., ssm.txt')
    parser.add_argument('-o', '--output_dir', required=True,
                       help='Output directory for bootstrapped SSM files. Subdirectories (bootstrap1, bootstrap2, etc.) will be created here.')
    parser.add_argument('-n', '--num_bootstraps', type=int, default=100,
                       help='Number of bootstrap iterations (default: 100)')
    args = parser.parse_args()

    # Ensure output directory exists
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # Read input SSM data
    print(f"Reading input SSM file: {args.input}")
    try:
        input_ssm_df = pd.read_csv(args.input, sep='\t')
    except FileNotFoundError:
        print(f"Error: Input SSM file not found at {args.input}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading SSM file: {e}")
        sys.exit(1)
    
    if input_ssm_df.empty:
        print("Input SSM file is empty. Exiting.")
        sys.exit(1)

    # Check for required columns
    required_cols = ['id', 'gene', 'a', 'd', 'mu_r', 'mu_v']
    missing_cols = [col for col in required_cols if col not in input_ssm_df.columns]
    if missing_cols:
        print(f"Error: Input SSM file is missing required columns: {', '.join(missing_cols)}")
        sys.exit(1)

    # Apply VAF pre-filtering before bootstrap processing
    filtered_ssm_df = apply_vaf_prefiltering(input_ssm_df, threshold=0.9)
    
    if filtered_ssm_df.empty:
        print("Error: All mutations were filtered out by VAF pre-filtering. Check input data quality.")
        sys.exit(1)
    
    # Check minimum mutation count for meaningful phylogenetic analysis
    min_mutations = 5
    if len(filtered_ssm_df) < min_mutations:
        print(f"Error: After VAF filtering, only {len(filtered_ssm_df)} mutations remain.")
        print(f"Pipeline requires at least {min_mutations} mutations for meaningful phylogenetic analysis.")
        print("This typically indicates:")
        print("  - Input data has too many high-VAF mutations (likely artifacts)")
        print("  - Sample has very few somatic mutations")
        print("  - VAF filtering threshold may be too strict for this dataset")
        print("Pipeline terminated to prevent poor-quality results.")
        sys.exit(1)
    
    print(f"✅ Quality check passed: {len(filtered_ssm_df)} mutations available for phylogenetic analysis")
    
    # Save filtered SSM file for downstream stages
    # Extract patient base directory from output path (output_dir is typically {patient_base}/initial/bootstraps)
    output_path = Path(args.output_dir)
    if output_path.name == 'bootstraps' and output_path.parent.name == 'initial':
        # Standard pipeline structure: {patient_base}/initial/bootstraps
        patient_initial_dir = output_path.parent
        filtered_ssm_path = patient_initial_dir / 'ssm_filtered.txt'
    else:
        # Fallback: save in same directory as bootstrap output
        filtered_ssm_path = output_path / 'ssm_filtered.txt'
    
    print(f"Saving filtered SSM file to: {filtered_ssm_path}")
    filtered_ssm_df.to_csv(filtered_ssm_path, sep='\t', index=False)
    
    # Process bootstraps using filtered data
    process_and_bootstrap_ssm(filtered_ssm_df, args.num_bootstraps, args.output_dir)

if __name__ == "__main__":
    main() 