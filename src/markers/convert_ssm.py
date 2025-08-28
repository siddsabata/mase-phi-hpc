#!/usr/bin/env python3
"""
Multi-sample SSM to DataFrame converter with flexible VAF filtering.

This version supports variable numbers of samples (2 to n) and provides
multiple filtering strategies for VAF-based mutation filtering.
"""

import pandas as pd
import argparse
import sys
from pathlib import Path


def parse_gene_string(gene_string):
    """Parse gene string in format: SYMBOL_CHR_POS_REF>ALT"""
    if pd.isna(gene_string) or not isinstance(gene_string, str):
        return {
            'symbol': 'Unknown',
            'chromosome': 'N/A',
            'position': 'N/A', 
            'ref_allele': 'N',
            'alt_allele': 'N'
        }
    
    parts = gene_string.split('_')
    
    if len(parts) >= 4:
        symbol = parts[0]
        chromosome = parts[1]
        position = parts[2]
        
        mutation_part = parts[3]
        if '>' in mutation_part:
            ref_allele, alt_allele = mutation_part.split('>', 1)
        else:
            ref_allele = 'N'
            alt_allele = 'N'
            
        return {
            'symbol': symbol,
            'chromosome': chromosome,
            'position': position,
            'ref_allele': ref_allele,
            'alt_allele': alt_allele
        }
    
    elif len(parts) == 1:
        return {
            'symbol': gene_string,
            'chromosome': 'N/A',
            'position': 'N/A',
            'ref_allele': 'N',
            'alt_allele': 'N'
        }
    
    else:
        return {
            'symbol': gene_string,
            'chromosome': 'N/A', 
            'position': 'N/A',
            'ref_allele': 'N',
            'alt_allele': 'N'
        }


def calculate_vaf(ref_count, total_depth):
    """Calculate Variant Allele Frequency (VAF)."""
    if total_depth == 0:
        return 0.0
    variant_count = total_depth - ref_count
    return variant_count / total_depth


def parse_ssm_counts(count_string):
    """Parse comma-separated count string from SSM file."""
    try:
        count_string = str(count_string).strip()
        if not count_string:
            return []
        counts = [int(x.strip()) for x in count_string.split(',') if x.strip()]
        return counts
    except (ValueError, AttributeError) as e:
        print(f"Warning: Could not parse count string '{count_string}': {e}")
        return []


def apply_vaf_filtering(df, strategy="any_high", threshold=0.9, specific_samples=None):
    """
    Apply VAF filtering with different strategies.
    
    Args:
        df: DataFrame with VAF columns (VAF_sample_0, VAF_sample_1, etc.)
        strategy: Filtering strategy
        threshold: VAF threshold for filtering
        specific_samples: List of sample indices for "specific_samples" strategy
    
    Returns:
        Filtered DataFrame
    """
    # Get VAF columns
    vaf_columns = [col for col in df.columns if col.startswith('VAF_sample_')]
    
    if not vaf_columns:
        print("Warning: No VAF columns found for filtering")
        return df
    
    original_count = len(df)
    
    if strategy == "any_high":
        # Filter if ANY sample VAF >= threshold
        mask = (df[vaf_columns] >= threshold).any(axis=1)
        df_filtered = df[~mask]
        
    elif strategy == "all_high":
        # Filter if ALL sample VAFs >= threshold
        mask = (df[vaf_columns] >= threshold).all(axis=1)
        df_filtered = df[~mask]
        
    elif strategy == "majority_high":
        # Filter if >50% of sample VAFs >= threshold
        high_count = (df[vaf_columns] >= threshold).sum(axis=1)
        mask = high_count > (len(vaf_columns) / 2)
        df_filtered = df[~mask]
        
    elif strategy == "specific_samples":
        # Filter based on specific sample indices
        if specific_samples is None:
            specific_samples = [0, 1]  # Default to first two samples
        
        specific_cols = [f'VAF_sample_{i}' for i in specific_samples if f'VAF_sample_{i}' in df.columns]
        if specific_cols:
            mask = (df[specific_cols] >= threshold).any(axis=1)
            df_filtered = df[~mask]
        else:
            print(f"Warning: Specified samples {specific_samples} not found in data")
            df_filtered = df
            
    else:
        print(f"Warning: Unknown filtering strategy '{strategy}'. No filtering applied.")
        df_filtered = df
    
    filtered_count = len(df_filtered)
    print(f"VAF filtering ({strategy}, threshold={threshold}): {original_count} → {filtered_count} mutations")
    
    return df_filtered


def convert_ssm_to_dataframe_multi(ssm_file_path, filter_strategy="any_high", filter_threshold=0.9, specific_samples=None):
    """
    Convert SSM file to DataFrame format with multi-sample support.
    
    Args:
        ssm_file_path: Path to SSM file
        filter_strategy: VAF filtering strategy
        filter_threshold: VAF threshold for filtering
        specific_samples: Sample indices for specific filtering
        
    Returns:
        pd.DataFrame: Converted DataFrame with all sample VAFs and metadata
    """
    # Read SSM file
    try:
        ssm_df = pd.read_csv(ssm_file_path, sep='\t')
    except Exception as e:
        print(f"Error reading SSM file {ssm_file_path}: {e}")
        sys.exit(1)
    
    # Validate required columns
    required_columns = ['id', 'gene', 'a', 'd', 'mu_r', 'mu_v']
    missing_columns = [col for col in required_columns if col not in ssm_df.columns]
    if missing_columns:
        print(f"Error: SSM file missing required columns: {missing_columns}")
        sys.exit(1)
    
    print(f"Processing {len(ssm_df)} mutations from SSM file...")
    
    # Determine number of samples from first mutation
    first_ref_counts = parse_ssm_counts(ssm_df.iloc[0]['a'])
    first_total_depths = parse_ssm_counts(ssm_df.iloc[0]['d'])
    num_samples = len(first_ref_counts)
    
    print(f"Detected {num_samples} samples per mutation")
    
    # Initialize output DataFrame
    output_data = []
    
    for idx, row in ssm_df.iterrows():
        # Parse gene information
        gene_info = parse_gene_string(row['gene'])
        
        # Parse read counts
        ref_counts = parse_ssm_counts(row['a'])
        total_depths = parse_ssm_counts(row['d'])
        
        # Validate count data
        if len(ref_counts) != len(total_depths):
            print(f"Warning: Mismatch in ref/total counts for mutation {row['id']}. Skipping.")
            continue
            
        if len(ref_counts) == 0:
            print(f"Warning: No valid count data for mutation {row['id']}. Skipping.")
            continue
        
        if len(ref_counts) != num_samples:
            print(f"Warning: Sample count mismatch for mutation {row['id']} (expected {num_samples}, got {len(ref_counts)}). Skipping.")
            continue
        
        # Calculate VAFs for all samples
        vafs = []
        for ref_count, total_depth in zip(ref_counts, total_depths):
            if total_depth < 0 or ref_count < 0 or ref_count > total_depth:
                print(f"Warning: Invalid counts (ref={ref_count}, total={total_depth}) for {row['id']}. Using VAF=0.")
                vafs.append(0.0)
            else:
                vaf = calculate_vaf(ref_count, total_depth)
                vafs.append(vaf)
        
        # Create output row with all metadata and VAFs
        output_row = {
            'Hugo_Symbol': gene_info['symbol'],
            'Reference_Allele': gene_info['ref_allele'],
            'Allele': gene_info['alt_allele'],
            'Chromosome': gene_info['chromosome'],
            'Start_Position': gene_info['position'],
        }
        
        # Add VAF columns for all samples
        for i, vaf in enumerate(vafs):
            output_row[f'VAF_sample_{i}'] = vaf
        
        output_data.append(output_row)
    
    # Create DataFrame
    output_df = pd.DataFrame(output_data)
    
    print(f"Successfully converted {len(output_df)} mutations to DataFrame format.")
    
    # Apply VAF filtering
    filtered_df = apply_vaf_filtering(output_df, filter_strategy, filter_threshold, specific_samples)
    
    return filtered_df


def main():
    """Main function to handle command line interface."""
    parser = argparse.ArgumentParser(
        description='Convert SSM files to DataFrame format with multi-sample support',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Filtering Strategies:
    any_high     - Filter if ANY sample VAF >= threshold (most restrictive)
    all_high     - Filter if ALL sample VAFs >= threshold (least restrictive)  
    majority_high - Filter if >50% of sample VAFs >= threshold
    specific_samples - Filter based on specific sample indices

Examples:
    # Convert with default filtering (any_high, threshold=0.9)
    python convert_ssm_to_dataframe_v2.py data/ssm.txt output.csv
    
    # Use specific samples for filtering
    python convert_ssm_to_dataframe_v2.py data/ssm.txt output.csv --filter-strategy specific_samples --specific-samples 0 2
    
    # Use majority-based filtering
    python convert_ssm_to_dataframe_v2.py data/ssm.txt output.csv --filter-strategy majority_high --threshold 0.8
        """
    )
    
    parser.add_argument('input_ssm', type=str,
                       help='Input SSM file path')
    
    parser.add_argument('output_csv', type=str, nargs='?',
                       help='Output CSV file path (required unless --info-only)')
    
    parser.add_argument('--filter-strategy', type=str, 
                       choices=['any_high', 'all_high', 'majority_high', 'specific_samples'],
                       default='any_high',
                       help='VAF filtering strategy (default: any_high)')
    
    parser.add_argument('--threshold', type=float, default=0.9,
                       help='VAF threshold for filtering (default: 0.9)')
    
    parser.add_argument('--specific-samples', type=int, nargs='+',
                       help='Sample indices for specific_samples strategy')
    
    parser.add_argument('--info-only', action='store_true',
                       help='Only show sample information, do not convert')
    
    args = parser.parse_args()
    
    # Validate input file
    if not Path(args.input_ssm).exists():
        print(f"Error: Input SSM file not found: {args.input_ssm}")
        sys.exit(1)
    
    # Show sample information if requested
    if args.info_only:
        try:
            ssm_df = pd.read_csv(args.input_ssm, sep='\t')
            print(f"SSM file: {args.input_ssm}")
            print(f"Number of mutations: {len(ssm_df)}")
            
            # Analyze sample counts
            if 'a' in ssm_df.columns and 'd' in ssm_df.columns:
                sample_counts = []
                for _, row in ssm_df.head().iterrows():
                    ref_counts = parse_ssm_counts(row['a'])
                    sample_counts.append(len(ref_counts))
                
                if sample_counts:
                    print(f"Number of samples per mutation: {sample_counts[0]} (based on first mutation)")
                    print(f"Available sample indices: 0 to {sample_counts[0]-1}")
                else:
                    print("Could not determine number of samples")
            
            print("\nFirst few mutations:")
            print(ssm_df[['id', 'gene']].head())
            
        except Exception as e:
            print(f"Error reading SSM file: {e}")
            sys.exit(1)
        
        return
    
    # Validate output file argument
    if not args.output_csv:
        print("Error: Output CSV file path is required (unless using --info-only)")
        sys.exit(1)
    
    # Convert SSM to DataFrame
    try:
        output_df = convert_ssm_to_dataframe_multi(
            args.input_ssm, 
            args.filter_strategy, 
            args.threshold,
            args.specific_samples
        )
        
        # Save to CSV
        output_df.to_csv(args.output_csv, index=False)
        print(f"Converted DataFrame saved to: {args.output_csv}")
        
        # Show summary statistics
        print("\nConversion Summary:")
        print(f"Total mutations after filtering: {len(output_df)}")
        
        # Show VAF ranges for all samples
        vaf_columns = [col for col in output_df.columns if col.startswith('VAF_sample_')]
        print(f"VAF ranges across {len(vaf_columns)} samples:")
        for col in vaf_columns:
            sample_idx = col.split('_')[-1]
            print(f"  Sample {sample_idx}: {output_df[col].min():.3f} - {output_df[col].max():.3f}")
        
    except Exception as e:
        print(f"Error during conversion: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 