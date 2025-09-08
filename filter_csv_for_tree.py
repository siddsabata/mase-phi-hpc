#!/usr/bin/env python3
"""
Filter longitudinal CSV to match tree-compatible mutations.

This script identifies which mutations are actually used in the phylogenetic trees
and filters the longitudinal CSV to only include those mutations.
"""

import pandas as pd
import sys
from pathlib import Path

def filter_csv_for_tree_mutations(ssm_file: str, csv_file: str, output_file: str, num_tree_mutations: int = 5):
    """
    Filter CSV to only include the first N mutations from the SSM file
    that correspond to the tree-compatible mutations.
    
    Args:
        ssm_file: Path to SSM file (filtered)
        csv_file: Path to longitudinal CSV file  
        output_file: Path for filtered CSV output
        num_tree_mutations: Number of mutations actually used in trees (default 5)
    """
    print(f"Reading SSM file: {ssm_file}")
    ssm_df = pd.read_csv(ssm_file, sep='\t')
    print(f"SSM contains {len(ssm_df)} mutations")
    
    # Get the gene names for the first N mutations (s0, s1, s2, ...)
    tree_gene_names = []
    for i in range(min(num_tree_mutations, len(ssm_df))):
        gene_name = ssm_df.iloc[i]['gene']
        tree_gene_names.append(gene_name)
        print(f"s{i} -> {gene_name}")
    
    print(f"\nTree-compatible genes: {tree_gene_names}")
    
    # Read longitudinal CSV
    print(f"\nReading CSV file: {csv_file}")
    csv_df = pd.read_csv(csv_file)
    print(f"Original CSV shape: {csv_df.shape}")
    print(f"Original genes: {sorted(csv_df['gene'].unique())}")
    
    # Filter CSV to only include tree-compatible genes
    filtered_df = csv_df[csv_df['gene'].isin(tree_gene_names)].copy()
    print(f"Filtered CSV shape: {filtered_df.shape}")
    print(f"Filtered genes: {sorted(filtered_df['gene'].unique())}")
    
    # Save filtered CSV
    filtered_df.to_csv(output_file, index=False)
    print(f"\nSaved filtered CSV to: {output_file}")
    
    # Summary
    print(f"\nSummary:")
    print(f"  Original mutations in CSV: {len(csv_df['gene'].unique())}")
    print(f"  Tree-compatible mutations: {len(tree_gene_names)}")
    print(f"  Filtered mutations in CSV: {len(filtered_df['gene'].unique())}")
    print(f"  Missing from tree: {set(csv_df['gene'].unique()) - set(tree_gene_names)}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_csv_for_tree.py <ssm_file> <csv_file> <output_file>")
        print("Example: python filter_csv_for_tree.py /path/to/ssm_filtered.txt /path/to/liquid.csv /path/to/liquid_filtered.csv")
        sys.exit(1)
    
    ssm_file = sys.argv[1]
    csv_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Check files exist
    if not Path(ssm_file).exists():
        print(f"Error: SSM file not found: {ssm_file}")
        sys.exit(1)
    
    if not Path(csv_file).exists():
        print(f"Error: CSV file not found: {csv_file}")
        sys.exit(1)
        
    filter_csv_for_tree_mutations(ssm_file, csv_file, output_file)