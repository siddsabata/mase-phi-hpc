#!/usr/bin/env python3
"""
Data loading for longitudinal cancer evolution analysis.

This module handles loading of tree distributions, tissue data from SSM files,
and longitudinal ddPCR data from CSV files.

Authors: TracerX Pipeline Development Team
"""

import logging
import pickle
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)


def load_tree_distributions(aggregation_dir: Path, method: str, logger: logging.Logger) -> Tuple[Dict, Dict]:
    """
    Load tree distribution data from aggregation stage outputs.
    
    Args:
        aggregation_dir: Path to aggregation results directory
        method: Phylogenetic method used (e.g., 'phylowgs')
        logger: Logger instance for reporting
        
    Returns:
        Tuple of (tree_distribution_summary, tree_distribution_full)
    """
    logger.info("Loading tree distribution data from aggregation stage...")
    
    # Load summary tree distribution
    summary_file = aggregation_dir / f'{method}_bootstrap_summary.pkl'
    with open(summary_file, 'rb') as f:
        tree_distribution_summary = pickle.load(f)
    logger.info(f"Loaded tree distribution summary: {len(tree_distribution_summary['tree_structure'])} trees")
    
    # Load full tree distribution data
    full_file = aggregation_dir / f'{method}_bootstrap_aggregation.pkl'
    with open(full_file, 'rb') as f:
        tree_distribution_full = pickle.load(f)
    logger.info(f"Loaded full tree distribution data: {len(tree_distribution_full['tree_structure'])} trees")
    
    return tree_distribution_summary, tree_distribution_full


def load_tissue_data_from_ssm(ssm_file: Path, logger: logging.Logger) -> Tuple[pd.DataFrame, Dict, List, List, Dict]:
    """
    Load and process tissue mutation data from SSM file format.
    
    This function processes SSM files to extract mutation data and create the necessary
    mappings for both optimization functions and marker validation.
    
    Args:
        ssm_file: Path to SSM file containing tissue mutation data
        logger: Logger instance for reporting
        
    Returns:
        Tuple of (tissue_dataframe, gene2idx_mapping, gene_name_list, gene_list, gene_name2idx_mapping)
    """
    logger.info(f"Loading tissue data from SSM file: {ssm_file}")
    
    # Read SSM file
    try:
        ssm_df = pd.read_csv(ssm_file, sep='\t')
        logger.info(f"SSM file loaded: {ssm_df.shape[0]} mutations, {ssm_df.shape[1]} columns")
    except Exception as e:
        logger.error(f"Error reading SSM file: {e}")
        raise
    
    # Validate SSM file format
    required_columns = ['id', 'gene', 'a', 'd', 'mu_r', 'mu_v']
    missing_columns = [col for col in required_columns if col not in ssm_df.columns]
    if missing_columns:
        logger.error(f"Missing required columns in SSM file: {missing_columns}")
        raise ValueError(f"Invalid SSM file format. Missing columns: {missing_columns}")
    
    # Process tissue data
    tissue_df = ssm_df.copy()
    
    # Create gene list and mapping from tissue data
    gene_list = [f's{idx}' for idx in range(len(tissue_df))]  # Format as 's0', 's1', etc.
    
    # Process gene names with duplicate handling (matching original logic)
    gene_name_list = []
    gene_count = {}
    
    for i in range(tissue_df.shape[0]):
        gene = tissue_df.iloc[i]['gene']  # Get gene name
        
        # Handle duplicate gene names
        if gene in gene_name_list:
            gene_count[gene] += 1
            gene = f"{gene}_{gene_count[gene]}"
        else:
            gene_count[gene] = 1
        
        # Handle non-string gene names (use genomic location if available)
        if not isinstance(gene, str):
            # Use ID as fallback
            gene = f"unknown_{tissue_df.iloc[i]['id']}"
        
        gene_name_list.append(gene)
    
    # Create DUAL mappings to support both use cases:
    # 1. gene2idx: Maps mutation IDs ('s0', 's1'...) to indices - for optimization functions
    # 2. gene_name2idx: Maps actual gene names to indices - for fixed marker validation
    gene2idx = {f's{idx}': idx for idx in range(len(gene_name_list))}
    gene_name2idx = {gene_name: idx for idx, gene_name in enumerate(gene_name_list)}
    
    logger.info(f"Created gene list: {len(gene_list)} genes in 's0', 's1'... format")
    logger.info(f"Processed gene names: {len(gene_name_list)} names with duplicate handling")
    logger.info(f"Created dual mappings: gene2idx for optimization, gene_name2idx for validation")
    
    return tissue_df, gene2idx, gene_name_list, gene_list, gene_name2idx


def load_longitudinal_data_from_csv(csv_file: Path, logger: logging.Logger) -> Dict[str, pd.DataFrame]:
    """
    Load and process longitudinal ddPCR data from CSV format.
    
    This function processes CSV files containing ddPCR measurements across multiple timepoints
    and formats them for Bayesian tree updating.
    
    Args:
        csv_file: Path to CSV file containing longitudinal ddPCR data
        logger: Logger instance for reporting
        
    Returns:
        Dictionary mapping timepoint dates to ddPCR DataFrames
    """
    logger.info(f"Loading longitudinal data from CSV file: {csv_file}")
    
    try:
        # Read CSV file
        longitudinal_df = pd.read_csv(csv_file)
        logger.info(f"Longitudinal CSV loaded: {longitudinal_df.shape[0]} rows, {longitudinal_df.shape[1]} columns")
        
        # Validate required columns
        required_columns = ['date', 'gene', 'mutant_droplets', 'total_droplets']
        missing_columns = [col for col in required_columns if col not in longitudinal_df.columns]
        if missing_columns:
            logger.error(f"Missing required columns in longitudinal CSV: {missing_columns}")
            raise ValueError(f"Invalid CSV format. Missing columns: {missing_columns}")
        
        # Group data by date/timepoint
        timepoint_data = {}
        unique_dates = sorted(longitudinal_df['date'].unique())
        logger.info(f"Found {len(unique_dates)} unique timepoints: {unique_dates}")
        
        for date in unique_dates:
            # Filter data for this timepoint
            timepoint_df = longitudinal_df[longitudinal_df['date'] == date].copy()
            
            # Create ddPCR-compatible format
            ddpcr_df = timepoint_df.set_index('gene')
            ddpcr_df['MutDOR'] = ddpcr_df['mutant_droplets']  # Mutant droplet count
            ddpcr_df['DOR'] = ddpcr_df['total_droplets']      # Total droplet count
            
            timepoint_data[date] = ddpcr_df
            logger.info(f"Processed timepoint {date}: {len(ddpcr_df)} markers")
        
        return timepoint_data
        
    except Exception as e:
        logger.error(f"Error processing longitudinal CSV file: {e}")
        raise