#!/usr/bin/env python3
"""
Bayesian tree updating for longitudinal cancer evolution analysis.

This module handles the core tree distribution updating using ddPCR measurements
and Bayesian inference to track cancer evolution over time.

Authors: TracerX Pipeline Development Team
"""

import logging
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple
from adjust_tree_distribution import adjust_tree_distribution_struct_bayesian, update_tree_distribution_bayesian

logger = logging.getLogger(__name__)


def process_ddpcr_measurements(selected_gene_names: List[str], timepoint_data: pd.DataFrame, 
                             timepoint: str, logger: logging.Logger) -> Tuple[List[Dict], List[int], List[int], Dict]:
    """
    Process ddPCR measurements for selected markers at a specific timepoint.
    
    Args:
        selected_gene_names: List of gene names to extract measurements for
        timepoint_data: DataFrame containing ddPCR data for this timepoint
        timepoint: Timepoint identifier for logging
        logger: Logger instance
        
    Returns:
        Tuple of (ddpcr_measurements, ddpcr_marker_counts, read_depth_list, marker_idx2gene)
    """
    ddpcr_measurements = []
    
    for gene_name in selected_gene_names:
        if gene_name in timepoint_data.index:
            mut_count = timepoint_data.loc[gene_name, 'MutDOR']  # Mutant droplets
            total_count = timepoint_data.loc[gene_name, 'DOR']   # Total droplets
            wt_count = total_count - mut_count  # Calculate WT count
            
            ddpcr_measurements.append({
                'gene': gene_name,
                'mut': mut_count,
                'WT': wt_count,
                'liquid_biopsy_sample': timepoint
            })
        else:
            logger.warning(f"Selected marker {gene_name} not found in ddPCR data for timepoint {timepoint}")
    
    if not ddpcr_measurements:
        logger.error(f"No ddPCR data found for selected markers at timepoint {timepoint}")
        raise ValueError(f"No ddPCR data available for timepoint {timepoint}")
    
    # Create DataFrame with ddPCR measurements
    df_ddpcr = pd.DataFrame(ddpcr_measurements)
    marker_idx2gene = {i: df_ddpcr["gene"].iloc[i] for i in range(len(df_ddpcr))}
    
    # Extract counts for Bayesian updating
    ddpcr_marker_counts = list(df_ddpcr["mut"])
    read_depth_list = list(df_ddpcr["mut"] + df_ddpcr["WT"])  # Total = mut + WT
    
    logger.info(f"ddPCR measurements: {len(ddpcr_measurements)} markers")
    logger.info(f"Mutant counts: {ddpcr_marker_counts}")
    logger.info(f"Read depths: {read_depth_list}")
    
    return ddpcr_measurements, ddpcr_marker_counts, read_depth_list, marker_idx2gene


def update_tree_distribution(current_tree_summary: Dict, ddpcr_marker_counts: List[int], 
                           read_depth_list: List[int], marker_idx2gene: Dict, 
                           logger: logging.Logger) -> Dict:
    """
    Update tree distribution using Bayesian inference with ddPCR measurements.
    
    This function takes the current tree distribution and updates the tree frequencies
    based on ddPCR measurements using Bayesian inference.
    
    Args:
        current_tree_summary: Current tree distribution summary
        ddpcr_marker_counts: List of mutant droplet counts
        read_depth_list: List of total droplet counts  
        marker_idx2gene: Mapping from marker index to gene name
        logger: Logger instance
        
    Returns:
        Updated tree distribution summary
    """
    logger.info("Updating tree distributions using Bayesian inference...")
    
    # Extract tree information from summary
    tree_list_summary = current_tree_summary['tree_structure']
    node_name_list_summary = current_tree_summary['node_dict_name']
    tree_freq_list_summary = current_tree_summary['freq']
    
    # Update tree distributions using Bayesian approach
    updated_tree_freq_list = adjust_tree_distribution_struct_bayesian(
        tree_list_summary, node_name_list_summary,
        tree_freq_list_summary, read_depth_list,
        ddpcr_marker_counts, marker_idx2gene)
    
    # Create updated tree distribution summary
    updated_tree_distribution_summary = update_tree_distribution_bayesian(
        current_tree_summary, updated_tree_freq_list)
    
    # Log the update results
    original_entropy = -np.sum([f * np.log(f + 1e-10) for f in tree_freq_list_summary if f > 0])
    updated_entropy = -np.sum([f * np.log(f + 1e-10) for f in updated_tree_freq_list if f > 0])
    
    logger.info(f"Tree frequency update completed")
    logger.info(f"Original entropy: {original_entropy:.4f}")
    logger.info(f"Updated entropy: {updated_entropy:.4f}")
    logger.info(f"Entropy change: {updated_entropy - original_entropy:.4f}")
    
    return updated_tree_distribution_summary


def prepare_tree_components_for_analysis(tree_distribution_summary: Dict, 
                                       logger: logging.Logger) -> Tuple[List, List, List, List]:
    """
    Extract and prepare tree components for marker selection analysis.
    
    This function prepares tree data in the format expected by the marker
    selection optimization functions.
    
    Args:
        tree_distribution_summary: Tree distribution summary from aggregation or previous update
        logger: Logger instance
        
    Returns:
        Tuple of (tree_list, node_list, tree_freq_list, clonal_freq_list)
    """
    # Extract tree components
    tree_list = tree_distribution_summary['tree_structure']
    node_list = tree_distribution_summary['node_dict']
    tree_freq_list = tree_distribution_summary['freq']
    
    # Recalculate clonal frequencies by averaging across samples
    clonal_freq_list = []
    for idx in range(len(tree_distribution_summary['vaf_frac'])):
        clonal_freq_dict = tree_distribution_summary['vaf_frac'][idx]
        clonal_freq_dict_new = {}
        for node, freqs in clonal_freq_dict.items():
            clonal_freq_dict_new[node] = [list(np.array(freqs).mean(axis=0))]
        clonal_freq_list.append(clonal_freq_dict_new)
    
    logger.info(f"Prepared tree components: {len(tree_list)} trees")
    
    return tree_list, node_list, tree_freq_list, clonal_freq_list