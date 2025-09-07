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
                           logger: logging.Logger) -> Tuple[Dict, Dict]:
    """
    Update tree distribution using Bayesian inference with ddPCR measurements.
    
    This function implements the corrected longitudinal updating algorithm from Mase-phi
    that updates tree frequencies AND recalculates clonal frequencies using the proper
    vaf_frac averaging approach.
    
    Args:
        current_tree_summary: Current tree distribution summary
        ddpcr_marker_counts: List of mutant droplet counts
        read_depth_list: List of total droplet counts  
        marker_idx2gene: Mapping from marker index to gene name
        logger: Logger instance
        
    Returns:
        Tuple of (updated_tree_distribution_summary, update_tracking_data)
    """
    logger.info("Updating tree distributions using corrected Bayesian inference...")
    
    # Extract tree information from summary
    tree_list_summary = current_tree_summary['tree_structure']
    node_name_list_summary = current_tree_summary['node_dict_name']
    tree_freq_list_summary = current_tree_summary['freq']
    
    # Store original state for tracking
    original_entropy = -np.sum([f * np.log(f + 1e-10) for f in tree_freq_list_summary if f > 0])
    
    # Update tree distributions using Bayesian approach
    updated_tree_freq_list = adjust_tree_distribution_struct_bayesian(
        tree_list_summary, node_name_list_summary,
        tree_freq_list_summary, read_depth_list,
        ddpcr_marker_counts, marker_idx2gene)
    
    # Create updated tree distribution summary with new frequencies
    updated_tree_distribution_summary = update_tree_distribution_bayesian(
        current_tree_summary, updated_tree_freq_list)
    
    # CORRECTED ALGORITHM: Recalculate clonal frequencies using vaf_frac averaging
    # This is the key missing piece from the original implementation
    clonal_freq_list = []
    for idx in range(len(updated_tree_distribution_summary['vaf_frac'])):
        clonal_freq_dict = updated_tree_distribution_summary['vaf_frac'][idx]
        clonal_freq_dict_new = {}
        for node, freqs in clonal_freq_dict.items():
            # Use the correct averaging method from Mase-phi lines 84-89
            clonal_freq_dict_new[node] = [list(np.array(freqs).mean(axis=0))]
        clonal_freq_list.append(clonal_freq_dict_new)
    
    # Calculate entropy change and other tracking metrics
    updated_entropy = -np.sum([f * np.log(f + 1e-10) for f in updated_tree_freq_list if f > 0])
    entropy_change = updated_entropy - original_entropy
    
    # Create comprehensive tracking data
    update_tracking_data = {
        'tree_frequencies_before': tree_freq_list_summary,
        'tree_frequencies_after': updated_tree_freq_list,
        'entropy_before': original_entropy,
        'entropy_after': updated_entropy,
        'entropy_change': entropy_change,
        'clonal_freq_list': clonal_freq_list,  # Recalculated clonal frequencies
        'ddpcr_measurements': {
            'marker_counts': ddpcr_marker_counts,
            'read_depths': read_depth_list,
            'markers': [marker_idx2gene[i] for i in range(len(ddpcr_marker_counts))]
        }
    }
    
    # Log the update results
    logger.info(f"Tree frequency update completed")
    logger.info(f"Original entropy: {original_entropy:.4f}")
    logger.info(f"Updated entropy: {updated_entropy:.4f}")
    logger.info(f"Entropy change: {entropy_change:.4f}")
    logger.info(f"Recalculated clonal frequencies for {len(clonal_freq_list)} trees")
    
    return updated_tree_distribution_summary, update_tracking_data


def prepare_tree_components_for_marker_selection(tree_distribution_summary: Dict, 
                                               clonal_freq_list: List[Dict],
                                               logger: logging.Logger) -> Tuple[List, List, List, List]:
    """
    Extract and prepare tree components for marker selection analysis.
    
    This function prepares tree data in the format expected by the marker
    selection optimization functions. It now takes pre-calculated clonal
    frequencies from the corrected update process.
    
    Args:
        tree_distribution_summary: Tree distribution summary from aggregation or previous update
        clonal_freq_list: Pre-calculated clonal frequencies from update process
        logger: Logger instance
        
    Returns:
        Tuple of (tree_list, node_list, tree_freq_list, clonal_freq_list)
    """
    # Extract tree components
    tree_list = tree_distribution_summary['tree_structure']
    node_list = tree_distribution_summary['node_dict']
    tree_freq_list = tree_distribution_summary['freq']
    
    # Convert clonal_freq_list to format expected by marker selection
    clonal_freq_list_processed = []
    for clonal_freq_dict in clonal_freq_list:
        clonal_freq_dict_processed = {}
        for node, freqs in clonal_freq_dict.items():
            # Ensure proper integer node keys and extract frequency values
            clonal_freq_dict_processed[int(node)] = freqs[0] if isinstance(freqs[0], list) else freqs
        clonal_freq_list_processed.append(clonal_freq_dict_processed)
    
    logger.info(f"Prepared tree components: {len(tree_list)} trees for marker selection")
    
    return tree_list, node_list, tree_freq_list, clonal_freq_list_processed