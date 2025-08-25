#!/usr/bin/env python3
"""
Clone frequency computation for longitudinal cancer evolution analysis.

This module computes clone-level frequencies from tree structures and VAF data,
implementing the root-correction approach as specified in the plan.
No copy-number adjustments are applied as per requirements.

Authors: TracerX Pipeline Development Team
"""

import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional

logger = logging.getLogger(__name__)


def compute_clone_frequencies(tree_distribution: Dict, vaf_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute clone frequencies from tree structure and VAF data.
    
    This function creates a tidy DataFrame with clone frequencies over time
    by aggregating mutation VAFs according to the tree structure.
    
    Args:
        tree_distribution: Tree distribution data from convergence
        vaf_df: DataFrame with VAF data [sample, time, mutation, vaf]
        
    Returns:
        DataFrame with columns [sample, time, clone_id, freq]
    """
    logger.info("Computing clone frequencies from tree structure and VAF data")
    
    # Get the best tree (highest frequency) for clone frequency calculation
    tree_frequencies = tree_distribution['freq']
    best_tree_idx = np.argmax(tree_frequencies)
    best_frequency = tree_frequencies[best_tree_idx]
    
    logger.info(f"Using best tree (index {best_tree_idx}, frequency {best_frequency:.3f}) for clone frequency calculation")
    
    # Extract tree components for the best tree
    tree_structure = tree_distribution['tree_structure'][best_tree_idx]
    node_dict = tree_distribution['node_dict'][best_tree_idx]
    node_dict_name = tree_distribution['node_dict_name'][best_tree_idx]
    
    # Create mutation to clone mapping
    mut_to_clone = {}
    for clone_id, mutations in node_dict_name.items():
        for mutation in mutations:
            mut_to_clone[mutation] = clone_id
    
    # Get unique samples and timepoints from VAF data
    samples = vaf_df['sample'].unique()
    timepoints = vaf_df['time'].unique()
    
    # Compute clone frequencies for each sample and timepoint
    clone_freq_data = []
    
    for sample in samples:
        for timepoint in timepoints:
            # Filter VAF data for this sample and timepoint
            sample_time_vaf = vaf_df[(vaf_df['sample'] == sample) & 
                                   (vaf_df['time'] == timepoint)]
            
            if sample_time_vaf.empty:
                continue
                
            # Compute frequency for each clone
            for clone_id in node_dict_name.keys():
                # Get mutations assigned to this clone
                clone_mutations = node_dict_name[clone_id]
                
                # Find VAF values for mutations in this clone
                clone_vafs = []
                for mutation in clone_mutations:
                    mut_vaf = sample_time_vaf[sample_time_vaf['mutation'] == mutation]
                    if not mut_vaf.empty:
                        clone_vafs.append(mut_vaf['vaf'].iloc[0])
                
                # Compute raw clone frequency as mean VAF of mutations in clone
                if clone_vafs:
                    clone_freq = np.mean(clone_vafs)
                else:
                    clone_freq = 0.0  # No mutations found for this clone
                
                clone_freq_data.append({
                    'sample': sample,
                    'time': timepoint,
                    'clone_id': clone_id,
                    'freq': clone_freq
                })
    
    # Create DataFrame
    freq_df = pd.DataFrame(clone_freq_data)
    
    logger.info(f"Computed clone frequencies for {len(samples)} samples, {len(timepoints)} timepoints, {len(node_dict_name)} clones")
    
    return freq_df


def compute_subtree_remainder(freq_df: pd.DataFrame, tree_distribution: Dict) -> pd.DataFrame:
    """
    Compute and append subtree remainder frequencies.
    
    This implements the root-correction approach:
    1. Calculate total unrooted subtree frequency
    2. Calculate sum of leaf clone frequencies  
    3. Compute remainder as difference
    4. Add remainder as pseudo-clone
    
    Args:
        freq_df: DataFrame with clone frequencies
        tree_distribution: Tree distribution data
        
    Returns:
        DataFrame with remainder pseudo-clone appended
    """
    logger.info("Computing subtree remainder frequencies for root correction")
    
    # Get the best tree structure for remainder calculation
    tree_frequencies = tree_distribution['freq']
    best_tree_idx = np.argmax(tree_frequencies)
    tree_structure = tree_distribution['tree_structure'][best_tree_idx]
    node_dict_name = tree_distribution['node_dict_name'][best_tree_idx]
    
    # Identify leaf clones (nodes with no children)
    leaf_clones = []
    for node_id in node_dict_name.keys():
        if node_id not in tree_structure:  # Node has no children
            leaf_clones.append(node_id)
    
    logger.info(f"Identified {len(leaf_clones)} leaf clones: {leaf_clones}")
    
    # For each sample and timepoint, compute remainder
    remainder_data = []
    
    for sample in freq_df['sample'].unique():
        for timepoint in freq_df['time'].unique():
            # Filter frequencies for this sample and timepoint
            sample_time_freq = freq_df[(freq_df['sample'] == sample) & 
                                     (freq_df['time'] == timepoint)]
            
            if sample_time_freq.empty:
                continue
            
            # Calculate total subtree frequency (sum of all clone frequencies)
            total_subtree_freq = sample_time_freq['freq'].sum()
            
            # Calculate sum of leaf clone frequencies
            leaf_freq_sum = sample_time_freq[sample_time_freq['clone_id'].isin(leaf_clones)]['freq'].sum()
            
            # Compute remainder
            remainder_freq = total_subtree_freq - leaf_freq_sum
            
            # Ensure remainder is non-negative
            remainder_freq = max(0.0, remainder_freq)
            
            remainder_data.append({
                'sample': sample,
                'time': timepoint,
                'clone_id': 'unrooted_remainder',
                'freq': remainder_freq
            })
    
    # Append remainder data to original DataFrame
    remainder_df = pd.DataFrame(remainder_data)
    result_df = pd.concat([freq_df, remainder_df], ignore_index=True)
    
    logger.info(f"Added remainder pseudo-clone for {len(remainder_data)} sample-timepoint combinations")
    
    return result_df


def create_mutation_to_clone_mapping(tree_distribution: Dict) -> Dict[str, int]:
    """
    Create mapping from mutation names to clone IDs.
    
    Args:
        tree_distribution: Tree distribution data
        
    Returns:
        Dictionary mapping mutation names to clone IDs
    """
    # Use the best tree for mapping
    tree_frequencies = tree_distribution['freq']
    best_tree_idx = np.argmax(tree_frequencies)
    node_dict_name = tree_distribution['node_dict_name'][best_tree_idx]
    
    mut_to_clone = {}
    for clone_id, mutations in node_dict_name.items():
        for mutation in mutations:
            mut_to_clone[mutation] = clone_id
    
    return mut_to_clone


def validate_clone_frequency_data(freq_df: pd.DataFrame) -> bool:
    """
    Validate clone frequency data for consistency.
    
    Args:
        freq_df: DataFrame with clone frequencies
        
    Returns:
        True if data is valid, False otherwise
    """
    required_columns = ['sample', 'time', 'clone_id', 'freq']
    
    # Check required columns
    for col in required_columns:
        if col not in freq_df.columns:
            logger.error(f"Missing required column: {col}")
            return False
    
    # Check for negative frequencies
    if (freq_df['freq'] < 0).any():
        logger.error("Found negative clone frequencies")
        return False
    
    # Check for reasonable frequency ranges (0-1)
    if (freq_df['freq'] > 1.0).any():
        logger.warning("Found clone frequencies > 1.0 (may be normal for VAF data)")
    
    logger.info("Clone frequency data validation passed")
    return True