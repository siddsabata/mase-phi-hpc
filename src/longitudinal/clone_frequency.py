#!/usr/bin/env python3
"""
Clone frequency computation for longitudinal cancer evolution analysis.

Implements the paper's algorithm (Fu et al.):
    "We infer the frequency of a clone using the mean VAF of markers appearing
    in the given clone minus the sum of mean VAFs of markers inferred for its
    child clones."

This requires bottom-up (post-order) traversal: process children before parents.
"""

import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Set
from collections import deque

logger = logging.getLogger(__name__)


def get_all_nodes(tree_structure: Dict) -> Set[int]:
    """Get all nodes in the tree (parents and children)."""
    nodes = set(tree_structure.keys())
    for children in tree_structure.values():
        nodes.update(children)
    return nodes


def find_root_node(tree_structure: Dict) -> int:
    """
    Find the root node (node that is parent but not a child of any other node).
    In PhyloWGS trees, this is typically node 0.
    """
    all_children = set()
    for children in tree_structure.values():
        all_children.update(children)

    parents = set(tree_structure.keys())
    roots = parents - all_children

    if len(roots) == 0:
        # If no explicit root, assume 0 is root
        return 0
    return min(roots)  # Return smallest root (typically 0)


def get_children(tree_structure: Dict, node_id: int) -> List[int]:
    """Get direct children of a node. Returns empty list if leaf."""
    return tree_structure.get(node_id, [])


def get_post_order_traversal(tree_structure: Dict, node_dict_name: Dict) -> List[int]:
    """
    Get nodes in post-order (leaves first, root last).
    This ensures children are processed before parents.

    Uses iterative approach to avoid recursion depth issues.
    """
    if not tree_structure and not node_dict_name:
        return []

    # Get all nodes from both tree_structure and node_dict_name
    all_nodes = get_all_nodes(tree_structure)
    all_nodes.update(node_dict_name.keys())

    if not all_nodes:
        return []

    root = find_root_node(tree_structure)

    # Iterative post-order traversal
    result = []
    visited = set()
    stack = [(root, False)]

    while stack:
        node, processed = stack.pop()

        if processed:
            result.append(node)
            continue

        if node in visited:
            continue

        visited.add(node)
        stack.append((node, True))

        # Add children to stack (will be processed before parent)
        children = get_children(tree_structure, node)
        for child in reversed(children):
            if child not in visited:
                stack.append((child, False))

    # Add any nodes not reachable from root (isolated nodes)
    for node in all_nodes:
        if node not in visited:
            result.append(node)

    return result


def compute_clone_frequencies_from_ddpcr(
    tree_distribution: Dict,
    ddpcr_data: pd.DataFrame,
    timepoint: str
) -> Dict[int, float]:
    """
    Compute clone frequencies from ddPCR data using the paper's algorithm.

    Paper's formula:
        clone_freq = mean(VAF of markers in clone) - sum(clone_freq of children)

    CRITICAL: Must be computed bottom-up (children before parents).

    Args:
        tree_distribution: Tree distribution with tree_structure, node_dict_name, freq
        ddpcr_data: DataFrame indexed by gene name with 'MutDOR' and 'DOR' columns
        timepoint: Current timepoint identifier (for logging)

    Returns:
        Dict mapping clone_id -> clone_frequency
    """
    # Get best tree (highest frequency)
    tree_frequencies = tree_distribution['freq']
    best_tree_idx = np.argmax(tree_frequencies)
    best_freq = tree_frequencies[best_tree_idx]

    tree_structure = tree_distribution['tree_structure'][best_tree_idx]
    node_dict_name = tree_distribution['node_dict_name'][best_tree_idx]

    logger.info(f"Computing clone frequencies for timepoint {timepoint}")
    logger.info(f"Using best tree (index {best_tree_idx}, frequency {best_freq:.2f}%)")

    # Get post-order traversal (leaves first)
    post_order = get_post_order_traversal(tree_structure, node_dict_name)
    logger.debug(f"Post-order traversal: {post_order}")

    # Compute VAF from ddPCR for each marker
    available_markers = set(ddpcr_data.index)

    # Dictionary to store computed frequencies
    computed_freqs = {}

    for node_id in post_order:
        # Get markers (mutations) assigned to this clone
        clone_mutations = node_dict_name.get(node_id, [])

        # Calculate VAF for each marker in this clone that's in ddPCR data
        vafs = []
        for mutation in clone_mutations:
            if mutation in available_markers:
                mut_count = ddpcr_data.loc[mutation, 'MutDOR']
                total_count = ddpcr_data.loc[mutation, 'DOR']
                if total_count > 0:
                    vaf = mut_count / total_count
                    vafs.append(vaf)

        # Raw VAF (cellular prevalence proxy) = mean of marker VAFs
        if vafs:
            raw_vaf = np.mean(vafs)
        else:
            raw_vaf = 0.0  # No markers for this clone in ddPCR data

        # CRITICAL: Subtract sum of child clone frequencies
        children = get_children(tree_structure, node_id)
        child_freq_sum = sum(computed_freqs.get(child, 0.0) for child in children)

        # Clone frequency = raw VAF - children's frequencies
        # NOTE: Can be negative (per paper - indicates measurement imprecision)
        clone_freq = raw_vaf - child_freq_sum

        computed_freqs[node_id] = clone_freq

        logger.debug(
            f"Node {node_id}: markers={len(clone_mutations)}, "
            f"measured={len(vafs)}, raw_vaf={raw_vaf:.4f}, "
            f"child_sum={child_freq_sum:.4f}, freq={clone_freq:.4f}"
        )

    return computed_freqs


def compute_clone_frequencies(
    tree_distribution: Dict,
    vaf_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Compute clone frequencies from tree structure and VAF data.

    This version accepts a DataFrame with columns [sample, time, mutation, vaf]
    and returns a tidy DataFrame with clone frequencies.

    Uses the paper's algorithm:
        clone_freq = mean(VAF of markers in clone) - sum(clone_freq of children)

    Args:
        tree_distribution: Tree distribution data from convergence
        vaf_df: DataFrame with VAF data [sample, time, mutation, vaf]

    Returns:
        DataFrame with columns [sample, time, clone_id, freq, raw_vaf, child_sum]
    """
    logger.info("Computing clone frequencies from tree structure and VAF data")

    # Get the best tree (highest frequency)
    tree_frequencies = tree_distribution['freq']
    best_tree_idx = np.argmax(tree_frequencies)
    best_frequency = tree_frequencies[best_tree_idx]

    logger.info(f"Using best tree (index {best_tree_idx}, frequency {best_frequency:.3f})")

    # Extract tree components for the best tree
    tree_structure = tree_distribution['tree_structure'][best_tree_idx]
    node_dict_name = tree_distribution['node_dict_name'][best_tree_idx]

    # Get post-order traversal (leaves first, root last)
    post_order = get_post_order_traversal(tree_structure, node_dict_name)

    # Get unique samples and timepoints from VAF data
    samples = vaf_df['sample'].unique()
    timepoints = vaf_df['time'].unique()

    # Compute clone frequencies for each sample and timepoint
    clone_freq_data = []

    for sample in samples:
        for timepoint in timepoints:
            # Filter VAF data for this sample and timepoint
            sample_time_vaf = vaf_df[
                (vaf_df['sample'] == sample) & (vaf_df['time'] == timepoint)
            ]

            if sample_time_vaf.empty:
                continue

            # Create mutation -> VAF lookup for this sample/timepoint
            mut_to_vaf = dict(zip(sample_time_vaf['mutation'], sample_time_vaf['vaf']))

            # Dictionary to store computed frequencies for this sample/timepoint
            computed_freqs = {}

            # Process in post-order (children before parents)
            for node_id in post_order:
                # Get mutations assigned to this clone
                clone_mutations = node_dict_name.get(node_id, [])

                # Find VAF values for mutations in this clone
                vafs = [mut_to_vaf[m] for m in clone_mutations if m in mut_to_vaf]

                # Raw VAF = mean of marker VAFs
                raw_vaf = np.mean(vafs) if vafs else 0.0

                # CRITICAL: Subtract child frequencies (paper's algorithm)
                children = get_children(tree_structure, node_id)
                child_freq_sum = sum(computed_freqs.get(c, 0.0) for c in children)

                # Clone frequency (can be negative per paper)
                clone_freq = raw_vaf - child_freq_sum
                computed_freqs[node_id] = clone_freq

                clone_freq_data.append({
                    'sample': sample,
                    'time': timepoint,
                    'clone_id': node_id,
                    'freq': clone_freq,
                    'raw_vaf': raw_vaf,
                    'child_sum': child_freq_sum,
                    'n_markers': len(clone_mutations),
                    'n_measured': len(vafs)
                })

    # Create DataFrame
    freq_df = pd.DataFrame(clone_freq_data)

    n_clones = len(node_dict_name)
    logger.info(
        f"Computed clone frequencies for {len(samples)} samples, "
        f"{len(timepoints)} timepoints, {n_clones} clones"
    )

    return freq_df


def convert_ddpcr_to_vaf_df(
    ddpcr_data: pd.DataFrame,
    timepoint: str,
    sample_name: str = 'liquid_biopsy'
) -> pd.DataFrame:
    """
    Convert ddPCR measurements to VAF DataFrame format.

    Args:
        ddpcr_data: DataFrame indexed by gene name with 'MutDOR' and 'DOR' columns
        timepoint: Timepoint identifier
        sample_name: Sample name (default 'liquid_biopsy')

    Returns:
        DataFrame with columns [sample, time, mutation, vaf]
    """
    records = []

    for gene_name in ddpcr_data.index:
        mut_count = ddpcr_data.loc[gene_name, 'MutDOR']
        total_count = ddpcr_data.loc[gene_name, 'DOR']

        vaf = mut_count / total_count if total_count > 0 else 0.0

        records.append({
            'sample': sample_name,
            'time': timepoint,
            'mutation': gene_name,
            'vaf': vaf
        })

    return pd.DataFrame(records)


def validate_clone_frequency_data(freq_df: pd.DataFrame) -> bool:
    """
    Validate clone frequency data for consistency.

    NOTE: Negative frequencies ARE allowed per paper (measurement imprecision).

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

    # Check for NaN values
    if freq_df['freq'].isna().any():
        logger.error("Found NaN clone frequencies")
        return False

    # Log warning for very negative values (likely error, but allowed)
    very_negative = freq_df[freq_df['freq'] < -0.1]
    if len(very_negative) > 0:
        logger.warning(
            f"Found {len(very_negative)} clone frequencies < -0.1 "
            f"(may indicate measurement imprecision per paper)"
        )

    # Log warning for frequencies > 1 (unusual but possible with noise)
    if (freq_df['freq'] > 1.0).any():
        logger.warning("Found clone frequencies > 1.0 (may indicate data issues)")

    logger.info("Clone frequency data validation passed")
    return True


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
