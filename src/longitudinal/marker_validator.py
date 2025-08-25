#!/usr/bin/env python3
"""
Fixed marker validation for longitudinal cancer evolution analysis.

This module handles validation of user-specified fixed markers against the
available dataset and converts them to the internal format used by the pipeline.

Authors: TracerX Pipeline Development Team
"""

import logging
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)


def validate_fixed_markers(user_specified_markers: List[str], gene_name2idx: Dict, gene_name_list: List[str],
                          logger: logging.Logger) -> Tuple[List[str], List[str]]:
    """
    Validate and convert user-specified fixed markers to internal format.
    
    This function checks that user-specified marker names exist in the dataset
    and converts them to the internal marker ID format used by the pipeline.
    
    Args:
        user_specified_markers: List of gene names specified by user
        gene_name2idx: Mapping from gene names to indices
        gene_name_list: List of all available gene names
        logger: Logger instance
        
    Returns:
        Tuple of (validated_marker_ids, validated_gene_names)
    """
    logger.info(f"Validating user-specified fixed markers: {user_specified_markers}")
    
    validated_gene_names = []
    validated_marker_ids = []
    missing_markers = []
    
    for gene_name in user_specified_markers:
        if gene_name in gene_name2idx:
            # Gene found in the dataset
            gene_idx = gene_name2idx[gene_name]
            marker_id = f"s{gene_idx}"  # Convert to marker ID format
            validated_gene_names.append(gene_name)
            validated_marker_ids.append(marker_id)
            logger.info(f"✓ Marker {gene_name} found (index: {gene_idx}, ID: {marker_id})")
        else:
            # Gene not found in dataset
            missing_markers.append(gene_name)
            logger.warning(f"✗ Marker {gene_name} not found in dataset")
    
    # Report validation results
    if validated_marker_ids:
        logger.info(f"Successfully validated {len(validated_marker_ids)} markers: {validated_gene_names}")
    
    if missing_markers:
        logger.warning(f"Missing {len(missing_markers)} markers from dataset: {missing_markers}")
        logger.info("Available gene names (first 10):")
        for i, gene_name in enumerate(gene_name_list[:10]):
            logger.info(f"  {i}: {gene_name}")
        if len(gene_name_list) > 10:
            logger.info(f"  ... and {len(gene_name_list) - 10} more")
    
    if not validated_marker_ids:
        raise ValueError("No valid markers found. Please check marker names against available dataset.")
    
    return validated_marker_ids, validated_gene_names