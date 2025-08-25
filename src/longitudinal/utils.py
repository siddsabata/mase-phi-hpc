#!/usr/bin/env python3
"""
Utilities for longitudinal cancer evolution analysis.

This module contains shared helper functions and utilities used 
across different components of the longitudinal analysis pipeline.

Authors: TracerX Pipeline Development Team
"""

import logging
import numpy as np
from typing import Dict, List, Tuple, Any

logger = logging.getLogger(__name__)


def calculate_tree_entropy(tree_frequencies: List[float]) -> float:
    """
    Calculate entropy of tree distribution.
    
    Args:
        tree_frequencies: List of tree frequencies
        
    Returns:
        Entropy value
    """
    return -sum([f * np.log(f + 1e-10) for f in tree_frequencies if f > 0])


def get_dominant_tree_frequency(tree_frequencies: List[float]) -> float:
    """
    Get the frequency of the most dominant tree.
    
    Args:
        tree_frequencies: List of tree frequencies
        
    Returns:
        Maximum frequency value
    """
    return max(tree_frequencies) if tree_frequencies else 0.0


def count_significant_trees(tree_frequencies: List[float], threshold: float = 0.01) -> int:
    """
    Count trees with frequency above threshold.
    
    Args:
        tree_frequencies: List of tree frequencies
        threshold: Minimum frequency threshold
        
    Returns:
        Number of trees above threshold
    """
    return sum(1 for f in tree_frequencies if f > threshold)


def validate_timepoint_data(timepoint_data: Dict, required_genes: List[str], 
                          timepoint: str, logger: logging.Logger) -> bool:
    """
    Validate that timepoint data contains required genes.
    
    Args:
        timepoint_data: DataFrame containing ddPCR data for timepoint
        required_genes: List of required gene names
        timepoint: Timepoint identifier for logging
        logger: Logger instance
        
    Returns:
        True if all required genes are present, False otherwise
    """
    missing_genes = []
    for gene in required_genes:
        if gene not in timepoint_data.index:
            missing_genes.append(gene)
    
    if missing_genes:
        logger.warning(f"Missing genes in timepoint {timepoint}: {missing_genes}")
        return False
    
    return True


def safe_division(numerator: float, denominator: float, default: float = 0.0) -> float:
    """
    Perform safe division with default value for zero denominator.
    
    Args:
        numerator: Numerator value
        denominator: Denominator value  
        default: Default value when denominator is zero
        
    Returns:
        Division result or default value
    """
    return numerator / denominator if denominator != 0 else default


def format_gene_list_for_display(gene_list: List[str], max_display: int = 5) -> str:
    """
    Format gene list for display with truncation if too long.
    
    Args:
        gene_list: List of gene names
        max_display: Maximum number of genes to display
        
    Returns:
        Formatted string representation
    """
    if len(gene_list) <= max_display:
        return ', '.join(gene_list)
    else:
        displayed = ', '.join(gene_list[:max_display])
        remaining = len(gene_list) - max_display
        return f"{displayed} (and {remaining} more)"


def extract_numeric_from_marker_id(marker_id: str) -> int:
    """
    Extract numeric index from marker ID (e.g., 's5' -> 5).
    
    Args:
        marker_id: Marker ID in format 'sN' where N is integer
        
    Returns:
        Numeric index
    """
    if marker_id.startswith('s') and marker_id[1:].isdigit():
        return int(marker_id[1:])
    else:
        raise ValueError(f"Invalid marker ID format: {marker_id}")


def create_result_metadata(patient_id: str, analysis_mode: str, 
                          n_timepoints: int, **kwargs) -> Dict[str, Any]:
    """
    Create standardized metadata for analysis results.
    
    Args:
        patient_id: Patient identifier
        analysis_mode: Analysis mode ('dynamic' or 'fixed')
        n_timepoints: Number of timepoints processed
        **kwargs: Additional metadata fields
        
    Returns:
        Dictionary containing result metadata
    """
    from datetime import datetime
    
    metadata = {
        'patient_id': patient_id,
        'analysis_mode': analysis_mode,
        'timestamp': datetime.now().isoformat(),
        'n_timepoints': n_timepoints,
        'pipeline_version': '2.0.0'
    }
    
    # Add any additional metadata
    metadata.update(kwargs)
    
    return metadata


def log_analysis_progress(logger: logging.Logger, current_step: int, 
                         total_steps: int, step_description: str):
    """
    Log analysis progress in a standardized format.
    
    Args:
        logger: Logger instance
        current_step: Current step number (1-based)
        total_steps: Total number of steps
        step_description: Description of current step
    """
    progress_pct = (current_step / total_steps) * 100
    logger.info(f"Progress: Step {current_step}/{total_steps} ({progress_pct:.1f}%) - {step_description}")


def summarize_marker_usage(all_selections: List[Dict], logger: logging.Logger) -> Dict[str, Any]:
    """
    Summarize marker usage across all timepoints for dynamic analysis.
    
    Args:
        all_selections: List of marker selection results
        logger: Logger instance
        
    Returns:
        Dictionary containing usage summary
    """
    if not all_selections:
        return {'total_unique_markers': 0, 'marker_frequency': {}, 'usage_patterns': {}}
    
    # Count marker usage
    marker_counts = {}
    for selection in all_selections:
        for marker in selection.get('selected_markers', []):
            marker_counts[marker] = marker_counts.get(marker, 0) + 1
    
    # Calculate statistics
    total_unique = len(marker_counts)
    most_used_marker = max(marker_counts.items(), key=lambda x: x[1]) if marker_counts else ('none', 0)
    
    usage_summary = {
        'total_unique_markers': total_unique,
        'marker_frequency': marker_counts,
        'most_used_marker': most_used_marker[0],
        'max_usage_count': most_used_marker[1],
        'avg_markers_per_timepoint': np.mean([len(sel.get('selected_markers', [])) for sel in all_selections])
    }
    
    logger.info(f"Marker usage summary: {total_unique} unique markers used")
    logger.info(f"Most frequently used: {most_used_marker[0]} ({most_used_marker[1]} times)")
    
    return usage_summary