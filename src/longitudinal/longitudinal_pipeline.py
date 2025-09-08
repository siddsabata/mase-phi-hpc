#!/usr/bin/env python3
"""
Unified longitudinal cancer evolution analysis pipeline.

This module implements the corrected longitudinal analysis pipeline that:
1. Updates tree frequencies using Bayesian inference with liquid biopsy data
2. Recalculates clonal frequencies using proper vaf_frac averaging
3. Runs marker selection at each timepoint using both optimization objectives
4. Tracks complete evolution in comprehensive JSON format

Authors: TracerX Pipeline Development Team
"""

import logging
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Any

# Import our corrected modules
from tree_updater import (
    process_ddpcr_measurements, 
    update_tree_distribution, 
    prepare_tree_components_for_marker_selection
)
from longitudinal_tracker import LongitudinalTracker, create_marker_selection_results_dict

# Import marker selection functions
from optimize_fraction import select_markers_tree_gp

logger = logging.getLogger(__name__)


def run_unified_longitudinal_analysis(args, logger: logging.Logger, 
                                     tree_distribution_summary: Dict, 
                                     tree_distribution_full: Dict,
                                     gene_list: List[str], gene2idx: Dict, 
                                     gene_name_list: List[str], 
                                     timepoint_data: Dict[str, pd.DataFrame],
                                     output_dir: Path, gene_name2idx: Dict, 
                                     mutation_id_to_gene: Dict) -> Dict:
    """
    Run unified longitudinal analysis pipeline with corrected updating algorithm.
    
    This replaces the separate fixed_analysis.py and dynamic_analysis.py with a single
    unified approach that runs marker selection at each timepoint.
    
    Args:
        args: Configuration arguments
        logger: Logger instance
        tree_distribution_summary: Tree distribution summary from aggregation
        tree_distribution_full: Full tree distribution from aggregation
        gene_list: List of gene IDs
        gene2idx: Mapping from gene IDs to indices
        gene_name_list: List of human-readable gene names
        timepoint_data: Dictionary mapping timepoint names to ddPCR data
        output_dir: Output directory for results
        gene_name2idx: Mapping from gene names to indices
        
    Returns:
        Dictionary containing comprehensive results summary
    """
    logger.info("Starting unified longitudinal analysis pipeline")
    
    # Initialize longitudinal tracker for comprehensive JSON storage
    tracker = LongitudinalTracker(
        patient_id=args.patient_id,
        output_dir=output_dir,
        initial_tree_distribution=tree_distribution_summary
    )
    
    # Create directories for intermediate results
    pipeline_dir = output_dir / 'longitudinal_pipeline'
    pipeline_dir.mkdir(exist_ok=True, parents=True)
    
    # Sort timepoints chronologically
    timepoint_names = list(timepoint_data.keys())
    try:
        # Try to parse as dates for proper sorting
        timepoint_dates = [pd.to_datetime(tp) for tp in timepoint_names]
        sorted_indices = np.argsort(timepoint_dates)
        sorted_timepoints = [timepoint_names[i] for i in sorted_indices]
    except:
        # Fall back to lexicographic sorting if date parsing fails
        sorted_timepoints = sorted(timepoint_names)
        logger.warning("Could not parse timepoints as dates, using lexicographic sorting")
    
    logger.info(f"Processing {len(sorted_timepoints)} timepoints: {sorted_timepoints}")
    
    # Initialize current tree distribution (start with aggregation results)
    current_tree_summary = tree_distribution_summary.copy()
    
    # Track results for each timepoint
    timepoint_results = []
    
    # Process each timepoint
    for order_idx, timepoint in enumerate(sorted_timepoints):
        logger.info(f"\n=== Processing Timepoint {order_idx}: {timepoint} ===")
        
        # Get ddPCR data for this timepoint
        current_ddpcr_data = timepoint_data[timepoint]
        logger.info(f"ddPCR data shape: {current_ddpcr_data.shape}")
        
        try:
            # Step 1: Use markers available in ddPCR CSV file for this timepoint
            logger.info("Step 1: Using markers from ddPCR CSV file for this timepoint")
            
            # Get all available markers from ddPCR data
            selected_markers = sorted(current_ddpcr_data.index.tolist())
            logger.info(f"Found {len(selected_markers)} markers in ddPCR data: {selected_markers}")
            
            # Step 2: Process ddPCR measurements for selected markers
            logger.info("Step 2: Processing ddPCR measurements for selected markers")
            
            ddpcr_measurements, ddpcr_marker_counts, read_depth_list, marker_idx2gene = \
                process_ddpcr_measurements(selected_markers, current_ddpcr_data, timepoint, logger)
            
            # Step 3: Update tree distribution using corrected Bayesian algorithm
            logger.info("Step 3: Updating tree distribution with corrected algorithm")
            
            updated_tree_summary, update_tracking_data = update_tree_distribution(
                current_tree_summary, ddpcr_marker_counts, read_depth_list, marker_idx2gene, logger)
            
            # Step 4: Run marker selection on the UPDATED trees (for potential future use)
            logger.info("Step 4: Running marker selection on updated trees for future timepoints")
            
            # Get clonal frequencies from the update tracking data (corrected algorithm)
            clonal_freq_list = update_tracking_data['clonal_freq_list']
            
            # Prepare tree components for marker selection using updated trees
            tree_list, node_list, tree_freq_list, clonal_freq_list_processed = \
                prepare_tree_components_for_marker_selection(
                    updated_tree_summary, clonal_freq_list, logger)
            
            # Run marker selection with fraction optimization (λ1=1, λ2=0)
            logger.info("Running fraction optimization (λ1=1, λ2=0) on updated trees")
            fraction_markers, fraction_obj_frac, fraction_obj_struct = select_markers_tree_gp(
                gene_list, args.n_markers, tree_list, node_list, clonal_freq_list_processed,
                gene2idx, tree_freq_list, read_depth=args.read_depth, 
                lam1=1, lam2=0, focus_sample_idx=0, mutation_id_to_gene=mutation_id_to_gene)
            
            # Run marker selection with structure optimization (λ1=0, λ2=1)
            logger.info("Running structure optimization (λ1=0, λ2=1) on updated trees")  
            structure_markers, structure_obj_frac, structure_obj_struct = select_markers_tree_gp(
                gene_list, args.n_markers, tree_list, node_list, clonal_freq_list_processed,
                gene2idx, tree_freq_list, read_depth=args.read_depth,
                lam1=0, lam2=1, focus_sample_idx=0, mutation_id_to_gene=mutation_id_to_gene)
            
            # Convert marker IDs to gene names
            fraction_gene_names = [gene_name_list[int(marker[1:])] for marker in fraction_markers]
            structure_gene_names = [gene_name_list[int(marker[1:])] for marker in structure_markers]
            
            logger.info(f"Fraction optimization selected: {fraction_gene_names}")
            logger.info(f"Structure optimization selected: {structure_gene_names}")
            
            # Step 5: Prepare marker selection results for tracking
            marker_selection_results = create_marker_selection_results_dict(
                {
                    'selected_markers': fraction_markers,
                    'obj_frac': fraction_obj_frac,
                    'obj_struct': fraction_obj_struct
                },
                {
                    'selected_markers': structure_markers,
                    'obj_frac': structure_obj_frac,
                    'obj_struct': structure_obj_struct
                },
                "post_update_selection"  # Indicates this selection happened after tree update
            )
            
            # Step 6: Add comprehensive tracking data
            tracker.add_timepoint_update(
                timepoint=timepoint,
                order_idx=order_idx,
                marker_selection_results=marker_selection_results,
                ddpcr_measurements=ddpcr_measurements,
                update_tracking_data=update_tracking_data,
                selected_markers=selected_markers  # The markers used for THIS timepoint's tree update
            )
            
            # Step 7: Save updated tree distribution for next iteration
            updated_tree_file = pipeline_dir / f'tree_distribution_updated_timepoint_{order_idx}.pkl'
            with open(updated_tree_file, 'wb') as f:
                pickle.dump(updated_tree_summary, f)
            
            logger.info(f"Saved updated tree distribution: {updated_tree_file}")
            
            # Step 8: Update clonal frequency list for next iteration
            # (Already obtained from update_tracking_data in Step 4)
            
            # Store results for this timepoint
            timepoint_result = {
                'timepoint': timepoint,
                'order_idx': order_idx,
                'selected_markers': selected_markers,
                'marker_selection_results': marker_selection_results,
                'tree_update_summary': {
                    'entropy_change': update_tracking_data['entropy_change'],
                    'significant_changes': len([i for i, (old, new) in enumerate(
                        zip(update_tracking_data['tree_frequencies_before'], 
                            update_tracking_data['tree_frequencies_after'])) 
                                               if abs(old - new) > 1.0])
                },
                'updated_tree_file': str(updated_tree_file)
            }
            timepoint_results.append(timepoint_result)
            
            # Update current tree summary for next iteration
            current_tree_summary = updated_tree_summary
            
            logger.info(f"Completed timepoint {timepoint} successfully")
            
        except Exception as e:
            logger.error(f"Failed to process timepoint {timepoint}: {e}")
            # Continue with next timepoint rather than failing completely
            continue
    
    # Step 8: Finalize tracking and generate comprehensive results
    logger.info("Finalizing longitudinal analysis results")
    
    final_json_path = tracker.finalize_tracking(current_tree_summary)
    
    # Save final tree distribution
    final_tree_file = pipeline_dir / f'tree_distribution_final.pkl'
    with open(final_tree_file, 'wb') as f:
        pickle.dump(current_tree_summary, f)
    
    # Create comprehensive results summary
    results_summary = {
        'patient_id': args.patient_id,
        'analysis_type': 'unified_longitudinal',
        'total_timepoints_processed': len(timepoint_results),
        'total_timepoints_available': len(sorted_timepoints),
        'timepoint_results': timepoint_results,
        'final_tree_distribution_file': str(final_tree_file),
        'comprehensive_json_file': str(final_json_path),
        'output_directory': str(output_dir),
        'parameters': {
            'n_markers': args.n_markers,
            'read_depth': args.read_depth,
            'method': args.method
        }
    }
    
    logger.info(f"Unified longitudinal analysis completed successfully")
    logger.info(f"Processed {len(timepoint_results)} timepoints")
    logger.info(f"Results saved to: {output_dir}")
    logger.info(f"Comprehensive JSON: {final_json_path}")
    
    return results_summary


def choose_optimal_markers(fraction_results: Dict, structure_results: Dict, 
                          strategy: str = "fraction") -> Tuple[List[str], str]:
    """
    Choose optimal markers from fraction and structure optimization results.
    
    For now, this implements a simple strategy selection, but could be extended
    with more sophisticated selection criteria.
    
    Args:
        fraction_results: Results from fraction optimization
        structure_results: Results from structure optimization  
        strategy: Selection strategy ("fraction", "structure", or "adaptive")
        
    Returns:
        Tuple of (selected_markers, selected_strategy)
    """
    if strategy == "fraction":
        return fraction_results['selected_markers'], "fraction_optimization"
    elif strategy == "structure":
        return structure_results['selected_markers'], "structure_optimization"
    elif strategy == "adaptive":
        # Could implement adaptive selection based on scores, convergence, etc.
        # For now, default to fraction
        return fraction_results['selected_markers'], "fraction_optimization"
    else:
        return fraction_results['selected_markers'], "fraction_optimization"