#!/usr/bin/env python3
"""
Fixed marker analysis for longitudinal cancer evolution.

This module implements fixed marker analysis where user-specified markers 
are tracked consistently across all timepoints.

Authors: TracerX Pipeline Development Team
"""

import logging
import pickle
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple
from marker_validator import validate_fixed_markers
from tree_updater import process_ddpcr_measurements, update_tree_distribution
from longitudinal_visualizer import create_visualization_plots, save_visualization_summary

logger = logging.getLogger(__name__)


def run_fixed_marker_analysis(args, logger: logging.Logger, tree_distribution_summary: Dict, 
                            tree_distribution_full: Dict, gene_list: List[str], gene2idx: Dict, 
                            gene_name_list: List[str], timepoint_data: Dict, 
                            output_dir: Path, gene_name2idx: Dict) -> Dict:
    """
    Run complete fixed marker analysis across all timepoints.
    
    This approach uses constant user-specified markers and tracks their VAF changes
    over time to update tree distributions using Bayesian inference.
    
    Args:
        args: Configuration arguments
        logger: Logger instance
        tree_distribution_summary: Initial tree distribution summary
        tree_distribution_full: Initial full tree distribution data
        gene_list: List of gene IDs (s0, s1, ...)
        gene2idx: Mapping from gene IDs to indices
        gene_name_list: List of actual gene names
        timepoint_data: Dictionary of timepoint ddPCR data
        output_dir: Output directory for results
        gene_name2idx: Mapping from gene names to indices
        
    Returns:
        Dictionary containing analysis results summary
    """
    logger.info("Starting fixed marker analysis")
    
    # Create output subdirectories
    fixed_trees_dir = output_dir / 'fixed_marker_analysis' / 'updated_trees'
    fixed_data_dir = output_dir / 'fixed_marker_analysis' / 'marker_data'
    fixed_trees_dir.mkdir(parents=True, exist_ok=True)
    fixed_data_dir.mkdir(parents=True, exist_ok=True)
    
    # Validate fixed markers
    try:
        fixed_marker_ids, fixed_gene_names = validate_fixed_markers(
            args.fixed_markers, gene_name2idx, gene_name_list, logger)
    except ValueError as e:
        logger.error(f"Fixed marker validation failed: {e}")
        raise
    
    logger.info(f"Validated fixed markers: {fixed_gene_names}")
    logger.info(f"Marker IDs: {fixed_marker_ids}")
    
    # Sort timepoints chronologically
    sorted_timepoints = sorted(timepoint_data.keys())
    logger.info(f"Processing {len(sorted_timepoints)} timepoints with fixed markers")
    
    # Initialize tracking
    analysis_summary = {
        'patient_id': args.patient_id,
        'analysis_mode': 'fixed',
        'fixed_gene_names': fixed_gene_names,
        'fixed_marker_ids': fixed_marker_ids,
        'timepoints_processed': []
    }
    
    # Track current tree distribution (starts with original)
    current_tree_summary = tree_distribution_summary
    
    # Process each timepoint sequentially
    for order_idx, timepoint in enumerate(sorted_timepoints):
        logger.info(f"Processing timepoint {order_idx + 1}/{len(sorted_timepoints)}: {timepoint}")
        
        # Get ddPCR data for current timepoint
        current_ddpcr_data = timepoint_data[timepoint]
        
        # Process ddPCR measurements for fixed markers
        try:
            ddpcr_measurements, ddpcr_marker_counts, read_depth_list, marker_idx2gene = \
                process_ddpcr_measurements(fixed_gene_names, current_ddpcr_data, timepoint, logger)
        except ValueError as e:
            logger.error(f"Failed to process ddPCR measurements for timepoint {timepoint}: {e}")
            continue
        
        # Update tree distributions using Bayesian approach
        updated_tree_summary = update_tree_distribution(
            current_tree_summary, ddpcr_marker_counts, read_depth_list, marker_idx2gene, logger)
        
        # Save updated tree distribution for next iteration
        updated_file = fixed_trees_dir / f'phylowgs_bootstrap_summary_updated_timepoint_{order_idx}.pkl'
        with open(updated_file, 'wb') as f:
            pickle.dump(updated_tree_summary, f)
        
        logger.info(f"Saved updated tree distribution: {updated_file}")
        
        # Save marker data for this timepoint
        marker_data = {
            'timepoint': timepoint,
            'order_idx': order_idx,
            'fixed_markers': fixed_gene_names,
            'ddpcr_measurements': ddpcr_measurements,
            'mutant_counts': ddpcr_marker_counts,
            'read_depths': read_depth_list
        }
        
        marker_file = fixed_data_dir / f'fixed_markers_timepoint_{order_idx}.json'
        with open(marker_file, 'w') as f:
            json.dump(marker_data, f, indent=2, default=str)
        
        logger.info(f"Saved marker data: {marker_file}")
        
        # Update tracking
        analysis_summary['timepoints_processed'].append(timepoint)
        
        # Set current tree for next iteration
        current_tree_summary = updated_tree_summary
    
    # Calculate final metrics
    logger.info("Calculating fixed approach performance metrics...")
    
    final_tree_freq = current_tree_summary['freq']
    tree_entropy = -sum([f * __import__('numpy').log(f + 1e-10) for f in final_tree_freq if f > 0])
    dominant_tree_freq = max(final_tree_freq)
    
    analysis_summary.update({
        'convergence_metrics': {
            'final_tree_entropy': tree_entropy,
            'dominant_tree_frequency': dominant_tree_freq,
            'convergence_timepoint': len(sorted_timepoints),
            'n_trees_remaining': sum(1 for f in final_tree_freq if f > 0.01)
        },
        'performance_metrics': {
            'markers_tracked': len(fixed_gene_names),
            'total_timepoints': len(sorted_timepoints),
            'successful_updates': len(analysis_summary['timepoints_processed'])
        }
    })
    
    # Save final tree distribution
    final_tree_file = fixed_trees_dir / f'phylowgs_bootstrap_summary_final_fixed.pkl'
    with open(final_tree_file, 'wb') as f:
        pickle.dump(current_tree_summary, f)
    
    logger.info(f"Saved final tree distribution: {final_tree_file}")
    
    # Clone frequency tracking (if enabled)
    clone_freq_file = None
    if hasattr(args, 'track_clone_freq') and args.track_clone_freq:
        logger.info("Computing clone frequencies from converged tree")
        
        # Import clone frequency modules
        from clone_frequency import compute_clone_frequencies, compute_subtree_remainder
        from clone_visualizer import plot_clone_trajectories, plot_clone_heatmap
        from output_manager import write_clone_frequencies
        
        # Prepare VAF data from ddPCR measurements
        vaf_data = []
        for timepoint in sorted_timepoints:
            current_ddpcr_data = timepoint_data[timepoint]
            for gene_name in fixed_gene_names:
                if gene_name in current_ddpcr_data.index:
                    mut_count = current_ddpcr_data.loc[gene_name, 'MutDOR']
                    total_count = current_ddpcr_data.loc[gene_name, 'DOR']
                    vaf = mut_count / total_count if total_count > 0 else 0
                    vaf_data.append({
                        'sample': 'sample_1',  # Single sample for fixed analysis
                        'time': timepoint,
                        'mutation': gene_name,
                        'vaf': vaf
                    })
        
        if vaf_data:
            import pandas as pd
            vaf_df = pd.DataFrame(vaf_data)
            
            # Compute clone frequencies
            try:
                freq_df = compute_clone_frequencies(current_tree_summary, vaf_df)
                freq_df = compute_subtree_remainder(freq_df, current_tree_summary)
                
                # Save clone frequencies
                clone_freq_file = write_clone_frequencies(freq_df, output_dir, 'fixed')
                
                # Generate clone frequency visualizations
                plot_clone_trajectories(freq_df, output_dir, 'fixed', args.patient_id)
                plot_clone_heatmap(freq_df, output_dir, 'fixed', args.patient_id)
                
                logger.info(f"Clone frequency tracking completed successfully")
                
            except Exception as e:
                logger.error(f"Error in clone frequency tracking: {e}")
                clone_freq_file = None
        else:
            logger.warning("No VAF data available for clone frequency computation")
    
    # Create results summary
    results_summary = {
        'analysis_summary': analysis_summary,
        'final_tree_distribution_file': str(final_tree_file),
        'clone_frequency_file': str(clone_freq_file) if clone_freq_file else None
    }
    
    # Save complete analysis results
    fixed_results_file = output_dir / 'fixed_marker_analysis' / 'fixed_marker_results.json'
    with open(fixed_results_file, 'w') as f:
        json.dump(results_summary, f, indent=2, default=str)
    
    # Generate visualization plots
    logger.info("Generating visualization plots for fixed marker analysis")
    try:
        plot_files = create_visualization_plots(
            'fixed', args.patient_id, output_dir, results_summary, timepoint_data, logger)
        
        # Save visualization summary
        viz_dir = output_dir / 'fixed_marker_analysis' / 'visualizations'
        save_visualization_summary(viz_dir, plot_files, 'fixed', args.patient_id)
        
        logger.info("Successfully generated all visualization plots")
    except Exception as e:
        logger.warning(f"Failed to generate visualization plots: {e}")
    
    logger.info(f"Fixed marker analysis completed successfully")
    logger.info(f"Results saved to: {output_dir / 'fixed_marker_analysis'}")
    
    return results_summary