#!/usr/bin/env python3
"""
Dynamic marker analysis for longitudinal cancer evolution.

This module implements dynamic marker analysis where optimal markers are 
selected at each timepoint based on the current tree distribution.

Authors: TracerX Pipeline Development Team
"""

import logging
import pickle
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple
from optimize_fraction import select_markers_tree_gp
from tree_updater import process_ddpcr_measurements, update_tree_distribution, prepare_tree_components_for_analysis
from longitudinal_visualizer import create_visualization_plots, save_visualization_summary

logger = logging.getLogger(__name__)


def run_dynamic_marker_analysis(args, logger: logging.Logger, tree_distribution_summary: Dict, 
                               tree_distribution_full: Dict, gene_list: List[str], gene2idx: Dict, 
                               gene_name_list: List[str], timepoint_data: Dict, 
                               output_dir: Path) -> Dict:
    """
    Run complete dynamic marker analysis across all timepoints.
    
    This approach selects optimal markers at each timepoint based on the current
    tree distribution and uses them to update tree distributions with Bayesian inference.
    
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
        
    Returns:
        Dictionary containing analysis results summary
    """
    logger.info("Starting dynamic marker selection analysis")
    
    # Create output subdirectories
    dynamic_trees_dir = output_dir / 'dynamic_marker_analysis' / 'updated_trees'
    dynamic_selections_dir = output_dir / 'dynamic_marker_analysis' / 'marker_selections'
    dynamic_trees_dir.mkdir(parents=True, exist_ok=True)
    dynamic_selections_dir.mkdir(parents=True, exist_ok=True)
    
    # Sort timepoints chronologically
    sorted_timepoints = sorted(timepoint_data.keys())
    logger.info(f"Processing {len(sorted_timepoints)} timepoints with dynamic marker selection")
    
    # Initialize tracking
    all_marker_selections = []
    analysis_summary = {
        'patient_id': args.patient_id,
        'analysis_mode': 'dynamic',
        'timepoints_processed': []
    }
    
    # Track current tree distribution (starts with original)
    current_tree_summary = tree_distribution_summary
    current_tree_full = tree_distribution_full
    
    # Process each timepoint sequentially
    for order_idx, timepoint in enumerate(sorted_timepoints):
        logger.info(f"Processing timepoint {order_idx + 1}/{len(sorted_timepoints)}: {timepoint}")
        
        # For the first timepoint, use original tree distributions
        if order_idx == 0:
            # Extract tree components for analysis
            tree_list, node_list, tree_freq_list, clonal_freq_list = \
                prepare_tree_components_for_analysis(current_tree_summary, logger)
            
            logger.info(f"Using original tree distributions: {len(tree_list)} trees")
        
        # For subsequent timepoints, load updated tree distribution from previous iteration
        else:
            updated_file = dynamic_trees_dir / f'phylowgs_bootstrap_summary_updated_timepoint_{order_idx-1}.pkl'
            
            if not updated_file.exists():
                logger.error(f"Updated tree distribution file not found: {updated_file}")
                raise FileNotFoundError(f"Missing updated tree distribution for timepoint {order_idx-1}")
            
            with open(updated_file, 'rb') as f:
                current_tree_summary = pickle.load(f)
            
            # Extract tree components
            tree_list, node_list, tree_freq_list, clonal_freq_list = \
                prepare_tree_components_for_analysis(current_tree_summary, logger)
            
            logger.info(f"Using updated tree distributions from previous timepoint: {len(tree_list)} trees")
        
        # Select optimal markers based on current tree structure
        logger.info(f"Selecting {args.n_markers} optimal markers for timepoint {timepoint}")
        
        try:
            selected_markers, obj_frac, obj_struct = select_markers_tree_gp(
                gene_list, args.n_markers, tree_list, node_list, clonal_freq_list,
                gene2idx, tree_freq_list, read_depth=args.read_depth,
                lam1=args.lambda1, lam2=args.lambda2, focus_sample_idx=args.focus_sample)
        except Exception as e:
            logger.error(f"Error in marker selection: {e}")
            import traceback
            logger.error(f"Full traceback:\n{traceback.format_exc()}")
            raise
        
        # Convert selected marker IDs to gene names
        selected_gene_names = [gene_name_list[int(marker_id[1:])] for marker_id in selected_markers]
        
        logger.info(f"Selected markers: {selected_markers}")
        logger.info(f"Selected gene names: {selected_gene_names}")
        
        # Get ddPCR data for current timepoint
        current_ddpcr_data = timepoint_data[timepoint]
        
        # Process ddPCR measurements for selected markers
        try:
            ddpcr_measurements, ddpcr_marker_counts, read_depth_list, marker_idx2gene = \
                process_ddpcr_measurements(selected_gene_names, current_ddpcr_data, timepoint, logger)
        except ValueError as e:
            logger.error(f"Failed to process ddPCR measurements for timepoint {timepoint}: {e}")
            continue
        
        # Update tree distributions using Bayesian approach
        updated_tree_summary = update_tree_distribution(
            current_tree_summary, ddpcr_marker_counts, read_depth_list, marker_idx2gene, logger)
        
        # Save updated tree distribution for next iteration
        updated_file = dynamic_trees_dir / f'phylowgs_bootstrap_summary_updated_timepoint_{order_idx}.pkl'
        with open(updated_file, 'wb') as f:
            pickle.dump(updated_tree_summary, f)
        
        logger.info(f"Saved updated tree distribution: {updated_file}")
        
        # Save marker selection results
        marker_selection_results = {
            'timepoint': timepoint,
            'order_idx': order_idx,
            'selected_markers': selected_markers,
            'selected_gene_names': selected_gene_names,
            'ddpcr_measurements': ddpcr_measurements,
            'objective_fraction': obj_frac,
            'objective_structure': obj_struct,
            'parameters': {
                'n_markers': args.n_markers,
                'read_depth': args.read_depth,
                'lambda1': args.lambda1,
                'lambda2': args.lambda2
            }
        }
        
        all_marker_selections.append(marker_selection_results)
        analysis_summary['timepoints_processed'].append(timepoint)
        
        marker_file = dynamic_selections_dir / f'marker_selection_timepoint_{order_idx}.json'
        with open(marker_file, 'w') as f:
            json.dump(marker_selection_results, f, indent=2, default=str)
        
        logger.info(f"Saved marker selection results: {marker_file}")
        
        # Set current tree for next iteration
        current_tree_summary = updated_tree_summary
    
    # Calculate final performance metrics
    logger.info("Calculating dynamic approach performance metrics...")
    
    final_tree_freq = current_tree_summary['freq']
    tree_entropy = -sum([f * __import__('numpy').log(f + 1e-10) for f in final_tree_freq if f > 0])
    dominant_tree_freq = max(final_tree_freq)
    
    all_markers_used = set()
    for selection in all_marker_selections:
        all_markers_used.update(selection['selected_markers'])
    
    analysis_summary.update({
        'convergence_metrics': {
            'final_tree_entropy': tree_entropy,
            'dominant_tree_frequency': dominant_tree_freq,
            'convergence_timepoint': len(sorted_timepoints),
            'n_trees_remaining': sum(1 for f in final_tree_freq if f > 0.01)
        },
        'performance_metrics': {
            'total_unique_markers': len(all_markers_used),
            'total_selections': len(all_marker_selections),
            'markers_per_timepoint': [len(sel['selected_markers']) for sel in all_marker_selections]
        }
    })
    
    # Save final tree distribution
    final_tree_file = dynamic_trees_dir / f'phylowgs_bootstrap_summary_final_dynamic.pkl'
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
        
        # Prepare VAF data from ddPCR measurements across all timepoints
        vaf_data = []
        for selection in all_marker_selections:
            timepoint = selection['timepoint']
            ddpcr_measurements = selection.get('ddpcr_measurements', [])
            
            for measurement in ddpcr_measurements:
                gene_name = measurement['gene']
                mut_count = measurement['mut']
                total_count = mut_count + measurement['WT']
                vaf = mut_count / total_count if total_count > 0 else 0
                vaf_data.append({
                    'sample': 'sample_1',  # Single sample for dynamic analysis
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
                clone_freq_file = write_clone_frequencies(freq_df, output_dir, 'dynamic')
                
                # Generate clone frequency visualizations
                plot_clone_trajectories(freq_df, output_dir, 'dynamic', args.patient_id)
                plot_clone_heatmap(freq_df, output_dir, 'dynamic', args.patient_id)
                
                logger.info(f"Clone frequency tracking completed successfully")
                
            except Exception as e:
                logger.error(f"Error in clone frequency tracking: {e}")
                clone_freq_file = None
        else:
            logger.warning("No VAF data available for clone frequency computation")
    
    # Create results summary
    results_summary = {
        'analysis_summary': analysis_summary,
        'all_marker_selections': all_marker_selections,
        'final_tree_distribution_file': str(final_tree_file),
        'clone_frequency_file': str(clone_freq_file) if clone_freq_file else None
    }
    
    # Save complete analysis results
    dynamic_results_file = output_dir / 'dynamic_marker_analysis' / 'dynamic_marker_results.json'
    with open(dynamic_results_file, 'w') as f:
        json.dump(results_summary, f, indent=2, default=str)
    
    # Generate visualization plots
    logger.info("Generating visualization plots for dynamic marker analysis")
    try:
        plot_files = create_visualization_plots(
            'dynamic', args.patient_id, output_dir, results_summary, timepoint_data, logger)
        
        # Save visualization summary
        viz_dir = output_dir / 'dynamic_marker_analysis' / 'visualizations'
        save_visualization_summary(viz_dir, plot_files, 'dynamic', args.patient_id)
        
        logger.info("Successfully generated all visualization plots")
    except Exception as e:
        logger.warning(f"Failed to generate visualization plots: {e}")
    
    logger.info(f"Dynamic marker analysis completed successfully")
    logger.info(f"Results saved to: {output_dir / 'dynamic_marker_analysis'}")
    
    return results_summary