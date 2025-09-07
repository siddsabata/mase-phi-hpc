#!/usr/bin/env python3
"""
Updated longitudinal visualization module for unified pipeline.

This module creates essential tracking plots for the unified longitudinal
analysis pipeline, focusing on tree convergence, marker evolution, and
clonal frequency tracking without complex heatmaps.

Authors: TracerX Pipeline Development Team  
"""

import logging
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Any, Optional

logger = logging.getLogger(__name__)


def create_unified_visualization_plots(patient_id: str, output_dir: Path, 
                                     comprehensive_json_path: Path,
                                     logger: logging.Logger) -> Dict[str, str]:
    """
    Create essential visualization plots for unified longitudinal analysis.
    
    Args:
        patient_id: Patient identifier
        output_dir: Output directory for plots
        comprehensive_json_path: Path to comprehensive JSON tracking file
        logger: Logger instance
        
    Returns:
        Dictionary mapping plot types to file paths
    """
    logger.info("Creating unified longitudinal visualization plots")
    
    # Create visualization directory
    viz_dir = output_dir / 'visualizations'
    viz_dir.mkdir(parents=True, exist_ok=True)
    
    # Load comprehensive tracking data
    with open(comprehensive_json_path, 'r') as f:
        tracking_data = json.load(f)
    
    plot_files = {}
    
    try:
        # 1. Tree Convergence Plot - Show how tree frequencies evolve over timepoints
        convergence_plot = create_tree_convergence_plot(
            patient_id, viz_dir, tracking_data, logger)
        plot_files['tree_convergence'] = convergence_plot
        
        # 2. Marker Selection Evolution - Track which markers are selected at each timepoint  
        marker_evolution_plot = create_marker_selection_evolution_plot(
            patient_id, viz_dir, tracking_data, logger)
        plot_files['marker_evolution'] = marker_evolution_plot
        
        # 3. Clonal Frequency Evolution - Line plots showing clone frequency changes
        clonal_freq_plot = create_clonal_frequency_evolution_plot(
            patient_id, viz_dir, tracking_data, logger)
        plot_files['clonal_frequency_evolution'] = clonal_freq_plot
        
        # 4. Entropy Evolution - Show convergence metrics over time
        entropy_plot = create_entropy_evolution_plot(
            patient_id, viz_dir, tracking_data, logger)
        plot_files['entropy_evolution'] = entropy_plot
        
        # 5. Save visualization summary
        summary_file = save_visualization_summary(viz_dir, plot_files, tracking_data)
        plot_files['summary'] = summary_file
        
        logger.info(f"Successfully created all unified visualization plots")
        logger.info(f"Plots saved in: {viz_dir}")
        
    except Exception as e:
        logger.error(f"Error creating unified visualization plots: {e}")
        raise
    
    return plot_files


def create_tree_convergence_plot(patient_id: str, viz_dir: Path, 
                                tracking_data: Dict, logger: logging.Logger) -> str:
    """
    Create tree convergence plot showing how tree frequencies change over timepoints.
    
    Args:
        patient_id: Patient identifier
        viz_dir: Visualization directory
        tracking_data: Comprehensive tracking data
        logger: Logger instance
        
    Returns:
        Path to created plot
    """
    logger.info("Creating tree convergence plot")
    
    timepoints = tracking_data['timepoints']
    if not timepoints:
        logger.warning("No timepoint data available for tree convergence plot")
        return ""
    
    # Extract timepoint names and tree frequency data
    timepoint_names = [tp['timepoint'] for tp in timepoints]
    
    # Get initial tree frequencies
    initial_freqs = tracking_data['initial_tree_distribution']['freq']
    num_trees = len(initial_freqs)
    
    # Prepare data for plotting
    tree_freq_data = []
    
    # Add initial state (timepoint 0)
    for tree_idx in range(num_trees):
        tree_freq_data.append({
            'timepoint': 'Initial',
            'tree_idx': f'Tree_{tree_idx}',
            'frequency': initial_freqs[tree_idx]
        })
    
    # Add data for each update timepoint
    for tp_idx, tp_data in enumerate(timepoints):
        frequencies_after = tp_data['tree_update']['frequencies_after']
        for tree_idx in range(min(num_trees, len(frequencies_after))):
            tree_freq_data.append({
                'timepoint': tp_data['timepoint'],
                'tree_idx': f'Tree_{tree_idx}',
                'frequency': frequencies_after[tree_idx]
            })
    
    # Create DataFrame and plot
    df_freq = pd.DataFrame(tree_freq_data)
    
    plt.figure(figsize=(12, 8))
    
    # Create line plot with different colors for each tree
    colors = plt.cm.Set1(np.linspace(0, 1, num_trees))
    
    for tree_idx in range(num_trees):
        tree_data = df_freq[df_freq['tree_idx'] == f'Tree_{tree_idx}']
        plt.plot(range(len(tree_data)), tree_data['frequency'], 
                marker='o', linewidth=2, markersize=6, 
                color=colors[tree_idx], label=f'Tree {tree_idx}')
    
    plt.xlabel('Timepoint', fontsize=12, fontweight='bold')
    plt.ylabel('Tree Frequency', fontsize=12, fontweight='bold') 
    plt.title(f'{patient_id} - Tree Frequency Convergence Over Time', 
              fontsize=14, fontweight='bold')
    
    # Set x-axis labels
    x_labels = ['Initial'] + timepoint_names
    plt.xticks(range(len(x_labels)), x_labels, rotation=45)
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    plot_path = viz_dir / f'{patient_id}_tree_convergence.png'
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved tree convergence plot: {plot_path}")
    return str(plot_path)


def create_marker_selection_evolution_plot(patient_id: str, viz_dir: Path,
                                         tracking_data: Dict, logger: logging.Logger) -> str:
    """
    Create marker selection evolution plot showing which markers are selected at each timepoint.
    
    Args:
        patient_id: Patient identifier
        viz_dir: Visualization directory
        tracking_data: Comprehensive tracking data
        logger: Logger instance
        
    Returns:
        Path to created plot
    """
    logger.info("Creating marker selection evolution plot")
    
    timepoints = tracking_data['timepoints']
    if not timepoints:
        logger.warning("No timepoint data available for marker evolution plot")
        return ""
    
    # Extract marker selection data
    timepoint_names = []
    fraction_markers = []
    structure_markers = []
    selected_markers = []
    
    for tp_data in timepoints:
        timepoint_names.append(tp_data['timepoint'])
        
        # Get markers from both optimization strategies
        marker_sel = tp_data['marker_selection']
        fraction_markers.append(marker_sel['fraction_optimization']['markers'])
        structure_markers.append(marker_sel['structure_optimization']['markers'])
        selected_markers.append(tp_data['selected_markers'])
    
    # Create subplot figure
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 12))
    
    # Plot 1: Selected markers over time
    all_markers = set()
    for markers in selected_markers:
        all_markers.update(markers)
    
    if all_markers:
        all_markers = sorted(list(all_markers))
        selection_matrix = np.zeros((len(all_markers), len(timepoint_names)))
        
        for t_idx, markers in enumerate(selected_markers):
            for marker in markers:
                if marker in all_markers:
                    m_idx = all_markers.index(marker)
                    selection_matrix[m_idx, t_idx] = 1
        
        # Create heatmap for selected markers
        sns.heatmap(selection_matrix, 
                   xticklabels=timepoint_names,
                   yticklabels=all_markers,
                   cmap='Blues', 
                   cbar=True,
                   ax=ax1)
        ax1.set_title(f'{patient_id} - Final Selected Markers Over Time', fontweight='bold')
        ax1.set_xlabel('Timepoint')
        ax1.set_ylabel('Markers')
    
    # Plot 2: Fraction vs Structure marker comparison
    strategy_data = []
    for t_idx, tp_name in enumerate(timepoint_names):
        frac_markers = fraction_markers[t_idx]
        struct_markers = structure_markers[t_idx]
        
        strategy_data.append({
            'timepoint': tp_name,
            'fraction_markers': len(frac_markers),
            'structure_markers': len(struct_markers),
            'overlap': len(set(frac_markers) & set(struct_markers))
        })
    
    df_strategy = pd.DataFrame(strategy_data)
    
    x_pos = np.arange(len(timepoint_names))
    width = 0.25
    
    ax2.bar(x_pos - width, df_strategy['fraction_markers'], width, 
           label='Fraction Optimization', color='orange', alpha=0.8)
    ax2.bar(x_pos, df_strategy['structure_markers'], width,
           label='Structure Optimization', color='green', alpha=0.8)
    ax2.bar(x_pos + width, df_strategy['overlap'], width,
           label='Overlap', color='purple', alpha=0.8)
    
    ax2.set_xlabel('Timepoint')
    ax2.set_ylabel('Number of Markers')
    ax2.set_title('Marker Selection Strategy Comparison', fontweight='bold')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(timepoint_names, rotation=45)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Marker consistency over time
    marker_consistency = []
    for t_idx in range(len(timepoint_names)):
        if t_idx == 0:
            consistency = 1.0  # First timepoint is always 100% consistent with itself
        else:
            prev_markers = set(selected_markers[t_idx-1])
            curr_markers = set(selected_markers[t_idx])
            if prev_markers:
                consistency = len(prev_markers & curr_markers) / len(prev_markers)
            else:
                consistency = 0.0
        marker_consistency.append(consistency)
    
    ax3.plot(range(len(timepoint_names)), marker_consistency, 
            marker='o', linewidth=2, markersize=8, color='red')
    ax3.set_xlabel('Timepoint')
    ax3.set_ylabel('Marker Consistency')
    ax3.set_title('Marker Selection Consistency Over Time', fontweight='bold')
    ax3.set_xticks(range(len(timepoint_names)))
    ax3.set_xticklabels(timepoint_names, rotation=45)
    ax3.set_ylim(0, 1.1)
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    plot_path = viz_dir / f'{patient_id}_marker_evolution.png'
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved marker evolution plot: {plot_path}")
    return str(plot_path)


def create_clonal_frequency_evolution_plot(patient_id: str, viz_dir: Path,
                                         tracking_data: Dict, logger: logging.Logger) -> str:
    """
    Create clonal frequency evolution plot showing how clone frequencies change over time.
    
    Args:
        patient_id: Patient identifier
        viz_dir: Visualization directory
        tracking_data: Comprehensive tracking data
        logger: Logger instance
        
    Returns:
        Path to created plot
    """
    logger.info("Creating clonal frequency evolution plot")
    
    timepoints = tracking_data['timepoints']
    if not timepoints:
        logger.warning("No timepoint data available for clonal frequency plot")
        return ""
    
    # Extract clonal frequency data
    clone_freq_data = []
    
    # Get initial clonal frequencies (calculate from initial tree distribution)
    initial_tree_dist = tracking_data['initial_tree_distribution']
    initial_clones = set()
    for vaf_frac in initial_tree_dist['vaf_frac']:
        initial_clones.update(vaf_frac.keys())
    
    # Add initial state
    for clone in initial_clones:
        clone_freq_data.append({
            'timepoint': 'Initial',
            'clone_id': f'clone_{clone}', 
            'frequency': 0.0  # Would need to calculate from initial data
        })
    
    # Add data for each timepoint
    for tp_data in timepoints:
        clonal_freqs_before = tp_data['clonal_frequencies']['before']
        clonal_freqs_after = tp_data['clonal_frequencies']['after']
        
        # Use 'after' frequencies for tracking evolution
        for clone_id, frequency in clonal_freqs_after.items():
            clone_freq_data.append({
                'timepoint': tp_data['timepoint'],
                'clone_id': clone_id,
                'frequency': frequency
            })
    
    if not clone_freq_data:
        logger.warning("No clonal frequency data to plot")
        return ""
    
    df_clones = pd.DataFrame(clone_freq_data)
    
    plt.figure(figsize=(12, 8))
    
    # Plot each clone as a separate line
    unique_clones = df_clones['clone_id'].unique()
    colors = plt.cm.Set3(np.linspace(0, 1, len(unique_clones)))
    
    for clone_idx, clone_id in enumerate(unique_clones):
        clone_data = df_clones[df_clones['clone_id'] == clone_id]
        plt.plot(range(len(clone_data)), clone_data['frequency'],
                marker='o', linewidth=2, markersize=6,
                color=colors[clone_idx], label=clone_id)
    
    plt.xlabel('Timepoint', fontsize=12, fontweight='bold')
    plt.ylabel('Clonal Frequency', fontsize=12, fontweight='bold')
    plt.title(f'{patient_id} - Clonal Frequency Evolution Over Time', 
              fontsize=14, fontweight='bold')
    
    # Set x-axis labels
    unique_timepoints = df_clones['timepoint'].unique()
    plt.xticks(range(len(unique_timepoints)), unique_timepoints, rotation=45)
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save plot
    plot_path = viz_dir / f'{patient_id}_clonal_frequency_evolution.png'
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved clonal frequency evolution plot: {plot_path}")
    return str(plot_path)


def create_entropy_evolution_plot(patient_id: str, viz_dir: Path,
                                tracking_data: Dict, logger: logging.Logger) -> str:
    """
    Create entropy evolution plot showing convergence metrics over time.
    
    Args:
        patient_id: Patient identifier  
        viz_dir: Visualization directory
        tracking_data: Comprehensive tracking data
        logger: Logger instance
        
    Returns:
        Path to created plot
    """
    logger.info("Creating entropy evolution plot")
    
    timepoints = tracking_data['timepoints']
    if not timepoints:
        logger.warning("No timepoint data available for entropy plot")
        return ""
    
    # Extract entropy data
    timepoint_names = []
    entropy_before = []
    entropy_after = []
    entropy_changes = []
    
    for tp_data in timepoints:
        timepoint_names.append(tp_data['timepoint'])
        entropy_before.append(tp_data['tree_update']['entropy_before'])
        entropy_after.append(tp_data['tree_update']['entropy_after'])
        entropy_changes.append(tp_data['tree_update']['entropy_change'])
    
    # Create subplot figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Plot 1: Entropy evolution
    x_pos = range(len(timepoint_names))
    ax1.plot(x_pos, entropy_before, marker='o', linewidth=2, markersize=6,
            color='blue', label='Entropy Before Update')
    ax1.plot(x_pos, entropy_after, marker='s', linewidth=2, markersize=6,
            color='red', label='Entropy After Update')
    
    ax1.set_xlabel('Timepoint')
    ax1.set_ylabel('Entropy')
    ax1.set_title(f'{patient_id} - Tree Distribution Entropy Evolution', fontweight='bold')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(timepoint_names, rotation=45)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Entropy changes
    colors = ['green' if change < 0 else 'orange' for change in entropy_changes]
    ax2.bar(x_pos, entropy_changes, color=colors, alpha=0.7)
    ax2.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    ax2.set_xlabel('Timepoint')
    ax2.set_ylabel('Entropy Change')
    ax2.set_title('Entropy Change Per Update (Green=Convergence, Orange=Divergence)', fontweight='bold')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(timepoint_names, rotation=45)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    plot_path = viz_dir / f'{patient_id}_entropy_evolution.png'
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Saved entropy evolution plot: {plot_path}")
    return str(plot_path)


def save_visualization_summary(viz_dir: Path, plot_files: Dict[str, str], 
                             tracking_data: Dict) -> str:
    """
    Save visualization summary JSON file.
    
    Args:
        viz_dir: Visualization directory
        plot_files: Dictionary of created plot files
        tracking_data: Comprehensive tracking data
        
    Returns:
        Path to summary file
    """
    summary_data = {
        'patient_id': tracking_data['patient_id'],
        'analysis_type': 'unified_longitudinal',
        'visualization_timestamp': pd.Timestamp.now().isoformat(),
        'total_timepoints': len(tracking_data['timepoints']),
        'plot_files': plot_files,
        'convergence_summary': tracking_data.get('convergence_summary', {}),
        'final_results_summary': tracking_data.get('final_results', {})
    }
    
    summary_path = viz_dir / 'visualization_summary.json'
    with open(summary_path, 'w') as f:
        json.dump(summary_data, f, indent=2, default=str)
    
    return str(summary_path)