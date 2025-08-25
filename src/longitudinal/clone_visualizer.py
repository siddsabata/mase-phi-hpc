#!/usr/bin/env python3
"""
Clone frequency visualization for longitudinal cancer evolution analysis.

This module provides visualization functions for clone frequency trajectories,
implementing the clone-level plotting system with remainder pseudo-clone support.

Authors: TracerX Pipeline Development Team
"""

import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, Optional, List

logger = logging.getLogger(__name__)


def plot_clone_trajectories(freq_df: pd.DataFrame, output_dir: Path, 
                          analysis_mode: str, patient_id: str,
                          palette: Optional[str] = None) -> Optional[str]:
    """
    Plot clone frequency trajectories over time.
    
    Creates a line plot showing how clone frequencies change over time,
    with the remainder pseudo-clone plotted in semi-transparent grey.
    
    Args:
        freq_df: DataFrame with clone frequencies [sample, time, clone_id, freq]
        output_dir: Output directory for plots
        analysis_mode: Analysis mode ('dynamic' or 'fixed')
        patient_id: Patient identifier
        palette: Color palette name (optional)
        
    Returns:
        Path to saved plot file, or None if failed
    """
    logger.info("Creating clone frequency trajectory plot")
    
    try:
        # Create output directory
        mode_dir = output_dir / f'{analysis_mode}_marker_analysis'
        viz_dir = mode_dir / 'visualizations'
        viz_dir.mkdir(parents=True, exist_ok=True)
        
        # Validate required columns
        required_cols = ['sample', 'time', 'clone_id', 'freq']
        if not all(col in freq_df.columns for col in required_cols):
            logger.error(f"Missing required columns in freq_df: {required_cols}")
            return None
        
        if freq_df.empty:
            logger.warning("Empty frequency DataFrame - cannot create plot")
            return None
        
        # Get unique clones and separate remainder
        all_clones = freq_df['clone_id'].unique()
        regular_clones = [c for c in all_clones if c != 'unrooted_remainder']
        has_remainder = 'unrooted_remainder' in all_clones
        
        logger.info(f"Plotting {len(regular_clones)} regular clones")
        if has_remainder:
            logger.info("Including remainder pseudo-clone")
        
        # Create figure with appropriate size
        plt.figure(figsize=(12, 8))
        
        # Set up color palette
        if palette is None:
            palette = 'Set2'
        
        # Get colors for regular clones
        if regular_clones:
            colors = sns.color_palette(palette, len(regular_clones))
        else:
            colors = []
        
        # Plot regular clones
        for i, clone_id in enumerate(regular_clones):
            clone_data = freq_df[freq_df['clone_id'] == clone_id]
            
            if not clone_data.empty:
                # Sort by time for proper line plotting
                clone_data = clone_data.sort_values('time')
                
                plt.plot(clone_data['time'], clone_data['freq'], 
                        marker='o', linestyle='-', color=colors[i], 
                        markersize=6, linewidth=2, label=f'Clone {clone_id}')
        
        # Plot remainder pseudo-clone (if present) in semi-transparent grey
        if has_remainder:
            remainder_data = freq_df[freq_df['clone_id'] == 'unrooted_remainder']
            if not remainder_data.empty:
                remainder_data = remainder_data.sort_values('time')
                plt.plot(remainder_data['time'], remainder_data['freq'], 
                        marker='s', linestyle='--', color='grey', 
                        markersize=4, linewidth=1.5, alpha=0.7, 
                        label='Remainder')
        
        # Customize plot
        plt.xlabel('Time', fontweight='bold', fontsize=12)
        plt.ylabel('Clone Frequency', fontweight='bold', fontsize=12)
        plt.title(f'{patient_id} - Clone Frequency Trajectories ({analysis_mode.title()} Analysis)',
                 fontsize=14, fontweight='bold')
        
        # Add grid for better readability
        plt.grid(True, alpha=0.3)
        
        # Add legend
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Format time axis if needed
        timepoints = sorted(freq_df['time'].unique())
        if len(timepoints) > 6:
            plt.xticks(rotation=45)
        
        # Set reasonable y-axis limits
        plt.ylim(0, max(freq_df['freq']) * 1.1)
        
        # Save plot
        plot_path = viz_dir / f'{patient_id}_{analysis_mode}_clone_trajectories.png'
        plt.tight_layout()
        plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        logger.info(f"Saved clone trajectory plot: {plot_path}")
        return str(plot_path)
        
    except Exception as e:
        logger.error(f"Error creating clone trajectory plot: {e}")
        return None


def plot_clone_heatmap(freq_df: pd.DataFrame, output_dir: Path, 
                      analysis_mode: str, patient_id: str) -> Optional[str]:
    """
    Create a heatmap of clone frequencies across time.
    
    Args:
        freq_df: DataFrame with clone frequencies
        output_dir: Output directory for plots
        analysis_mode: Analysis mode
        patient_id: Patient identifier
        
    Returns:
        Path to saved heatmap file, or None if failed
    """
    logger.info("Creating clone frequency heatmap")
    
    try:
        # Create output directory
        mode_dir = output_dir / f'{analysis_mode}_marker_analysis'
        viz_dir = mode_dir / 'visualizations'
        viz_dir.mkdir(parents=True, exist_ok=True)
        
        if freq_df.empty:
            logger.warning("Empty frequency DataFrame - cannot create heatmap")
            return None
        
        # Pivot data for heatmap
        # Use first sample if multiple samples exist
        if freq_df['sample'].nunique() > 1:
            first_sample = freq_df['sample'].iloc[0]
            freq_df = freq_df[freq_df['sample'] == first_sample]
            logger.info(f"Using first sample for heatmap: {first_sample}")
        
        pivot_df = freq_df.pivot(index='clone_id', columns='time', values='freq')
        
        # Create heatmap
        plt.figure(figsize=(10, 8))
        
        # Create heatmap with custom colormap
        ax = sns.heatmap(pivot_df, 
                        cmap='YlOrRd',  # Yellow-Orange-Red colormap
                        annot=True,     # Show values
                        fmt='.3f',      # Format numbers
                        cbar_kws={'label': 'Clone Frequency'})
        
        # Customize plot
        plt.title(f'{patient_id} - Clone Frequency Heatmap ({analysis_mode.title()} Analysis)',
                 fontsize=14, fontweight='bold')
        plt.xlabel('Time', fontweight='bold', fontsize=12)
        plt.ylabel('Clone ID', fontweight='bold', fontsize=12)
        
        # Rotate time labels if needed
        if len(pivot_df.columns) > 6:
            plt.xticks(rotation=45)
        
        # Save heatmap
        heatmap_path = viz_dir / f'{patient_id}_{analysis_mode}_clone_heatmap.png'
        plt.tight_layout()
        plt.savefig(heatmap_path, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        logger.info(f"Saved clone frequency heatmap: {heatmap_path}")
        return str(heatmap_path)
        
    except Exception as e:
        logger.error(f"Error creating clone frequency heatmap: {e}")
        return None


def create_clone_frequency_summary(freq_df: pd.DataFrame, output_dir: Path,
                                 analysis_mode: str, patient_id: str) -> Optional[str]:
    """
    Create summary statistics for clone frequencies.
    
    Args:
        freq_df: DataFrame with clone frequencies
        output_dir: Output directory for results
        analysis_mode: Analysis mode
        patient_id: Patient identifier
        
    Returns:
        Path to saved summary file, or None if failed
    """
    logger.info("Creating clone frequency summary statistics")
    
    try:
        # Create output directory
        mode_dir = output_dir / f'{analysis_mode}_marker_analysis'
        mode_dir.mkdir(parents=True, exist_ok=True)
        
        if freq_df.empty:
            logger.warning("Empty frequency DataFrame - cannot create summary")
            return None
        
        # Calculate summary statistics
        summary_stats = {}
        
        # Overall statistics
        summary_stats['overall'] = {
            'total_clones': freq_df['clone_id'].nunique(),
            'total_timepoints': freq_df['time'].nunique(),
            'total_samples': freq_df['sample'].nunique(),
            'mean_frequency': freq_df['freq'].mean(),
            'median_frequency': freq_df['freq'].median(),
            'max_frequency': freq_df['freq'].max(),
            'min_frequency': freq_df['freq'].min()
        }
        
        # Per-clone statistics
        clone_stats = {}
        for clone_id in freq_df['clone_id'].unique():
            clone_data = freq_df[freq_df['clone_id'] == clone_id]
            clone_stats[str(clone_id)] = {
                'mean_freq': clone_data['freq'].mean(),
                'median_freq': clone_data['freq'].median(),
                'max_freq': clone_data['freq'].max(),
                'min_freq': clone_data['freq'].min(),
                'std_freq': clone_data['freq'].std(),
                'n_observations': len(clone_data)
            }
        
        summary_stats['per_clone'] = clone_stats
        
        # Per-timepoint statistics
        timepoint_stats = {}
        for timepoint in freq_df['time'].unique():
            timepoint_data = freq_df[freq_df['time'] == timepoint]
            timepoint_stats[str(timepoint)] = {
                'n_clones': timepoint_data['clone_id'].nunique(),
                'mean_freq': timepoint_data['freq'].mean(),
                'total_freq': timepoint_data['freq'].sum(),
                'dominant_clone': timepoint_data.loc[timepoint_data['freq'].idxmax(), 'clone_id']
            }
        
        summary_stats['per_timepoint'] = timepoint_stats
        
        # Save summary
        summary_file = mode_dir / 'clone_frequency_summary.json'
        import json
        with open(summary_file, 'w') as f:
            json.dump(summary_stats, f, indent=2, default=str)
        
        logger.info(f"Saved clone frequency summary: {summary_file}")
        return str(summary_file)
        
    except Exception as e:
        logger.error(f"Error creating clone frequency summary: {e}")
        return None


def validate_frequency_data_for_plotting(freq_df: pd.DataFrame) -> bool:
    """
    Validate frequency data is suitable for plotting.
    
    Args:
        freq_df: DataFrame with clone frequencies
        
    Returns:
        True if data is valid for plotting, False otherwise
    """
    required_cols = ['sample', 'time', 'clone_id', 'freq']
    
    # Check required columns
    if not all(col in freq_df.columns for col in required_cols):
        logger.error(f"Missing required columns: {required_cols}")
        return False
    
    # Check for empty DataFrame
    if freq_df.empty:
        logger.error("Empty frequency DataFrame")
        return False
    
    # Check for null values
    if freq_df[required_cols].isnull().any().any():
        logger.error("Found null values in frequency data")
        return False
    
    # Check for reasonable frequency values
    if (freq_df['freq'] < 0).any():
        logger.error("Found negative frequencies")
        return False
    
    logger.info("Frequency data validation passed for plotting")
    return True
