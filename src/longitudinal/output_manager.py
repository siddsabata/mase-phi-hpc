#!/usr/bin/env python3
"""
Output management for longitudinal cancer evolution analysis.

This module provides standardized output directory structure and 
result saving functionality for both dynamic and fixed analysis modes.

Authors: TracerX Pipeline Development Team
"""

import logging
import json
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Dict, Any

logger = logging.getLogger(__name__)


def setup_logging(output_dir: Path, patient_id: str) -> logging.Logger:
    """
    Set up comprehensive logging for the longitudinal analysis pipeline.
    
    Args:
        output_dir: Directory where log files will be stored
        patient_id: Patient identifier for log file naming
        
    Returns:
        Configured logger instance
    """
    # Create logs directory
    log_dir = output_dir / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Create timestamp for log file
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = log_dir / f'{patient_id}_longitudinal_analysis_{timestamp}.log'
    
    # Configure logging format
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    
    logger = logging.getLogger('longitudinal_analysis')
    logger.info(f"Starting longitudinal analysis for patient {patient_id}")
    logger.info(f"Log file: {log_file}")
    
    return logger


def create_standard_output_structure(output_dir: Path, analysis_mode: str) -> Dict[str, Path]:
    """
    Create standardized output directory structure for analysis results.
    
    Args:
        output_dir: Base output directory
        analysis_mode: Analysis mode ('dynamic' or 'fixed')
        
    Returns:
        Dictionary mapping directory types to paths
    """
    # Create main output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create mode-specific subdirectories
    mode_dir = output_dir / f'{analysis_mode}_marker_analysis'
    
    directories = {
        'base': output_dir,
        'mode_base': mode_dir,
        'updated_trees': mode_dir / 'updated_trees',
        'marker_data': mode_dir / 'marker_data',
        'logs': output_dir / 'logs'
    }
    
    # Create all directories
    for dir_path in directories.values():
        dir_path.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Created output directory structure for {analysis_mode} analysis")
    logger.info(f"Base directory: {output_dir}")
    logger.info(f"Mode directory: {mode_dir}")
    
    return directories


def save_analysis_summary(output_dir: Path, summary_data: Dict[str, Any], 
                         analysis_mode: str) -> Path:
    """
    Save analysis summary in standardized format.
    
    Args:
        output_dir: Output directory
        summary_data: Analysis summary data
        analysis_mode: Analysis mode ('dynamic' or 'fixed')
        
    Returns:
        Path to saved summary file
    """
    summary_file = output_dir / 'analysis_summary.json'
    
    # Add standard metadata
    summary_with_metadata = {
        'analysis_metadata': {
            'analysis_mode': analysis_mode,
            'timestamp': datetime.now().isoformat(),
            'pipeline_version': '2.0.0',
            'format_version': '1.0'
        },
        **summary_data
    }
    
    with open(summary_file, 'w') as f:
        json.dump(summary_with_metadata, f, indent=2, default=str)
    
    logger.info(f"Saved analysis summary: {summary_file}")
    return summary_file


def save_final_report(output_dir: Path, args, results_summary: Dict, 
                     timepoint_data: Dict) -> Path:
    """
    Generate and save a comprehensive final report.
    
    Args:
        output_dir: Output directory
        args: Configuration arguments
        results_summary: Analysis results summary
        timepoint_data: Timepoint data dictionary
        
    Returns:
        Path to saved report file
    """
    report_data = {
        'patient_id': args.patient_id,
        'analysis_timestamp': datetime.now().isoformat(),
        'analysis_mode': args.analysis_mode,
        'configuration': {
            'n_markers': args.n_markers,
            'read_depth': args.read_depth,
            'method': args.method,
            'lambda1': args.lambda1,
            'lambda2': args.lambda2
        },
        'input_files': {
            'aggregation_dir': str(args.aggregation_dir),
            'ssm_file': str(args.ssm_file),
            'longitudinal_data': str(args.longitudinal_data)
        },
        'analysis_results': {
            'total_timepoints': len(timepoint_data),
            'timepoints_processed': sorted(timepoint_data.keys()),
            'output_directory': str(output_dir),
            'analysis_mode': args.analysis_mode
        },
        'summary': results_summary.get('analysis_summary', {})
    }
    
    report_file = output_dir / 'final_report.json'
    with open(report_file, 'w') as f:
        json.dump(report_data, f, indent=2, default=str)
    
    logger.info(f"Generated final report: {report_file}")
    logger.info(f"Analysis completed for patient {args.patient_id}")
    
    return report_file


def write_clone_frequencies(freq_df: pd.DataFrame, output_dir: Path, 
                          analysis_mode: str) -> Path:
    """
    Write clone frequency data to CSV file.
    
    This function saves the clone frequency DataFrame in a standardized format
    alongside other analysis results.
    
    Args:
        freq_df: DataFrame with clone frequencies [sample, time, clone_id, freq]
        output_dir: Output directory for results
        analysis_mode: Analysis mode ('dynamic' or 'fixed')
        
    Returns:
        Path to saved clone frequencies file
    """
    logger.info("Writing clone frequencies to CSV file")
    
    # Create mode-specific output directory
    mode_dir = output_dir / f'{analysis_mode}_marker_analysis'
    mode_dir.mkdir(parents=True, exist_ok=True)
    
    # Define output file path
    clone_freq_file = mode_dir / 'clone_frequencies.csv'
    
    # Sort data for consistent output
    freq_df_sorted = freq_df.sort_values(['sample', 'time', 'clone_id'])
    
    # Save to CSV with standard formatting
    freq_df_sorted.to_csv(clone_freq_file, index=False, float_format='%.6f')
    
    logger.info(f"Saved clone frequencies: {clone_freq_file}")
    logger.info(f"Clone frequency data shape: {freq_df_sorted.shape}")
    logger.info(f"Samples: {freq_df_sorted['sample'].nunique()}")
    logger.info(f"Timepoints: {freq_df_sorted['time'].nunique()}")
    logger.info(f"Clones: {freq_df_sorted['clone_id'].nunique()}")
    
    return clone_freq_file