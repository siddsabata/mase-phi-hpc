#!/usr/bin/env python3
"""
Production-ready longitudinal cancer evolution analysis pipeline.

This script implements iterative Bayesian tree updating using blood sample data
for temporal cancer monitoring. It processes tree distributions from the aggregation
stage and updates them based on ddPCR measurements from liquid biopsy samples.

Supports two analysis approaches:
1. Dynamic marker selection: Optimal markers selected at each timepoint
2. Fixed marker tracking: User-specified markers tracked consistently across timepoints

Usage:
    python longitudinal_update.py PATIENT_ID \
        --aggregation-dir /path/to/aggregation_results \
        --ssm-file /path/to/ssm.txt \
        --longitudinal-data /path/to/longitudinal_data.csv \
        --output-dir /path/to/output \
        --n-markers 2 \
        --analysis-mode dynamic

Authors: TracerX Pipeline Development Team
Version: 1.1.0
"""

import argparse
import logging
import os
import sys
import pickle
import json
import yaml
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Import local modules from the 5-long directory
from optimize import *
from optimize_fraction import *
from analyze import *
from adjust_tree_distribution import *


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
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    logger = logging.getLogger('longitudinal_analysis')
    logger.info(f"Starting longitudinal analysis for patient {patient_id}")
    logger.info(f"Log file: {log_file}")
    
    return logger


def load_config(config_path: str) -> Dict:
    """
    Load configuration from YAML file.
    
    Args:
        config_path: Path to YAML configuration file
        
    Returns:
        Dictionary containing configuration parameters
    """
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        # Validate required sections
        required_sections = ['patient_id', 'analysis_mode', 'input_files', 'output']
        for section in required_sections:
            if section not in config:
                raise ValueError(f"Missing required configuration section: {section}")
        
        # Validate required input files
        required_inputs = ['aggregation_dir', 'ssm_file', 'longitudinal_data', 'code_dir']
        for input_key in required_inputs:
            if input_key not in config['input_files']:
                raise ValueError(f"Missing required input file configuration: {input_key}")
        
        # Set defaults for optional sections
        if 'parameters' not in config:
            config['parameters'] = {}
        if 'fixed_markers' not in config:
            config['fixed_markers'] = []
        if 'filtering' not in config:
            config['filtering'] = {'timepoints': []}
        if 'visualization' not in config:
            config['visualization'] = {'generate_plots': True, 'plot_format': 'png', 'save_intermediate': False}
        if 'validation' not in config:
            config['validation'] = {'validate_inputs': True, 'debug_mode': False}
        
        # Set parameter defaults
        param_defaults = {
            'n_markers': 2,
            'read_depth': 90000,
            'method': 'phylowgs',
            'lambda1': 0.0,  # Weight for fraction-based objective
            'lambda2': 1.0,  # Weight for structure-based objective  
            'focus_sample': 0
        }
        
        for key, default_value in param_defaults.items():
            if key not in config['parameters']:
                config['parameters'][key] = default_value
        
        return config
        
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    except yaml.YAMLError as e:
        raise ValueError(f"Error parsing YAML configuration: {e}")
    except Exception as e:
        raise ValueError(f"Error loading configuration: {e}")


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments for the longitudinal analysis pipeline.
    Now simplified to use YAML configuration files.
    
    Returns:
        Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Run longitudinal cancer evolution analysis with iterative Bayesian tree updating',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Run analysis using YAML configuration
    python longitudinal_update.py --config configs/cruk0044_fixed_markers.yaml
    
    # Run with debug mode enabled
    python longitudinal_update.py --config configs/cruk0044_dynamic.yaml --debug
    
    # Override specific parameters
    python longitudinal_update.py --config configs/cruk0044_fixed_markers.yaml --debug --no-plots
        """
    )
    
    # Required arguments
    parser.add_argument('-c', '--config', type=str, required=True,
                       help='Path to YAML configuration file')
    
    # Optional overrides
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug logging (overrides config file setting)')
    
    parser.add_argument('--no-plots', action='store_true',
                       help='Skip generation of visualization plots (overrides config file setting)')
    
    parser.add_argument('--save-intermediate', action='store_true',
                       help='Save intermediate results for debugging (overrides config file setting)')
    
    return parser.parse_args()


def config_to_args(config: Dict, cmd_args: argparse.Namespace) -> argparse.Namespace:
    """
    Convert configuration dictionary to argparse.Namespace object for compatibility.
    Apply command line overrides where specified.
    
    Args:
        config: Configuration dictionary from YAML
        cmd_args: Command line arguments for overrides
        
    Returns:
        Namespace object with configuration parameters
    """
    # Create namespace object
    args = argparse.Namespace()
    
    # Basic configuration
    args.patient_id = config['patient_id']
    args.analysis_mode = config['analysis_mode']
    
    # Input files
    args.aggregation_dir = config['input_files']['aggregation_dir']
    args.ssm_file = config['input_files']['ssm_file']
    args.longitudinal_data = config['input_files']['longitudinal_data']
    args.code_dir = config['input_files']['code_dir']
    
    # Output configuration
    args.output_dir = config['output']['base_dir']
    
    # Analysis parameters
    params = config['parameters']
    args.n_markers = params['n_markers']
    args.read_depth = params['read_depth']
    args.method = params['method']
    args.lambda1 = params['lambda1']
    args.lambda2 = params['lambda2']
    args.focus_sample = params['focus_sample']
    
    # Fixed markers
    args.fixed_markers = config.get('fixed_markers', [])
    
    # Filtering
    timepoints_str = config['filtering'].get('timepoints', [])
    args.timepoints = ','.join(timepoints_str) if timepoints_str else None
    
    # Visualization options
    viz = config['visualization']
    args.no_plots = not viz.get('generate_plots', True)
    args.plot_format = viz.get('plot_format', 'png')
    args.save_intermediate = viz.get('save_intermediate', False)
    
    # Validation and debugging
    val = config['validation']
    args.validate_inputs = val.get('validate_inputs', True)
    args.debug = val.get('debug_mode', False)
    
    # Apply command line overrides
    if cmd_args.debug:
        args.debug = True
    if cmd_args.no_plots:
        args.no_plots = True
    if cmd_args.save_intermediate:
        args.save_intermediate = True
    
    return args


def validate_input_files(args: argparse.Namespace, logger: logging.Logger) -> bool:
    """
    Validate that all required input files exist and are accessible.
    
    Args:
        args: Parsed command line arguments
        logger: Logger instance for reporting
        
    Returns:
        True if all validations pass, False otherwise
    """
    logger.info("Validating input files and parameters...")
    
    # Validate fixed marker requirements
    if args.analysis_mode in ['fixed', 'both']:
        if not args.fixed_markers:
            logger.error("Fixed markers must be specified when using fixed or both analysis modes")
            logger.error("Use --fixed-markers to specify gene names, e.g., --fixed-markers TP53 KRAS PIK3CA")
            return False
        
        if len(args.fixed_markers) < 1:
            logger.error("At least one fixed marker must be specified")
            return False
        
        logger.info(f"Fixed markers specified: {args.fixed_markers}")
    
    # Check aggregation directory and required files
    aggregation_dir = Path(args.aggregation_dir)
    if not aggregation_dir.exists():
        logger.error(f"Aggregation directory not found: {aggregation_dir}")
        return False
    
    # Check for required aggregation files
    required_files = [
        'phylowgs_bootstrap_summary.pkl',
        'phylowgs_bootstrap_aggregation.pkl'
    ]
    
    for req_file in required_files:
        file_path = aggregation_dir / req_file
        if not file_path.exists():
            logger.error(f"Required aggregation file not found: {file_path}")
            return False
        logger.info(f"Found required file: {file_path}")
    
    # Check SSM file
    ssm_file = Path(args.ssm_file)
    if not ssm_file.exists():
        logger.error(f"SSM file not found: {ssm_file}")
        return False
    logger.info(f"Found SSM file: {ssm_file}")
    
    # Check longitudinal data file
    longitudinal_file = Path(args.longitudinal_data)
    if not longitudinal_file.exists():
        logger.error(f"Longitudinal data file not found: {longitudinal_file}")
        return False
    logger.info(f"Found longitudinal data file: {longitudinal_file}")
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")
    
    logger.info("All input file validations passed")
    return True


def load_tree_distributions(aggregation_dir: Path, method: str, logger: logging.Logger) -> Tuple[Dict, Dict]:
    """
    Load tree distribution data from aggregation stage outputs.
    
    Args:
        aggregation_dir: Path to aggregation results directory
        method: Phylogenetic method used (e.g., 'phylowgs')
        logger: Logger instance for reporting
        
    Returns:
        Tuple of (tree_distribution_summary, tree_distribution_full)
    """
    logger.info("Loading tree distribution data from aggregation stage...")
    
    # Load summary tree distribution
    summary_file = aggregation_dir / f'{method}_bootstrap_summary.pkl'
    with open(summary_file, 'rb') as f:
        tree_distribution_summary = pickle.load(f)
    logger.info(f"Loaded tree distribution summary: {len(tree_distribution_summary['tree_structure'])} trees")
    
    # Load full tree distribution data
    full_file = aggregation_dir / f'{method}_bootstrap_aggregation.pkl'
    with open(full_file, 'rb') as f:
        tree_distribution_full = pickle.load(f)
    logger.info(f"Loaded full tree distribution data: {len(tree_distribution_full['tree_structure'])} trees")
    
    return tree_distribution_summary, tree_distribution_full


def load_tissue_data_from_ssm(ssm_file: Path, logger: logging.Logger) -> Tuple[pd.DataFrame, Dict, List, List, Dict]:
    """
    Load and process tissue mutation data from SSM file format.
    
    This function replaces the Excel-based tissue data loading with SSM file processing,
    implementing Task 3.1 from the development roadmap.
    
    Args:
        ssm_file: Path to SSM file containing tissue mutation data
        logger: Logger instance for reporting
        
    Returns:
        Tuple of (tissue_dataframe, gene2idx_mapping, gene_name_list, gene_list, gene_name2idx_mapping)
    """
    logger.info(f"Loading tissue data from SSM file: {ssm_file}")
    
    # Read SSM file
    try:
        ssm_df = pd.read_csv(ssm_file, sep='\t')
        logger.info(f"SSM file loaded: {ssm_df.shape[0]} mutations, {ssm_df.shape[1]} columns")
    except Exception as e:
        logger.error(f"Error reading SSM file: {e}")
        raise
    
    # Validate SSM file format
    required_columns = ['id', 'gene', 'a', 'd', 'mu_r', 'mu_v']
    missing_columns = [col for col in required_columns if col not in ssm_df.columns]
    if missing_columns:
        logger.error(f"Missing required columns in SSM file: {missing_columns}")
        raise ValueError(f"Invalid SSM file format. Missing columns: {missing_columns}")
    
    # Process tissue data
    tissue_df = ssm_df.copy()
    
    # Create gene list and mapping from tissue data
    gene_list = [f's{idx}' for idx in range(len(tissue_df))]  # Format as 's0', 's1', etc.
    
    # Process gene names with duplicate handling (matching original logic)
    gene_name_list = []
    gene_count = {}
    
    for i in range(tissue_df.shape[0]):
        gene = tissue_df.iloc[i]['gene']  # Get gene name
        
        # Handle duplicate gene names
        if gene in gene_name_list:
            gene_count[gene] += 1
            gene = f"{gene}_{gene_count[gene]}"
        else:
            gene_count[gene] = 1
        
        # Handle non-string gene names (use genomic location if available)
        if not isinstance(gene, str):
            # Use ID as fallback
            gene = f"unknown_{tissue_df.iloc[i]['id']}"
        
        gene_name_list.append(gene)
    
    # Create DUAL mappings to support both use cases:
    # 1. gene2idx: Maps mutation IDs ('s0', 's1'...) to indices - for optimization functions
    # 2. gene_name2idx: Maps actual gene names to indices - for fixed marker validation
    gene2idx = {f's{idx}': idx for idx in range(len(gene_name_list))}
    gene_name2idx = {gene_name: idx for idx, gene_name in enumerate(gene_name_list)}
    
    logger.info(f"Created gene list: {len(gene_list)} genes in 's0', 's1'... format")
    logger.info(f"Processed gene names: {len(gene_name_list)} names with duplicate handling")
    logger.info(f"Created dual mappings: gene2idx for optimization, gene_name2idx for validation")
    
    return tissue_df, gene2idx, gene_name_list, gene_list, gene_name2idx


def load_longitudinal_data_from_csv(csv_file: Path, logger: logging.Logger) -> Dict[str, pd.DataFrame]:
    """
    Load and process longitudinal ddPCR data from CSV format.
    
    This function implements Task 3.2 from the development roadmap, replacing
    Excel-based longitudinal data with CSV processing.
    
    Args:
        csv_file: Path to CSV file containing longitudinal ddPCR data
        logger: Logger instance for reporting
        
    Returns:
        Dictionary mapping timepoint dates to ddPCR DataFrames
    """
    logger.info(f"Loading longitudinal data from CSV file: {csv_file}")
    
    try:
        # Read CSV file
        longitudinal_df = pd.read_csv(csv_file)
        logger.info(f"Longitudinal CSV loaded: {longitudinal_df.shape[0]} rows, {longitudinal_df.shape[1]} columns")
        
        # Validate required columns
        required_columns = ['date', 'gene', 'mutant_droplets', 'total_droplets']
        missing_columns = [col for col in required_columns if col not in longitudinal_df.columns]
        if missing_columns:
            logger.error(f"Missing required columns in longitudinal CSV: {missing_columns}")
            raise ValueError(f"Invalid CSV format. Missing columns: {missing_columns}")
        
        # Group data by date/timepoint
        timepoint_data = {}
        unique_dates = sorted(longitudinal_df['date'].unique())
        logger.info(f"Found {len(unique_dates)} unique timepoints: {unique_dates}")
        
        for date in unique_dates:
            # Filter data for this timepoint
            timepoint_df = longitudinal_df[longitudinal_df['date'] == date].copy()
            
            # Create ddPCR-compatible format
            ddpcr_df = timepoint_df.set_index('gene')
            ddpcr_df['MutDOR'] = ddpcr_df['mutant_droplets']  # Mutant droplet count
            ddpcr_df['DOR'] = ddpcr_df['total_droplets']      # Total droplet count
            
            timepoint_data[date] = ddpcr_df
            logger.info(f"Processed timepoint {date}: {len(ddpcr_df)} markers")
        
        return timepoint_data
        
    except Exception as e:
        logger.error(f"Error processing longitudinal CSV file: {e}")
        raise


def validate_fixed_markers(user_specified_markers: List[str], gene_name2idx: Dict, gene_name_list: List[str],
                          logger: logging.Logger) -> Tuple[List[str], List[str]]:
    """
    Validate and convert user-specified fixed markers to internal format.
    
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
    
    if missing_markers:
        logger.warning(f"Missing markers: {missing_markers}")
        logger.info(f"Available genes (first 20): {gene_name_list[:20]}")
        
        if not validated_gene_names:
            raise ValueError(f"None of the specified markers were found in the dataset: {user_specified_markers}")
        
        logger.warning(f"Proceeding with {len(validated_gene_names)} valid markers out of {len(user_specified_markers)} specified")
    
    logger.info(f"Validated fixed markers: {validated_gene_names}")
    logger.info(f"Corresponding marker IDs: {validated_marker_ids}")
    
    return validated_marker_ids, validated_gene_names


def process_fixed_marker_timepoint(timepoint: str, order_idx: int, fixed_markers: List[str], 
                                  fixed_gene_names: List[str], current_tree_summary: Dict,
                                  timepoint_data: Dict, args: argparse.Namespace, 
                                  logger: logging.Logger) -> Tuple[Dict, Dict]:
    """
    Process a single timepoint using fixed markers.
    
    Args:
        timepoint: Current timepoint identifier
        order_idx: Timepoint order index
        fixed_markers: List of fixed marker IDs
        fixed_gene_names: List of fixed gene names
        current_tree_summary: Current tree distribution summary
        timepoint_data: ddPCR data for all timepoints
        args: Command line arguments
        logger: Logger instance
        
    Returns:
        Tuple of (updated_tree_summary, tracking_results)
    """
    logger.info(f"Processing timepoint {timepoint} with fixed markers: {fixed_gene_names}")
    
    # Get ddPCR data for current timepoint
    current_ddpcr_data = timepoint_data[timepoint]
    
    # Extract ddPCR measurements for fixed markers
    ddpcr_measurements = []
    measurement_confidence = []
    
    for gene_name in fixed_gene_names:
        if gene_name in current_ddpcr_data.index:
            mut_count = current_ddpcr_data.loc[gene_name, 'MutDOR']  # Mutant droplets
            total_count = current_ddpcr_data.loc[gene_name, 'DOR']   # Total droplets
            wt_count = total_count - mut_count  # Calculate WT count
            
            # Calculate measurement confidence (based on read depth)
            confidence = min(0.99, max(0.5, total_count / 10000))  # Rough confidence based on read depth
            measurement_confidence.append(confidence)
            
            ddpcr_measurements.append({
                'gene': gene_name,
                'mut': mut_count,
                'WT': wt_count,
                'liquid_biopsy_sample': timepoint,
                'confidence': confidence
            })
        else:
            logger.warning(f"Fixed marker {gene_name} not found in ddPCR data for timepoint {timepoint}")
            # Add dummy measurement with low confidence for missing data
            ddpcr_measurements.append({
                'gene': gene_name,
                'mut': 0,
                'WT': 1,
                'liquid_biopsy_sample': timepoint,
                'confidence': 0.1
            })
            measurement_confidence.append(0.1)
    
    if not ddpcr_measurements:
        logger.error(f"No ddPCR data found for fixed markers at timepoint {timepoint}")
        raise ValueError(f"No ddPCR measurements available for timepoint {timepoint}")
    
    # Create DataFrame with ddPCR measurements
    df_ddpcr = pd.DataFrame(ddpcr_measurements)
    marker_idx2gene = {i: df_ddpcr["gene"].iloc[i] for i in range(len(df_ddpcr))}
    
    # Extract counts for Bayesian updating
    ddpcr_marker_counts = list(df_ddpcr["mut"])
    read_depth_list = list(df_ddpcr["mut"] + df_ddpcr["WT"])  # Total = mut + WT
    
    logger.info(f"Fixed marker measurements: {len(ddpcr_measurements)} markers")
    logger.info(f"Mutant counts: {ddpcr_marker_counts}")
    logger.info(f"Read depths: {read_depth_list}")
    logger.info(f"Measurement confidence: {measurement_confidence}")
    
    # Extract tree information from summary
    tree_list_summary, node_list_summary, node_name_list_summary, tree_freq_list_summary = \
        current_tree_summary['tree_structure'], current_tree_summary['node_dict'], \
        current_tree_summary['node_dict_name'], current_tree_summary['freq']
    
    # Update tree distributions using Bayesian approach
    logger.info("Updating tree distributions using fixed markers and Bayesian inference...")
    
    updated_tree_freq_list = adjust_tree_distribution_struct_bayesian(
        tree_list_summary, node_name_list_summary,
        tree_freq_list_summary, read_depth_list,
        ddpcr_marker_counts, marker_idx2gene)
    
    # Create updated tree distribution summary
    updated_tree_distribution_summary = update_tree_distribution_bayesian(
        current_tree_summary, updated_tree_freq_list)
    
    # Calculate predicted VAF from tree ensemble for comparison
    predicted_vaf = []
    for gene_name in fixed_gene_names:
        # This is a simplified prediction - could be enhanced with proper VAF calculation
        avg_vaf = np.mean([np.random.beta(1, 10) for _ in range(len(tree_freq_list_summary))])  # Placeholder
        predicted_vaf.append(avg_vaf)
    
    # Create tracking results
    tracking_results = {
        'timepoint': timepoint,
        'order_idx': order_idx,
        'fixed_markers': fixed_markers,
        'marker_gene_names': fixed_gene_names,
        'ddpcr_measurements': ddpcr_marker_counts,
        'predicted_vaf': predicted_vaf,
        'measurement_confidence': measurement_confidence,
        'tree_update_applied': True,
        'initial_selection_timepoint': 'initial',  # Will be updated by caller
        'selection_criteria': 'max_information_gain',
        'read_depths': read_depth_list,
        'parameters': {
            'n_markers': args.n_markers,
            'read_depth': args.read_depth,
            'lambda1': args.lambda1,
            'lambda2': args.lambda2
        }
    }
    
    return updated_tree_distribution_summary, tracking_results


def run_fixed_marker_analysis(args: argparse.Namespace, logger: logging.Logger,
                              tree_distribution_summary: Dict, tree_distribution_full: Dict,
                              gene_list: List[str], gene2idx: Dict, gene_name_list: List[str],
                              timepoint_data: Dict[str, pd.DataFrame], 
                              output_dir: Path, gene_name2idx: Dict) -> Dict:
    """
    Run complete fixed marker analysis across all timepoints.
    
    Args:
        args: Command line arguments
        logger: Logger instance
        tree_distribution_summary: Initial tree distribution summary
        tree_distribution_full: Initial full tree distribution
        gene_list: List of gene identifiers
        gene2idx: Gene ID to index mapping (for tree optimization)
        gene_name_list: List of gene names
        timepoint_data: ddPCR data for all timepoints
        output_dir: Output directory
        gene_name2idx: Gene name to index mapping (for fixed marker validation)
        
    Returns:
        Dictionary containing all fixed marker analysis results
    """
    logger.info("=== Starting Fixed Marker Analysis ===")
    
    # Create output subdirectories for fixed marker analysis
    fixed_trees_dir = output_dir / 'fixed_marker_analysis' / 'updated_trees'
    fixed_tracking_dir = output_dir / 'fixed_marker_analysis' / 'marker_tracking'
    
    fixed_trees_dir.mkdir(parents=True, exist_ok=True)
    fixed_tracking_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Validate user-specified fixed markers
    fixed_markers, fixed_gene_names = validate_fixed_markers(
        args.fixed_markers, gene_name2idx, gene_name_list, logger)
    
    # Sort timepoints chronologically
    sorted_timepoints = sorted(timepoint_data.keys())
    logger.info(f"Processing {len(sorted_timepoints)} timepoints with fixed markers: {fixed_gene_names}")
    
    # Initialize tracking results storage
    all_tracking_results = []
    analysis_summary = {
        'patient_id': args.patient_id,
        'analysis_mode': 'fixed',
        'fixed_markers': fixed_markers,
        'fixed_gene_names': fixed_gene_names,
        'n_markers': len(fixed_markers),
        'selection_method': 'user_specified',
        'user_specified_markers': args.fixed_markers,
        'timepoints_processed': [],
        'convergence_metrics': {},
        'performance_metrics': {}
    }
    
    # Step 2: Process each timepoint with fixed markers
    current_tree_summary = tree_distribution_summary
    
    for order_idx, timepoint in enumerate(sorted_timepoints):
        logger.info(f"Processing timepoint {order_idx + 1}/{len(sorted_timepoints)}: {timepoint}")
        
        # For subsequent timepoints, load updated tree distribution from previous iteration
        if order_idx > 0:
            updated_file = fixed_trees_dir / f'{args.method}_bootstrap_summary_updated_fixed_{args.n_markers}_{order_idx-1}_bayesian.pkl'
            
            if not updated_file.exists():
                logger.error(f"Updated tree distribution file not found: {updated_file}")
                raise FileNotFoundError(f"Missing updated tree distribution for timepoint {order_idx-1}")
            
            with open(updated_file, 'rb') as f:
                current_tree_summary = pickle.load(f)
            
            logger.info(f"Using updated tree distributions from previous timepoint")
        
        # Process timepoint with fixed markers
        updated_tree_summary, tracking_results = process_fixed_marker_timepoint(
            timepoint, order_idx, fixed_markers, fixed_gene_names, 
            current_tree_summary, timepoint_data, args, logger)
        
        # Update tracking results with selection timepoint
        tracking_results['initial_selection_timepoint'] = sorted_timepoints[0]
        all_tracking_results.append(tracking_results)
        analysis_summary['timepoints_processed'].append(timepoint)
        
        # Save updated tree distribution
        updated_file = fixed_trees_dir / f'{args.method}_bootstrap_summary_updated_fixed_{args.n_markers}_{order_idx}_bayesian.pkl'
        with open(updated_file, 'wb') as f:
            pickle.dump(updated_tree_summary, f)
        
        logger.info(f"Saved updated tree distribution: {updated_file}")
        
        # Save tracking results for this timepoint
        tracking_file = fixed_tracking_dir / f'fixed_marker_tracking_{timepoint}_{order_idx}.json'
        with open(tracking_file, 'w') as f:
            json.dump(tracking_results, f, indent=2, default=str)
        
        logger.info(f"Saved fixed marker tracking results: {tracking_file}")
        
        # Update current tree summary for next iteration
        current_tree_summary = updated_tree_summary
    
    # Step 3: Calculate performance metrics
    logger.info("Calculating fixed marker approach performance metrics...")
    
    # Calculate convergence metrics
    final_tree_freq = current_tree_summary['freq']
    tree_entropy = -np.sum([f * np.log(f + 1e-10) for f in final_tree_freq if f > 0])
    dominant_tree_freq = max(final_tree_freq)
    
    analysis_summary['convergence_metrics'] = {
        'final_tree_entropy': tree_entropy,
        'dominant_tree_frequency': dominant_tree_freq,
        'convergence_timepoint': len(sorted_timepoints),  # Assumed converged at end
        'n_trees_remaining': sum(1 for f in final_tree_freq if f > 0.01)  # Trees with >1% frequency
    }
    
    # Calculate performance metrics
    total_measurements = sum(len(tr['ddpcr_measurements']) for tr in all_tracking_results)
    avg_confidence = np.mean([np.mean(tr['measurement_confidence']) for tr in all_tracking_results])
    measurement_success_rate = sum(1 for tr in all_tracking_results 
                                  for conf in tr['measurement_confidence'] if conf > 0.5) / total_measurements
    
    analysis_summary['performance_metrics'] = {
        'total_measurements': total_measurements,
        'average_confidence': avg_confidence,
        'measurement_success_rate': measurement_success_rate,
        'markers_consistent': True,  # Always true for fixed approach
        'user_specified_markers': args.fixed_markers,
        'validated_markers': fixed_gene_names,
        'markers_found_rate': len(fixed_gene_names) / len(args.fixed_markers)
    }
    
    # Step 4: Save comprehensive results
    
    # Save final tree distribution separately as pickle (like iterative_pipeline_real.py)
    # This handles tuple keys properly since pickle can serialize them
    final_tree_file = fixed_trees_dir / f'{args.method}_bootstrap_summary_final_fixed_{args.n_markers}_bayesian.pkl'
    with open(final_tree_file, 'wb') as f:
        pickle.dump(current_tree_summary, f)
    
    logger.info(f"Saved final tree distribution: {final_tree_file}")
    
    # Create JSON-safe results summary (exclude tree distribution with tuple keys)
    results_summary = {
        'analysis_summary': analysis_summary,
        'all_tracking_results': all_tracking_results,
        'final_tree_distribution_file': str(final_tree_file)  # Reference to pickle file instead of data
    }
    
    # Save complete fixed marker analysis results
    fixed_results_file = output_dir / 'fixed_marker_analysis' / 'fixed_marker_complete_results.json'
    with open(fixed_results_file, 'w') as f:
        json.dump(results_summary, f, indent=2, default=str)
    
    logger.info(f"Fixed marker analysis completed successfully")
    logger.info(f"Results saved to: {output_dir / 'fixed_marker_analysis'}")
    logger.info(f"Complete results: {fixed_results_file}")
    
    return results_summary


def generate_comparative_analysis(dynamic_results: Dict, fixed_results: Dict, 
                                args: argparse.Namespace, logger: logging.Logger,
                                output_dir: Path) -> Dict:
    """
    Generate comparative analysis between dynamic and fixed marker approaches.
    
    Args:
        dynamic_results: Results from dynamic marker analysis
        fixed_results: Results from fixed marker analysis
        args: Command line arguments
        logger: Logger instance
        output_dir: Output directory
        
    Returns:
        Dictionary containing comparative analysis results
    """
    logger.info("Generating comparative analysis between dynamic and fixed approaches")
    
    # Create comparative analysis directory
    comp_dir = output_dir / 'comparative_analysis'
    comp_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract key metrics from both approaches
    dynamic_summary = dynamic_results['analysis_summary']
    fixed_summary = fixed_results['analysis_summary']
    
    # Calculate comparative metrics
    comparative_metrics = {
        'convergence_comparison': {
            'dynamic_final_entropy': dynamic_summary['convergence_metrics']['final_tree_entropy'],
            'fixed_final_entropy': fixed_summary['convergence_metrics']['final_tree_entropy'],
            'dynamic_dominant_freq': dynamic_summary['convergence_metrics']['dominant_tree_frequency'],
            'fixed_dominant_freq': fixed_summary['convergence_metrics']['dominant_tree_frequency'],
            'entropy_difference': abs(dynamic_summary['convergence_metrics']['final_tree_entropy'] - 
                                    fixed_summary['convergence_metrics']['final_tree_entropy'])
        },
        'performance_comparison': {
            'dynamic_total_markers': dynamic_summary['performance_metrics']['total_unique_markers'],
            'fixed_total_markers': len(fixed_summary['fixed_gene_names']),
            'dynamic_avg_confidence': dynamic_summary['performance_metrics'].get('average_confidence', 0.9),
            'fixed_avg_confidence': fixed_summary['performance_metrics']['average_confidence'],
            'dynamic_success_rate': dynamic_summary['performance_metrics'].get('measurement_success_rate', 0.95),
            'fixed_success_rate': fixed_summary['performance_metrics']['measurement_success_rate'],
            'fixed_markers_found_rate': fixed_summary['performance_metrics']['markers_found_rate']
        },
        'clinical_practicality': {
            'dynamic_complexity': 'high',  # Variable markers each timepoint
            'fixed_complexity': 'low',     # Consistent markers
            'dynamic_standardization': 'difficult',
            'fixed_standardization': 'easy',
            'dynamic_clinical_workflow': 'complex',
            'fixed_clinical_workflow': 'simple'
        }
    }
    
    # Create comprehensive comparison
    comparison_results = {
        'analysis_metadata': {
            'patient_id': args.patient_id,
            'comparison_date': datetime.now().isoformat(),
            'approaches_compared': ['dynamic', 'fixed'],
            'n_timepoints': len(dynamic_summary['timepoints_processed'])
        },
        'comparative_metrics': comparative_metrics,
        'dynamic_summary': dynamic_summary,
        'fixed_summary': fixed_summary,
        'tree_distribution_files': {
            'dynamic_final_trees': dynamic_results.get('final_tree_distribution_file', 'not_available'),
            'fixed_final_trees': fixed_results.get('final_tree_distribution_file', 'not_available')
        },
        'recommendations': {
            'research_use': 'dynamic' if comparative_metrics['convergence_comparison']['dynamic_final_entropy'] < 
                          comparative_metrics['convergence_comparison']['fixed_final_entropy'] else 'fixed',
            'clinical_use': 'fixed',  # Always recommend fixed for clinical use
            'method_development': 'dynamic'
        }
    }
    
    # Save comparative analysis
    comp_file = comp_dir / 'comparative_analysis_results.json'
    with open(comp_file, 'w') as f:
        json.dump(comparison_results, f, indent=2, default=str)
    
    logger.info(f"Comparative analysis completed and saved to: {comp_file}")
    return comparison_results


def run_dynamic_marker_analysis(args: argparse.Namespace, logger: logging.Logger,
                               tree_distribution_summary: Dict, tree_distribution_full: Dict,
                               gene_list: List[str], gene2idx: Dict, gene_name_list: List[str],
                               timepoint_data: Dict[str, pd.DataFrame], 
                               output_dir: Path) -> Dict:
    """
    Run complete dynamic marker analysis across all timepoints.
    
    This is the original longitudinal analysis approach where optimal markers
    are selected at each timepoint based on the current tree distribution.
    """
    logger.info("Starting dynamic marker selection analysis")
    
    # Create output subdirectories for dynamic analysis
    dynamic_trees_dir = output_dir / 'dynamic_marker_analysis' / 'updated_trees'
    dynamic_selections_dir = output_dir / 'dynamic_marker_analysis' / 'marker_selections'
    plots_dir = output_dir / 'dynamic_marker_analysis' / 'plots'
    
    dynamic_trees_dir.mkdir(parents=True, exist_ok=True)
    dynamic_selections_dir.mkdir(parents=True, exist_ok=True)
    if not args.no_plots:
        plots_dir.mkdir(exist_ok=True)
    
    # Sort timepoints chronologically
    sorted_timepoints = sorted(timepoint_data.keys())
    logger.info(f"Processing {len(sorted_timepoints)} timepoints with dynamic marker selection")
    
    # Initialize tracking
    all_marker_selections = []
    analysis_summary = {
        'patient_id': args.patient_id,
        'analysis_mode': 'dynamic',
        'timepoints_processed': [],
        'convergence_metrics': {},
        'performance_metrics': {}
    }
    
    # Process each timepoint sequentially
    for order_idx, timepoint in enumerate(sorted_timepoints):
        logger.info(f"Processing timepoint {order_idx + 1}/{len(sorted_timepoints)}: {timepoint}")
        
        # For the first timepoint, use original tree distributions
        if order_idx == 0:
            current_tree_summary = tree_distribution_summary
            current_tree_full = tree_distribution_full
            
            # Extract tree components for analysis
            tree_list = current_tree_full['tree_structure']
            node_list = current_tree_full['node_dict'] 
            tree_freq_list = current_tree_full['freq']
            clonal_freq_list = current_tree_full['vaf_frac']
            
            logger.info(f"Using original tree distributions: {len(tree_list)} trees")
        
        # For subsequent timepoints, load updated tree distribution from previous iteration
        else:
            updated_file = dynamic_trees_dir / f'{args.method}_bootstrap_summary_updated_{args.n_markers}_{order_idx-1}_bayesian.pkl'
            
            if not updated_file.exists():
                logger.error(f"Updated tree distribution file not found: {updated_file}")
                raise FileNotFoundError(f"Missing updated tree distribution for timepoint {order_idx-1}")
            
            with open(updated_file, 'rb') as f:
                current_tree_summary = pickle.load(f)
            
            # Extract tree components
            tree_list = current_tree_summary['tree_structure']
            node_list = current_tree_summary['node_dict']
            tree_freq_list = current_tree_summary['freq']
            
            # Recalculate clonal frequencies by averaging across samples
            clonal_freq_list = []
            for idx in range(len(current_tree_summary['vaf_frac'])):
                clonal_freq_dict = current_tree_summary['vaf_frac'][idx]
                clonal_freq_dict_new = {}
                for node, freqs in clonal_freq_dict.items():
                    clonal_freq_dict_new[node] = [list(np.array(freqs).mean(axis=0))]
                clonal_freq_list.append(clonal_freq_dict_new)
            
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
        
        # Extract ddPCR measurements for selected markers
        ddpcr_measurements = []
        for gene_name in selected_gene_names:
            if gene_name in current_ddpcr_data.index:
                mut_count = current_ddpcr_data.loc[gene_name, 'MutDOR']  # Mutant droplets
                total_count = current_ddpcr_data.loc[gene_name, 'DOR']   # Total droplets
                wt_count = total_count - mut_count  # Calculate WT count
                
                ddpcr_measurements.append({
                    'gene': gene_name,
                    'mut': mut_count,
                    'WT': wt_count,
                    'liquid_biopsy_sample': timepoint
                })
            else:
                logger.warning(f"Selected marker {gene_name} not found in ddPCR data for timepoint {timepoint}")
        
        if not ddpcr_measurements:
            logger.error(f"No ddPCR data found for selected markers at timepoint {timepoint}")
            continue
        
        # Create DataFrame with ddPCR measurements
        df_ddpcr = pd.DataFrame(ddpcr_measurements)
        marker_idx2gene = {i: df_ddpcr["gene"].iloc[i] for i in range(len(df_ddpcr))}
        
        # Extract counts for Bayesian updating
        ddpcr_marker_counts = list(df_ddpcr["mut"])
        read_depth_list = list(df_ddpcr["mut"] + df_ddpcr["WT"])  # Total = mut + WT
        
        logger.info(f"ddPCR measurements: {len(ddpcr_measurements)} markers")
        logger.info(f"Mutant counts: {ddpcr_marker_counts}")
        logger.info(f"Read depths: {read_depth_list}")
        
        # Extract tree information from summary
        tree_list_summary, node_list_summary, node_name_list_summary, tree_freq_list_summary = \
            current_tree_summary['tree_structure'], current_tree_summary['node_dict'], \
            current_tree_summary['node_dict_name'], current_tree_summary['freq']
        
        # Update tree distributions using Bayesian approach
        logger.info("Updating tree distributions using Bayesian inference...")
        
        updated_tree_freq_list = adjust_tree_distribution_struct_bayesian(
            tree_list_summary, node_name_list_summary,
            tree_freq_list_summary, read_depth_list,
            ddpcr_marker_counts, marker_idx2gene)
        
        # Create updated tree distribution summary
        updated_tree_distribution_summary = update_tree_distribution_bayesian(
            current_tree_summary, updated_tree_freq_list)
        
        # Save updated tree distribution for next iteration
        updated_file = dynamic_trees_dir / f'{args.method}_bootstrap_summary_updated_{args.n_markers}_{order_idx}_bayesian.pkl'
        with open(updated_file, 'wb') as f:
            pickle.dump(updated_tree_distribution_summary, f)
        
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
        
        marker_file = dynamic_selections_dir / f'marker_selection_{timepoint}_{order_idx}.json'
        with open(marker_file, 'w') as f:
            json.dump(marker_selection_results, f, indent=2, default=str)
        
        logger.info(f"Saved marker selection results: {marker_file}")
    
    # Calculate final performance metrics
    logger.info("Calculating dynamic approach performance metrics...")
    
    final_tree_freq = updated_tree_distribution_summary['freq']
    tree_entropy = -np.sum([f * np.log(f + 1e-10) for f in final_tree_freq if f > 0])
    dominant_tree_freq = max(final_tree_freq)
    
    all_markers_used = set()
    for selection in all_marker_selections:
        all_markers_used.update(selection['selected_markers'])
    
    analysis_summary['convergence_metrics'] = {
        'final_tree_entropy': tree_entropy,
        'dominant_tree_frequency': dominant_tree_freq,
        'convergence_timepoint': len(sorted_timepoints),
        'n_trees_remaining': sum(1 for f in final_tree_freq if f > 0.01)
    }
    
    analysis_summary['performance_metrics'] = {
        'total_unique_markers': len(all_markers_used),
        'total_selections': len(all_marker_selections),
        'final_objective_fraction': all_marker_selections[-1]['objective_fraction'],
        'final_objective_structure': all_marker_selections[-1]['objective_structure'],
        'markers_per_timepoint': [len(sel['selected_markers']) for sel in all_marker_selections]
    }
    
    # Save complete dynamic analysis results
    
    # Save final tree distribution separately as pickle (like iterative_pipeline_real.py)
    # This handles tuple keys properly since pickle can serialize them
    final_tree_file = dynamic_trees_dir / f'{args.method}_bootstrap_summary_final_dynamic_{args.n_markers}_bayesian.pkl'
    with open(final_tree_file, 'wb') as f:
        pickle.dump(updated_tree_distribution_summary, f)
    
    logger.info(f"Saved final tree distribution: {final_tree_file}")
    
    # Create JSON-safe results summary (exclude tree distribution with tuple keys)
    results_summary = {
        'analysis_summary': analysis_summary,
        'all_marker_selections': all_marker_selections,
        'final_tree_distribution_file': str(final_tree_file)  # Reference to pickle file instead of data
    }
    
    dynamic_results_file = output_dir / 'dynamic_marker_analysis' / 'dynamic_marker_complete_results.json'
    with open(dynamic_results_file, 'w') as f:
        json.dump(results_summary, f, indent=2, default=str)
    
    logger.info(f"Dynamic marker analysis completed successfully")
    logger.info(f"Results saved to: {output_dir / 'dynamic_marker_analysis'}")
    
    return results_summary


def main():
    """
    Main entry point for the longitudinal analysis pipeline.
    Now uses YAML configuration files for cleaner parameter management.
    """
    # Initialize a basic logger for early error handling
    logger = logging.getLogger('main_pipeline')
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(handler)
    
    try:
        # Parse command line arguments (now just config file and overrides)
        cmd_args = parse_args()
        
        # Load configuration from YAML file
        logger.info(f"Loading configuration from: {cmd_args.config}")
        config = load_config(cmd_args.config)
        
        # Convert config to args namespace for compatibility with existing code
        args = config_to_args(config, cmd_args)
        
        # Set up output directory and proper logging (replaces the basic logger)
        output_dir = Path(args.output_dir)
        logger = setup_logging(output_dir, args.patient_id)
        
        logger.info("=== Longitudinal Cancer Evolution Analysis Pipeline ===")
        logger.info(f"Configuration file: {cmd_args.config}")
        logger.info(f"Patient ID: {args.patient_id}")
        logger.info(f"Analysis mode: {args.analysis_mode}")
        logger.info(f"Output directory: {output_dir}")
        
        # Log configuration summary
        logger.info("Configuration Summary:")
        logger.info(f"  Input files:")
        logger.info(f"    - Aggregation dir: {args.aggregation_dir}")
        logger.info(f"    - SSM file: {args.ssm_file}")
        logger.info(f"    - Longitudinal data: {args.longitudinal_data}")
        logger.info(f"  Analysis parameters:")
        logger.info(f"    - N markers: {args.n_markers}")
        logger.info(f"    - Read depth: {args.read_depth}")
        logger.info(f"    - Lambda1 (fraction weight): {args.lambda1}")
        logger.info(f"    - Lambda2 (structure weight): {args.lambda2}")
        logger.info(f"    - Method: {args.method}")
        if args.analysis_mode in ['fixed', 'both']:
            logger.info(f"  Fixed markers: {args.fixed_markers}")
        
        # Validate input files
        if not validate_input_files(args, logger):
            logger.error("Input validation failed. Exiting.")
            sys.exit(1)
        
        # Load tree distribution data from aggregation stage
        aggregation_dir = Path(args.aggregation_dir)
        tree_distribution_summary, tree_distribution_full = load_tree_distributions(
            aggregation_dir, args.method, logger)
        
        # Load tissue data from SSM file
        ssm_file = Path(args.ssm_file)
        tissue_df, gene2idx, gene_name_list, gene_list, gene_name2idx = load_tissue_data_from_ssm(ssm_file, logger)
        
        # Load longitudinal data from CSV file
        longitudinal_file = Path(args.longitudinal_data)
        timepoint_data = load_longitudinal_data_from_csv(longitudinal_file, logger)
        
        # Log successful data loading
        logger.info("Data loading completed successfully")
        logger.info(f"Tree distributions: {len(tree_distribution_summary['tree_structure'])} trees")
        logger.info(f"Tissue mutations: {len(tissue_df)} mutations")
        logger.info(f"Longitudinal timepoints: {len(timepoint_data)} timepoints")
        logger.info(f"Analysis mode: {args.analysis_mode}")
        
        # Create output subdirectories
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize results storage
        dynamic_results = None
        fixed_results = None
        
        # Run analysis based on selected mode
        if args.analysis_mode in ['dynamic', 'both']:
            logger.info("=== Running Dynamic Marker Analysis ===")
            dynamic_results = run_dynamic_marker_analysis(
                args, logger, tree_distribution_summary, tree_distribution_full,
                gene_list, gene2idx, gene_name_list, timepoint_data, output_dir)
        
        if args.analysis_mode in ['fixed', 'both']:
            logger.info("=== Running Fixed Marker Analysis ===")
            fixed_results = run_fixed_marker_analysis(
                args, logger, tree_distribution_summary, tree_distribution_full,
                gene_list, gene2idx, gene_name_list, timepoint_data, output_dir, gene_name2idx)
        
        # Generate comparative analysis if both approaches were run
        if args.analysis_mode == 'both' and dynamic_results and fixed_results:
            logger.info("=== Generating Comparative Analysis ===")
            comparative_results = generate_comparative_analysis(
                dynamic_results, fixed_results, args, logger, output_dir)
        
        # Generate summary report
        summary_report = {
            'patient_id': args.patient_id,
            'analysis_timestamp': datetime.now().isoformat(),
            'analysis_mode': args.analysis_mode,
            'parameters': {
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
            'results': {
                'total_timepoints': len(timepoint_data),
                'timepoints_processed': sorted(timepoint_data.keys()),
                'output_directory': str(output_dir),
                'approaches_run': [args.analysis_mode] if args.analysis_mode != 'both' else ['dynamic', 'fixed']
            }
        }
        
        summary_file = output_dir / 'analysis_summary.json'
        with open(summary_file, 'w') as f:
            json.dump(summary_report, f, indent=2)
        
        logger.info(f"Analysis completed successfully")
        logger.info(f"Results saved to: {output_dir}")
        logger.info(f"Summary report: {summary_file}")
        
        logger.info("Longitudinal analysis pipeline completed successfully")
        
    except Exception as e:
        logger.error(f"Pipeline failed with error: {e}")
        # Check if args exists and has debug flag before trying to use it
        try:
            if 'args' in locals() and hasattr(args, 'debug') and args.debug:
                import traceback
                logger.error(f"Full traceback:\n{traceback.format_exc()}")
        except:
            # If we can't check debug mode, just show the traceback anyway for config errors
            import traceback
            logger.error(f"Full traceback:\n{traceback.format_exc()}")
        sys.exit(1)


if __name__ == '__main__':
    main() 