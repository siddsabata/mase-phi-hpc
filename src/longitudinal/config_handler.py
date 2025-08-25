#!/usr/bin/env python3
"""
Configuration handling for longitudinal cancer evolution analysis.

This module handles YAML configuration loading, validation, and conversion
to the expected format for the longitudinal analysis pipeline.

Authors: TracerX Pipeline Development Team
"""

import argparse
import yaml
import logging
from pathlib import Path
from typing import Dict, Any

logger = logging.getLogger(__name__)


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
            'focus_sample': 0,
            'track_clone_freq': True  # Enable clone frequency tracking
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
    Simplified to use YAML configuration files.
    
    Returns:
        Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        description='Run longitudinal cancer evolution analysis with iterative Bayesian tree updating',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Run fixed marker analysis using YAML configuration
    python longitudinal_main.py --config configs/cruk0044_fixed_markers.yaml
    
    # Run dynamic marker analysis using YAML configuration  
    python longitudinal_main.py --config configs/cruk0044_dynamic.yaml
    
    # Run with debug mode enabled
    python longitudinal_main.py --config configs/cruk0044_fixed_markers.yaml --debug
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
    
    parser.add_argument('--plot-clonal-freq', action='store_true',
                       help='Enable clone frequency tracking and plotting (overrides config file setting)')
    
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
    
    # Validate analysis mode (only 'dynamic' or 'fixed' allowed now)
    if args.analysis_mode not in ['dynamic', 'fixed']:
        raise ValueError(f"Invalid analysis mode: {args.analysis_mode}. Must be 'dynamic' or 'fixed'")
    
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
    args.track_clone_freq = params['track_clone_freq']
    
    # Fixed markers (only used in fixed mode)
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
    if hasattr(cmd_args, 'plot_clonal_freq') and cmd_args.plot_clonal_freq:
        args.track_clone_freq = True
    
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
    if args.analysis_mode == 'fixed':
        if not args.fixed_markers:
            logger.error("Fixed markers must be specified when using fixed analysis mode")
            logger.error("Add fixed_markers to YAML config, e.g., fixed_markers: ['TP53', 'KRAS', 'PIK3CA']")
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