#!/usr/bin/env python3
"""
Main entry point for longitudinal cancer evolution analysis.

This script orchestrates the modular longitudinal analysis pipeline,
supporting both dynamic and fixed marker approaches with clean,
focused functionality.

Usage:
    python longitudinal_main.py --config configs/cruk0044_fixed_markers.yaml
    python longitudinal_main.py --config configs/cruk0044_dynamic.yaml

Authors: TracerX Pipeline Development Team
Version: 2.0.0
"""

import sys
import logging
from pathlib import Path

# Import our modular components
from config_handler import parse_args, load_config, config_to_args, validate_input_files
from data_loader import load_tree_distributions, load_tissue_data_from_ssm, load_longitudinal_data_from_csv
from fixed_analysis import run_fixed_marker_analysis
from dynamic_analysis import run_dynamic_marker_analysis
from output_manager import setup_logging, save_final_report


def main():
    """
    Main entry point for the longitudinal analysis pipeline.
    
    This simplified version focuses on core functionality:
    - Load configuration from YAML files
    - Validate inputs
    - Load data
    - Run either fixed or dynamic analysis (not both)
    - Save results in standardized format
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
        
        logger.info("=== Longitudinal Cancer Evolution Analysis Pipeline v2.0 ===")
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
        logger.info(f"    - Method: {args.method}")
        if args.analysis_mode == 'fixed':
            logger.info(f"  Fixed markers: {args.fixed_markers}")
        
        # Validate input files
        if not validate_input_files(args, logger):
            logger.error("Input validation failed. Exiting.")
            sys.exit(1)
        
        # Load tree distribution data from aggregation stage
        logger.info("=== Loading Data ===")
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
        
        # Run analysis based on selected mode (simplified - only one mode at a time)
        logger.info(f"=== Running {args.analysis_mode.title()} Analysis ===")
        
        if args.analysis_mode == 'fixed':
            results_summary = run_fixed_marker_analysis(
                args, logger, tree_distribution_summary, tree_distribution_full,
                gene_list, gene2idx, gene_name_list, timepoint_data, 
                output_dir, gene_name2idx)
        
        elif args.analysis_mode == 'dynamic':
            results_summary = run_dynamic_marker_analysis(
                args, logger, tree_distribution_summary, tree_distribution_full,
                gene_list, gene2idx, gene_name_list, timepoint_data, output_dir)
        
        else:
            logger.error(f"Invalid analysis mode: {args.analysis_mode}")
            sys.exit(1)
        
        # Generate final report
        logger.info("=== Generating Final Report ===")
        report_file = save_final_report(output_dir, args, results_summary, timepoint_data)
        
        logger.info(f"Analysis completed successfully")
        logger.info(f"Results saved to: {output_dir}")
        logger.info(f"Final report: {report_file}")
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