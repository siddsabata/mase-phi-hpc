#!/usr/bin/env python3
"""
YAML configuration parser for TracerX Marker Selection Pipeline
Extracts configuration values and exports them as shell variables
"""

import yaml
import sys
import os

def parse_config(config_file):
    """Parse YAML configuration and export as shell variables."""
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    # Extract configuration values with defaults
    patient_id = config.get('patient_id', 'UNKNOWN')
    
    # Input configuration
    input_config = config.get('input', {})
    ssm_file = input_config.get('ssm_file', '')
    code_dir = input_config.get('code_dir', '')
    
    # Output configuration
    output_config = config.get('output', {})
    base_dir = output_config.get('base_dir', '')
    
    # Bootstrap configuration
    bootstrap_config = config.get('bootstrap', {})
    num_bootstraps = bootstrap_config.get('num_bootstraps', 100)
    
    # PhyloWGS configuration
    phylowgs_config = config.get('phylowgs', {})
    num_chains = phylowgs_config.get('num_chains', 5)
    array_limit = phylowgs_config.get('parallel_limit', 10)
    
    # Marker selection configuration
    marker_config = config.get('marker_selection', {})
    read_depth = marker_config.get('read_depth', 1500)
    # Note: VAF filtering parameters removed - filtering now handled in bootstrap stage
    
    # HPC configuration
    hpc_config = config.get('hpc', {})
    
    # Construct patient-specific base directory and filtered SSM file path
    patient_base_dir = f"{base_dir}/{patient_id}"
    filtered_ssm_file = f"{patient_base_dir}/initial/ssm_filtered.txt"
    
    # Print shell variable exports
    print(f'export PATIENT_ID="{patient_id}"')
    print(f'export INPUT_SSM_FILE="{ssm_file}"')
    print(f'export FILTERED_SSM_FILE="{filtered_ssm_file}"')
    print(f'export CODE_DIR="{code_dir}"')
    print(f'export PATIENT_BASE_DIR="{patient_base_dir}"')
    print(f'export NUM_BOOTSTRAPS="{num_bootstraps}"')
    print(f'export NUM_CHAINS="{num_chains}"')
    print(f'export ARRAY_LIMIT="{array_limit}"')
    print(f'export READ_DEPTH="{read_depth}"')
    # Note: Filter strategy and threshold exports removed - handled in bootstrap stage
    
    # HPC settings for each step
    for step in ['bootstrap', 'phylowgs', 'aggregation', 'marker_selection']:
        step_config = hpc_config.get(step, {})
        step_upper = step.upper()
        print(f'export {step_upper}_PARTITION="{step_config.get("partition", "pool1")}"')
        print(f'export {step_upper}_CPUS="{step_config.get("cpus_per_task", 1)}"')
        print(f'export {step_upper}_MEMORY="{step_config.get("memory", "8G")}"')
        print(f'export {step_upper}_WALLTIME="{step_config.get("walltime", "02:00:00")}"')
        print(f'export {step_upper}_CONDA_ENV="{step_config.get("conda_env", "base")}"')
    
    # Special handling for modules
    marker_modules = hpc_config.get('marker_selection', {}).get('modules', [])
    if marker_modules:
        print(f'export MARKER_SELECTION_MODULES="{" ".join(marker_modules)}"')

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python parse_config.py <config.yaml>", file=sys.stderr)
        sys.exit(1)
    
    try:
        parse_config(sys.argv[1])
    except Exception as e:
        print(f"Error parsing configuration: {e}", file=sys.stderr)
        sys.exit(1)