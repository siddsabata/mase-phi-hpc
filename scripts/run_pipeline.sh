#!/bin/bash
# Master Pipeline Script for TracerX Marker Selection Pipeline (Steps 1-4)
# This script runs on the login node and submits individual SLURM jobs for each pipeline step
# Usage: bash run_pipeline.sh config.yaml [--dry-run]
# Example: bash run_pipeline.sh configs/analysis/test_analysis.yaml
# Example: bash run_pipeline.sh configs/analysis/test_analysis.yaml --dry-run

set -e  # Exit on any error

# --- Setup conda environment for YAML parsing ---
echo "Using conda environment for YAML parsing..."
cd "$( dirname "${BASH_SOURCE[0]}" )/.."  # Change to repo root

# --- Input Validation ---
if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <config.yaml> [--dry-run]"
    echo "Example: bash $0 configs/standard_analysis.yaml"
    echo "Example: bash $0 configs/test_analysis.yaml --dry-run"
    exit 1
fi

CONFIG_FILE=$1
DRY_RUN=${2:-""}

# Validate dry-run option
if [ -n "$DRY_RUN" ] && [ "$DRY_RUN" != "--dry-run" ]; then
    echo "Error: Invalid option '$DRY_RUN'. Only '--dry-run' is supported."
    exit 1
fi

# --- Get script directory for absolute paths ---
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
echo "Pipeline script directory: ${SCRIPT_DIR}"

# --- Validate Configuration File ---
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file not found: $CONFIG_FILE"
    exit 1
fi

# Convert to absolute path
if [[ ! "$CONFIG_FILE" = /* ]]; then
    CONFIG_FILE="${SCRIPT_DIR}/${CONFIG_FILE}"
    echo "Converted config file path to absolute: ${CONFIG_FILE}"
fi

# --- Load Configuration Using Python ---
echo "Loading configuration from: ${CONFIG_FILE}"

# Create temporary Python script to parse YAML and export variables
TEMP_VARS_FILE=$(mktemp)
cat > /tmp/parse_config.py << 'EOF'
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
    array_limit = phylowgs_config.get('parallel_limit', 10)
    
    # Marker selection configuration
    marker_config = config.get('marker_selection', {})
    read_depth = marker_config.get('read_depth', 1500)
    filter_strategy = marker_config.get('filter_strategy', 'any_high')
    filter_threshold = marker_config.get('filter_threshold', 0.9)
    
    # HPC configuration
    hpc_config = config.get('hpc', {})
    
    # Print shell variable exports
    print(f'export PATIENT_ID="{patient_id}"')
    print(f'export INPUT_SSM_FILE="{ssm_file}"')
    print(f'export CODE_DIR="{code_dir}"')
    print(f'export PATIENT_BASE_DIR="{base_dir}"')
    print(f'export NUM_BOOTSTRAPS="{num_bootstraps}"')
    print(f'export ARRAY_LIMIT="{array_limit}"')
    print(f'export READ_DEPTH="{read_depth}"')
    print(f'export FILTER_STRATEGY="{filter_strategy}"')
    print(f'export FILTER_THRESHOLD="{filter_threshold}"')
    
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
EOF

# Parse configuration and source the variables
conda run -n mase_phi_hpc python /tmp/parse_config.py "$CONFIG_FILE" > "$TEMP_VARS_FILE"
if [ $? -ne 0 ]; then
    echo "Error: Failed to parse configuration file"
    rm -f "$TEMP_VARS_FILE" /tmp/parse_config.py
    exit 1
fi

source "$TEMP_VARS_FILE"
rm -f "$TEMP_VARS_FILE" /tmp/parse_config.py

# --- Validate Required Configuration ---
if [ -z "$PATIENT_ID" ] || [ -z "$INPUT_SSM_FILE" ] || [ -z "$CODE_DIR" ] || [ -z "$PATIENT_BASE_DIR" ]; then
    echo "Error: Missing required configuration parameters:"
    echo "  PATIENT_ID: '$PATIENT_ID'"
    echo "  INPUT_SSM_FILE: '$INPUT_SSM_FILE'"
    echo "  CODE_DIR: '$CODE_DIR'"
    echo "  PATIENT_BASE_DIR: '$PATIENT_BASE_DIR'"
    exit 1
fi

# --- Convert Paths to Absolute ---
if [[ ! "$INPUT_SSM_FILE" = /* ]]; then
    INPUT_SSM_FILE="${SCRIPT_DIR}/${INPUT_SSM_FILE}"
    echo "Converted input SSM file path to absolute: ${INPUT_SSM_FILE}"
fi

if [[ ! "$CODE_DIR" = /* ]]; then
    CODE_DIR="${SCRIPT_DIR}/${CODE_DIR}"
    echo "Converted code directory path to absolute: ${CODE_DIR}"
fi

if [[ ! "$PATIENT_BASE_DIR" = /* ]]; then
    PATIENT_BASE_DIR="${SCRIPT_DIR}/${PATIENT_BASE_DIR}"
    echo "Converted patient base directory path to absolute: ${PATIENT_BASE_DIR}"
fi

# --- Normalize PATIENT_BASE_DIR (remove trailing slashes) ---
PATIENT_BASE_DIR=$(echo "${PATIENT_BASE_DIR}" | sed 's:/*$::')

# --- Validate Input Files/Directories ---
if [ ! -f "$INPUT_SSM_FILE" ]; then
    echo "Error: Input SSM file not found: $INPUT_SSM_FILE"
    exit 1
fi

if [ ! -d "$CODE_DIR" ]; then
    echo "Error: Code directory not found: $CODE_DIR"
    exit 1
fi

# --- Define Directory Structure ---
INITIAL_DIR="${PATIENT_BASE_DIR}/initial"
BOOTSTRAP_STAGE_OUTPUT_DIR="${INITIAL_DIR}"
BOOTSTRAPS_DATA_DIR="${BOOTSTRAP_STAGE_OUTPUT_DIR}/bootstraps"
AGGREGATION_RESULTS_DIR="${BOOTSTRAP_STAGE_OUTPUT_DIR}/aggregation_results"
MARKERS_DIR="${BOOTSTRAP_STAGE_OUTPUT_DIR}/markers"
LOG_DIR="${INITIAL_DIR}/logs"

# --- Create Required Directories ---
if [ "$DRY_RUN" != "--dry-run" ]; then
    mkdir -p "${PATIENT_BASE_DIR}"
    mkdir -p "${INITIAL_DIR}"
    mkdir -p "${BOOTSTRAPS_DATA_DIR}"
    mkdir -p "${AGGREGATION_RESULTS_DIR}"
    mkdir -p "${MARKERS_DIR}"
    mkdir -p "${LOG_DIR}"
fi

# --- Logging Setup ---
MASTER_LOG="${LOG_DIR}/pipeline_master.log"
if [ "$DRY_RUN" != "--dry-run" ]; then
    exec > >(tee -a "${MASTER_LOG}") 2>&1
fi

echo "=== Pipeline Start: $(date) ==="
echo "Configuration File: ${CONFIG_FILE}"
echo "Patient ID: ${PATIENT_ID}"
echo "Input SSM: ${INPUT_SSM_FILE}"
echo "Patient Base Directory: ${PATIENT_BASE_DIR}"
echo "Initial Processing Directory: ${INITIAL_DIR}"
echo "Number of Bootstraps: ${NUM_BOOTSTRAPS}"
echo "Read Depth: ${READ_DEPTH}"
echo "Code Directory: ${CODE_DIR}"
echo "Log Directory: ${LOG_DIR}"
echo "Dry Run Mode: ${DRY_RUN:-false}"
echo "----------------------------------------"

# --- Function to Submit Job with YAML Configuration ---
submit_job_yaml() {
    local job_name=$1
    local dependency=$2
    local script_path=$3
    shift 3
    local args=("$@")

    local script_containing_dir
    script_containing_dir=$(dirname "$script_path")
    
    local dependency_flag=""
    if [ ! -z "$dependency" ]; then
        dependency_flag="--dependency=afterok:${dependency}"
    fi
    
    # Get step-specific HPC configuration
    local step_upper=$(echo "${job_name}" | tr '[:lower:]' '[:upper:]')
    local partition_var="${step_upper}_PARTITION"
    local cpus_var="${step_upper}_CPUS"
    local memory_var="${step_upper}_MEMORY"
    local walltime_var="${step_upper}_WALLTIME"
    
    local partition=${!partition_var:-pool1}
    local cpus=${!cpus_var:-1}
    local memory=${!memory_var:-8G}
    local walltime=${!walltime_var:-02:00:00}
    
    echo "Submitting ${job_name} with configuration:" >&2
    echo "  Partition: ${partition}" >&2
    echo "  CPUs: ${cpus}" >&2
    echo "  Memory: ${memory}" >&2
    echo "  Walltime: ${walltime}" >&2
    echo "  Script: ${script_path}" >&2
    
    # Build sbatch command array
    local sbatch_cmd_array=(
        sbatch --parsable
        --chdir="${script_containing_dir}"
        --job-name="${PATIENT_ID}_${job_name}"
        --partition="${partition}"
        --cpus-per-task="${cpus}"
        --mem="${memory}"
        --time="${walltime}"
        --output="${LOG_DIR}/slurm_${job_name}_%j.out"
        --error="${LOG_DIR}/slurm_${job_name}_%j.err"
    )
    
    # Add dependency if specified
    if [ ! -z "$dependency_flag" ]; then
        # For jobs that depend on array jobs (like aggregation depending on phylowgs),
        # use 'afterany' instead of 'afterok' to allow partial bootstrap failures
        if [ "$job_name" = "aggregation" ] || [ "$job_name" = "marker_selection" ]; then
            # These jobs depend on array jobs where some elements may fail (bootstraps)
            # Use afterany to proceed when array job completes (regardless of individual element success)
            local modified_dependency=$(echo "$dependency_flag" | sed 's/afterok/afterany/')
            sbatch_cmd_array+=("$modified_dependency")
        else
            sbatch_cmd_array+=("$dependency_flag")
        fi
    fi
    
    # Add array configuration for phylowgs
    if [ "$job_name" == "phylowgs" ]; then
        local array_max=$((NUM_BOOTSTRAPS - 1))
        sbatch_cmd_array+=("--array=0-${array_max}%${ARRAY_LIMIT}")
    fi
    
    sbatch_cmd_array+=("${script_path}")
    sbatch_cmd_array+=("${args[@]}")
    
    if [ "$DRY_RUN" == "--dry-run" ]; then
        echo "DRY RUN: Would execute: ${sbatch_cmd_array[*]}" >&2
        echo "12345"  # Return fake job ID for dry run
        return 0
    fi
    
    # Execute sbatch command
    local job_id_output
    job_id_output=$("${sbatch_cmd_array[@]}")
    local sbatch_status=$?

    if [ ${sbatch_status} -ne 0 ] || [ -z "${job_id_output}" ]; then
        echo "Error: Failed to submit ${job_name} job. sbatch command exited with status ${sbatch_status}." >&2
        echo "Script path attempted: ${script_path}" >&2
        exit 1 
    fi
    
    local job_id="${job_id_output}"
    echo "${job_name} job submitted with ID: ${job_id}" >&2
    echo "${job_id}"
}

# --- Pipeline Execution ---

# Step 1: Bootstrap
echo "Starting Bootstrap Stage..."
BOOTSTRAP_JOB_ID=$(submit_job_yaml "bootstrap" "" \
    "${CODE_DIR}/src/bootstrap/bootstrap.sh" \
    "${INPUT_SSM_FILE}" "${BOOTSTRAP_STAGE_OUTPUT_DIR}" "${CODE_DIR}" "${NUM_BOOTSTRAPS}")

# Step 2: PhyloWGS (Array Job)
echo "Starting PhyloWGS Stage..."
PHYLOWGS_JOB_ID=$(submit_job_yaml "phylowgs" "${BOOTSTRAP_JOB_ID}" \
    "${CODE_DIR}/src/phylowgs/phylowgs.sh" \
    "${BOOTSTRAPS_DATA_DIR}" "${CODE_DIR}")

# Step 3: Aggregation
echo "Starting Aggregation Stage..."
AGGREGATION_JOB_ID=$(submit_job_yaml "aggregation" "${PHYLOWGS_JOB_ID}" \
    "${CODE_DIR}/src/aggregation/aggregation.sh" \
    "${PATIENT_ID}" "${BOOTSTRAPS_DATA_DIR}" "${AGGREGATION_RESULTS_DIR}" "${CODE_DIR}")

# Step 4: Marker Selection
echo "Starting Marker Selection Stage..."
MARKER_SELECTION_JOB_ID=$(submit_job_yaml "marker_selection" "${AGGREGATION_JOB_ID}" \
    "${CODE_DIR}/src/markers/marker_selection.sh" \
    "${PATIENT_ID}" "${AGGREGATION_RESULTS_DIR}" "${INPUT_SSM_FILE}" "${CODE_DIR}" \
    "${READ_DEPTH}" "${FILTER_STRATEGY}" "${FILTER_THRESHOLD}")

# --- Final Status ---
echo "----------------------------------------"
if [ "$DRY_RUN" == "--dry-run" ]; then
    echo "DRY RUN COMPLETE: Pipeline configuration validated successfully!"
    echo "Configuration file: ${CONFIG_FILE}"
    echo "All job submissions would succeed with provided configuration."
else
    echo "Pipeline submitted successfully!"
    echo "Final marker selection job ID: ${MARKER_SELECTION_JOB_ID}"
    echo "Monitor progress with: squeue -u $USER"
    echo "Check logs in: ${LOG_DIR}"
    
    # Save Job IDs for Reference
    cat > "${LOG_DIR}/job_ids.txt" << EOF
Pipeline Job IDs for ${PATIENT_ID} (Config: ${CONFIG_FILE}):
Bootstrap: ${BOOTSTRAP_JOB_ID}
PhyloWGS: ${PHYLOWGS_JOB_ID}
Aggregation: ${AGGREGATION_JOB_ID}
Marker Selection: ${MARKER_SELECTION_JOB_ID}
EOF
fi

echo "=== Pipeline Submission Complete: $(date) ==="