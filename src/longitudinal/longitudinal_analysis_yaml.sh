#!/bin/bash
# Longitudinal cancer evolution analysis using YAML configuration
# Usage: sbatch longitudinal_analysis_yaml.sh <config_yaml_file> [additional_flags] [slurm_log_suffix]

set -e

# --- Argument Parsing and Validation ---
if [ "$#" -lt 1 ] || [ "$#" -gt 3 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: sbatch $0 <config_yaml_file> [additional_flags] [slurm_log_suffix]"
    echo ""
    echo "Examples:"
    echo "  sbatch $0 configs/cruk0044_longitudinal.yaml"
    echo "  sbatch $0 configs/cruk0044_longitudinal.yaml '--debug'"
    echo "  sbatch $0 configs/cruk0044_longitudinal.yaml '--debug --no-plots' 'test_run'"
    echo ""
    echo "Configuration files should be in YAML format with all required parameters."
    echo "The code directory is now specified in the YAML config file."
    echo "Uses unified longitudinal pipeline (v3.0) - handles all marker selection approaches"
    exit 1
fi

CONFIG_FILE=$1
ADDITIONAL_FLAGS=${2:-""} # Optional additional flags like --debug, --no-plots
LOG_SUFFIX=${3:-""} # Optional suffix for log files

# Validate required input files and directories
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file '$CONFIG_FILE' not found."
    exit 1
fi

# Extract key parameters from config for validation and logging setup
# This is a simple extraction - the Python script will do full validation
PATIENT_ID=$(grep "^patient_id:" "$CONFIG_FILE" | sed 's/patient_id: *"\?\([^"]*\)"\?/\1/' | tr -d '"')
OUTPUT_BASE=$(grep "base_dir:" "$CONFIG_FILE" | sed 's/.*base_dir: *"\?\([^"]*\)"\?/\1/' | tr -d '"')
CODE_DIR=$(grep "code_dir:" "$CONFIG_FILE" | sed 's/.*code_dir: *"\?\([^"]*\)"\?/\1/' | tr -d '"')

# Expand tilde paths
OUTPUT_BASE="${OUTPUT_BASE/#\~/$HOME}"
CODE_DIR="${CODE_DIR/#\~/$HOME}"

if [ ! -d "$CODE_DIR" ]; then
    echo "Error: Code directory '$CODE_DIR' (from config) not found."
    exit 1
fi

if [ -z "$PATIENT_ID" ]; then
    echo "Error: Could not extract patient_id from config file."
    exit 1
fi

if [ -z "$OUTPUT_BASE" ]; then
    echo "Error: Could not extract output base_dir from config file."
    exit 1
fi

if [ -z "$CODE_DIR" ]; then
    echo "Error: Could not extract code_dir from config file."
    exit 1
fi

# No analysis mode validation needed - unified pipeline handles both approaches

echo "--- Longitudinal YAML Analysis Script v3.0 Start ---"
echo "Configuration file: $CONFIG_FILE"
echo "Patient ID: $PATIENT_ID"
echo "Code directory: $CODE_DIR (from config)"
echo "Output base directory: $OUTPUT_BASE"
echo "Additional flags: $ADDITIONAL_FLAGS"

# Load Gurobi module
echo "Loading Gurobi 11.0.2 module..."
module load gurobi1102
if [ $? -ne 0 ]; then
    echo "Error: Failed to load gurobi1102 module. Exiting."
    exit 1
fi
echo "Gurobi module loaded successfully."

# --- Setup log directory ---
# Create logs directory in the output base directory
LOG_DIR="${OUTPUT_BASE}/logs"
mkdir -p "${LOG_DIR}"

echo "--- Unified Longitudinal Analysis Script v3.0 Execution ---"
echo "Job ID: $SLURM_JOB_ID"
echo "Configuration file: $CONFIG_FILE"
echo "Patient ID: $PATIENT_ID"
echo "Code directory: $CODE_DIR"
echo "Output base directory: $OUTPUT_BASE"
echo "Additional flags: $ADDITIONAL_FLAGS"
echo "---------------------------------------"

# --- Environment Setup ---
echo "Using conda environment mase_phi_hpc..."
cd "$CODE_DIR"

# Verify Gurobi is accessible from Python
echo "Verifying Gurobi is accessible from Python..."
conda run -n mase_phi_hpc python -c "import gurobipy; print(f'Gurobi version: {gurobipy.gurobi.version()}')"
if [ $? -ne 0 ]; then
    echo "Error: Failed to import gurobipy or access Gurobi."
    exit 1
fi
echo "Gurobi verification successful."

# Verify required Python packages
echo "Verifying required Python packages..."
conda run -n mase_phi_hpc python -c "import pandas, numpy, matplotlib, yaml; print('Core packages: OK')"
if [ $? -ne 0 ]; then
    echo "Error: Failed to import required Python packages."
    exit 1
fi
echo "Python package verification successful."

# --- Script Execution ---
LONGITUDINAL_SCRIPT_PATH="${CODE_DIR}/src/longitudinal/longitudinal_main.py"

echo "DEBUG: Longitudinal script path: $LONGITUDINAL_SCRIPT_PATH"
echo "DEBUG: Script exists: $(ls -la "$LONGITUDINAL_SCRIPT_PATH" 2>/dev/null || echo "NOT FOUND")"

if [ ! -f "$LONGITUDINAL_SCRIPT_PATH" ]; then
    echo "Error: Longitudinal analysis Python script not found at $LONGITUDINAL_SCRIPT_PATH"
    echo "Note: This script now uses longitudinal_main.py (v3.0 unified architecture)"
    exit 1
fi

echo "Running unified longitudinal analysis with YAML configuration (v3.0)..."
echo "Command: python $LONGITUDINAL_SCRIPT_PATH --config $CONFIG_FILE $ADDITIONAL_FLAGS"

# Execute the Python script with YAML configuration
conda run -n mase_phi_hpc python "$LONGITUDINAL_SCRIPT_PATH" --config "$CONFIG_FILE" $ADDITIONAL_FLAGS

SCRIPT_EXIT_CODE=$?
if [ $SCRIPT_EXIT_CODE -eq 0 ]; then
    echo "Unified longitudinal analysis completed successfully for patient ${PATIENT_ID}."
else
    echo "Error: Unified longitudinal analysis failed for patient ${PATIENT_ID} with exit code $SCRIPT_EXIT_CODE."
    exit $SCRIPT_EXIT_CODE 
fi

# --- Output Summary ---
echo "=== UNIFIED LONGITUDINAL ANALYSIS COMPLETED (v3.0) ==="
echo "Patient: ${PATIENT_ID}"
echo "Configuration: ${CONFIG_FILE}"
echo "Results directory: ${OUTPUT_BASE}"
echo "Logs directory: ${LOG_DIR}"
echo ""
echo "Key output files will be in subdirectories of: ${OUTPUT_BASE}"
echo "  - longitudinal_pipeline/: Analysis results"
echo "  - logs/: Detailed execution logs"
echo "  - Comprehensive JSON tracking in results directory"
echo ""
echo "Primary SLURM job log is in the submission directory (slurm-$SLURM_JOB_ID.out)."
echo "--- Unified Longitudinal Analysis Script v3.0 End ---" 