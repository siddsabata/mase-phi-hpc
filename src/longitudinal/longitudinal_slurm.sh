#!/bin/bash
# SLURM execution script for unified longitudinal analysis
# This script runs inside SLURM with proper #SBATCH directives
# Called by longitudinal_analysis_yaml.sh

#SBATCH --job-name=longitudinal_analysis
#SBATCH --partition=pool1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=longitudinal_analysis_%j.out
#SBATCH --error=longitudinal_analysis_%j.err

set -e

# Arguments passed from submission script
CONFIG_FILE=$1
ADDITIONAL_FLAGS=${2:-""}

echo "=== SLURM Longitudinal Analysis Execution ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Configuration: $CONFIG_FILE"
echo "Additional flags: $ADDITIONAL_FLAGS"
echo "Node: $SLURMD_NODENAME"
echo "================================="

# Extract basic parameters for logging
PATIENT_ID=$(grep "^patient_id:" "$CONFIG_FILE" | sed 's/patient_id: *"\?\([^"]*\)"\?/\1/' | tr -d '"')
CODE_DIR=$(grep "code_dir:" "$CONFIG_FILE" | sed 's/.*code_dir: *"\?\([^"]*\)"\?/\1/' | tr -d '"')

# Expand tilde paths
CODE_DIR="${CODE_DIR/#\~/$HOME}"

echo "Patient ID: $PATIENT_ID"
echo "Code directory: $CODE_DIR"

# Validate code directory exists
if [ ! -d "$CODE_DIR" ]; then
    echo "Error: Code directory '$CODE_DIR' not found."
    exit 1
fi

# Load Gurobi module
echo "Loading Gurobi 11.0.2 module..."
module load gurobi1102
if [ $? -ne 0 ]; then
    echo "Error: Failed to load gurobi1102 module. Exiting."
    exit 1
fi
echo "Gurobi module loaded successfully."

# Change to code directory
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

# Execute the longitudinal analysis
LONGITUDINAL_SCRIPT_PATH="${CODE_DIR}/src/longitudinal/longitudinal_main.py"

echo "DEBUG: Longitudinal script path: $LONGITUDINAL_SCRIPT_PATH"
echo "DEBUG: Script exists: $(ls -la "$LONGITUDINAL_SCRIPT_PATH" 2>/dev/null || echo "NOT FOUND")"

if [ ! -f "$LONGITUDINAL_SCRIPT_PATH" ]; then
    echo "Error: Longitudinal analysis Python script not found at $LONGITUDINAL_SCRIPT_PATH"
    exit 1
fi

echo "Running unified longitudinal analysis..."
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

echo "=== Longitudinal Analysis Job Complete ==="