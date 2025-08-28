#!/bin/bash
# Multi-sample marker selection using pre-filtered SSM data
# Usage: sbatch marker_selection.sh <patient_id> <aggregation_directory> <ssm_file_path> <code_directory> [read_depth]

set -e

# Load required modules from configuration
if [ -n "${MARKER_SELECTION_MODULES}" ]; then
    echo "Loading modules: ${MARKER_SELECTION_MODULES}"
    for module in ${MARKER_SELECTION_MODULES}; do
        echo "Loading module: $module"
        module load "$module"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to load module '$module'. Exiting."
            exit 1
        fi
        echo "Module $module loaded successfully."
    done
else
    echo "Warning: No modules specified in MARKER_SELECTION_MODULES. Proceeding without loading modules."
fi

# --- Argument Parsing and Validation ---
if [ "$#" -lt 4 ] || [ "$#" -gt 5 ]; then
  echo "Error: Incorrect number of arguments."
  echo "Usage: sbatch $0 <patient_id> <aggregation_directory> <ssm_file_path> <code_directory> [read_depth]"
  echo "Example: sbatch $0 CRUK0001 /path/to/data/CRUK0001/initial/aggregation_results /path/to/data/CRUK0001/initial/ssm_filtered.txt /path/to/tracerx-mp 1500"
  echo "Note: SSM file should be pre-filtered from bootstrap stage (VAF filtering applied there)"
  exit 1
fi

PATIENT_ID=$1
AGGREGATION_DIR=$2
SSM_FILE=$3
CODE_DIR=$4
READ_DEPTH=${5:-1500} # Default to 1500 if not provided

if [ ! -d "$AGGREGATION_DIR" ]; then
  echo "Error: Aggregation directory '$AGGREGATION_DIR' not found."
  exit 1
fi

if [ ! -f "$SSM_FILE" ]; then
  echo "Error: SSM file '$SSM_FILE' not found."
  exit 1
fi

if [ ! -d "$CODE_DIR" ]; then
  echo "Error: Code directory '$CODE_DIR' not found."
  exit 1
fi

echo "--- Marker Selection Script Execution ---"

# --- Setup directories and paths ---
# Get the parent directory of the aggregation directory (should be the 'initial' directory)
INITIAL_DIR="$(dirname "$AGGREGATION_DIR")"
MARKERS_DIR="${INITIAL_DIR}/markers"

# Create the markers directory if it doesn't exist
mkdir -p "${MARKERS_DIR}"
echo "Job ID: $SLURM_JOB_ID"
echo "Patient ID: ${PATIENT_ID}"
echo "Aggregation Directory: ${AGGREGATION_DIR}"
echo "Markers Directory: ${MARKERS_DIR}"
echo "SSM File (pre-filtered): ${SSM_FILE}"
echo "Code Directory: ${CODE_DIR}"
echo "Read Depth: ${READ_DEPTH}"
echo "Note: VAF filtering (threshold=0.9, any_high) applied in bootstrap stage"
echo "---------------------------------------"

# --- Environment Setup ---
echo "Setting up conda environment..."
cd "$CODE_DIR"
# No need to sync conda environment as it should already exist

# Verify Gurobi is accessible from Python
echo "Verifying Gurobi is accessible from Python..."
conda run -n mase_phi_hpc python -c "import gurobipy; print(f'Gurobi version: {gurobipy.gurobi.version()}')"
if [ $? -ne 0 ]; then
  echo "Error: Failed to import gurobipy or access Gurobi. Check that the module is properly loaded and gurobipy is installed."
  exit 1
fi
echo "Gurobi verification successful."

# --- Script Paths and Execution ---
# Use the absolute path to the run_data_multi_sample.py script based on the provided code directory
MARKER_SCRIPT_PATH="${CODE_DIR}/src/markers/run_data_multi_sample.py"

# Debug: Print the actual script path being used
echo "DEBUG: Marker script path: $MARKER_SCRIPT_PATH"
echo "DEBUG: Script exists check: $(ls -la "$MARKER_SCRIPT_PATH" 2>/dev/null || echo "NOT FOUND")"

if [ ! -f "$MARKER_SCRIPT_PATH" ]; then
  echo "Error: Multi-sample marker selection Python script not found at $MARKER_SCRIPT_PATH. Exiting."
  exit 1
fi

echo "Running multi-sample Python marker selection script: $MARKER_SCRIPT_PATH"

# Note: Filtering parameters removed - SSM file is already pre-filtered from bootstrap stage

conda run -n mase_phi_hpc python "$MARKER_SCRIPT_PATH" "${PATIENT_ID}" \
  -a "${AGGREGATION_DIR}" \
  -s "${SSM_FILE}" \
  -r "${READ_DEPTH}" \
  -o "${MARKERS_DIR}"

SCRIPT_EXIT_CODE=$?
if [ $SCRIPT_EXIT_CODE -eq 0 ]; then
  echo "Marker selection completed successfully for patient ${PATIENT_ID}."
else
  echo "Error: Marker selection Python script failed for patient ${PATIENT_ID} with exit code $SCRIPT_EXIT_CODE."
  # Consider exiting the sbatch script with the Python script's error code
  exit $SCRIPT_EXIT_CODE
fi

echo "Detailed Python script output is in: ${MARKERS_DIR}"
echo "Primary SLURM job log is in the submission directory (slurm-$SLURM_JOB_ID.out)."
echo "--- Marker Selection Script End ---"

