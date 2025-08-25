#!/bin/bash
#
# Multi-sample marker selection SLURM script
# Uses run_data_multi_sample.py with flexible VAF filtering strategies
#
#SBATCH --job-name=marker_selection_multi
# SLURM will use default log files (e.g., slurm-%j.out in submission dir for initial output).
#SBATCH --partition=pool1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

set -e

# load gurobi
echo "Loading Gurobi 11.0.2 module..."
module load gurobi1102
if [ $? -ne 0 ]; then
    echo "Error: Failed to load gurobi1102 module. Exiting."
    exit 1
fi
echo "Gurobi module loaded successfully."

# --- Argument Parsing and Validation ---
if [ "$#" -lt 4 ] || [ "$#" -gt 8 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: sbatch $0 <patient_id> <aggregation_directory> <ssm_file_path> <code_directory> [read_depth] [filter_strategy] [filter_threshold] [filter_samples]"
    echo "Example: sbatch $0 CRUK0001 /path/to/data/CRUK0001/initial/aggregation_results /path/to/data/CRUK0001/ssm.txt /path/to/tracerx-mp 1500 any_high 0.9"
    echo "Filter strategies: any_high, all_high, majority_high, specific_samples"
    exit 1
fi

PATIENT_ID=$1
AGGREGATION_DIR=$2
SSM_FILE=$3
CODE_DIR=$4
READ_DEPTH=${5:-1500} # Default to 1500 if not provided
FILTER_STRATEGY=${6:-any_high} # Default to any_high if not provided
FILTER_THRESHOLD=${7:-0.9} # Default to 0.9 if not provided
FILTER_SAMPLES=${8:-""} # Optional specific samples for specific_samples strategy

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

echo "--- Marker Selection Script Start (initial output to SLURM default log) ---"

# --- Setup directories and paths ---
# Get the parent directory of the aggregation directory (should be the 'initial' directory)
INITIAL_DIR="$(dirname "$AGGREGATION_DIR")"
MARKERS_DIR="${INITIAL_DIR}/markers"
LOG_DIR="${INITIAL_DIR}/logs"

# Create the markers directory if it doesn't exist
mkdir -p "${MARKERS_DIR}"
mkdir -p "${LOG_DIR}"

# Redirect subsequent script output to files in LOG_DIR/
exec > "$LOG_DIR/marker_selection_execution.log" 2> "$LOG_DIR/marker_selection_execution.err"

# From this point, all echo and command output goes to the files defined above.
echo "--- Marker Selection Script Execution (output redirected) ---"
echo "Job ID: $SLURM_JOB_ID"
echo "Patient ID: ${PATIENT_ID}"
echo "Aggregation Directory: ${AGGREGATION_DIR}"
echo "Markers Directory: ${MARKERS_DIR}"
echo "SSM File: ${SSM_FILE}"
echo "Code Directory: ${CODE_DIR}"
echo "Read Depth: ${READ_DEPTH}"
echo "Filter Strategy: ${FILTER_STRATEGY}"
echo "Filter Threshold: ${FILTER_THRESHOLD}"
echo "Filter Samples: ${FILTER_SAMPLES}"
echo "Log Directory: ${LOG_DIR}"
echo "---------------------------------------"

# --- Environment Activation ---
# Assuming the same environment as aggregation.sh, or a suitable one for run_data.py
echo "Activating conda environment: markers_env" 
source ~/miniconda3/bin/activate markers_env
if [ $? -ne 0 ]; then
    echo "Error: Failed to activate conda environment 'markers_env'. Exiting."
    exit 1
fi
echo "Conda environment activated."

# Verify Gurobi is accessible from Python
echo "Verifying Gurobi is accessible from Python..."
python -c "import gurobipy; print(f'Gurobi version: {gurobipy.gurobi.version()}')"
if [ $? -ne 0 ]; then
    echo "Error: Failed to import gurobipy or access Gurobi. Check that the module is properly loaded and gurobipy is installed."
    exit 1
fi
echo "Gurobi verification successful."

# --- Script Paths and Execution ---
# Use the absolute path to the step4_run_data_multi_sample.py script based on the provided code directory
MARKER_SCRIPT_PATH="${CODE_DIR}/src/markers/step4_run_data_multi_sample.py"

# Debug: Print the actual script path being used
echo "DEBUG: Marker script path: $MARKER_SCRIPT_PATH"
echo "DEBUG: Script exists check: $(ls -la "$MARKER_SCRIPT_PATH" 2>/dev/null || echo "NOT FOUND")"

if [ ! -f "$MARKER_SCRIPT_PATH" ]; then
    echo "Error: Multi-sample marker selection Python script not found at $MARKER_SCRIPT_PATH. Exiting."
    exit 1
fi

echo "Running multi-sample Python marker selection script: $MARKER_SCRIPT_PATH"

# Build the command with required arguments using shorthand options
CMD="python \"$MARKER_SCRIPT_PATH\" \"${PATIENT_ID}\" \
    -a \"${AGGREGATION_DIR}\" \
    -s \"${SSM_FILE}\" \
    -r \"${READ_DEPTH}\" \
    -o \"${MARKERS_DIR}\" \
    -f \"${FILTER_STRATEGY}\" \
    -t \"${FILTER_THRESHOLD}\""

# Add filter samples if specified (for specific_samples strategy)
if [ -n "$FILTER_SAMPLES" ]; then
    CMD="$CMD --filter-samples $FILTER_SAMPLES"
fi

echo "Executing command: $CMD"
eval $CMD

SCRIPT_EXIT_CODE=$?
if [ $SCRIPT_EXIT_CODE -eq 0 ]; then
    echo "Marker selection completed successfully for patient ${PATIENT_ID}."
else
    echo "Error: Marker selection Python script failed for patient ${PATIENT_ID} with exit code $SCRIPT_EXIT_CODE."
    # Consider exiting the sbatch script with the Python script's error code
    exit $SCRIPT_EXIT_CODE 
fi

echo "Detailed Python script output is in: ${MARKERS_DIR}"
echo "Script execution logs are in: ${LOG_DIR} (marker_selection_execution.log/err)"
echo "Primary SLURM job log (initial messages) is in the submission directory (e.g., slurm-$SLURM_JOB_ID.out)."
echo "--- Marker Selection Script End (redirected output) ---" 