#!/bin/bash
#SBATCH --job-name=aggregation
# SLURM will use default log files (e.g., slurm-%j.out in submission dir).
#SBATCH --partition=pool1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

set -e

# --- Argument Parsing and Validation ---
if [ "$#" -ne 4 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: sbatch $0 <patient_id> <bootstrap_parent_directory> <output_directory> <code_directory>"
    echo "Example: sbatch $0 CRUK0001 /path/to/data/initial/bootstraps /path/to/data/initial/aggregation_results /path/to/tracerx-mp"
    exit 1
fi

PATIENT_ID=$1
BOOTSTRAP_PARENT_DIR=$2 # This is the directory containing bootstrapN folders
OUTPUT_DIR=$3           # Explicit output directory for aggregation results
CODE_DIR=$4
NUM_BOOTSTRAPS=100      # Hardcoded as per previous request

# --- Validate Input Directories ---
if [ ! -d "$BOOTSTRAP_PARENT_DIR" ]; then
    echo "Error: Bootstrap parent directory '$BOOTSTRAP_PARENT_DIR' not found."
    exit 1
fi

if [ ! -d "$CODE_DIR" ]; then
    echo "Error: Code directory '$CODE_DIR' not found."
    exit 1
fi

echo "--- Aggregation Script Execution ---"
echo "Job ID: $SLURM_JOB_ID"
echo "Patient ID: ${PATIENT_ID}"
echo "Bootstrap Parent Directory: ${BOOTSTRAP_PARENT_DIR}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "Code Directory: ${CODE_DIR}"
echo "Number of Bootstraps (hardcoded): ${NUM_BOOTSTRAPS}"
echo "---------------------------------------"

# --- Environment Setup ---
echo "Using conda environment mase_phi_hpc..."
cd "$CODE_DIR"

# --- Script Paths and Execution ---
# Use the absolute path to the step3_aggregate.py script based on the provided code directory
PROCESS_SCRIPT_PATH="${CODE_DIR}/src/aggregation/step3_aggregate.py"

if [ ! -f "$PROCESS_SCRIPT_PATH" ]; then
    echo "Error: Aggregation Python script not found at $PROCESS_SCRIPT_PATH. Exiting."
    exit 1
fi

bootstrap_list=$(seq -s ' ' 1 $NUM_BOOTSTRAPS)

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "Running Python aggregation script: $PROCESS_SCRIPT_PATH"
echo "Input bootstrap data from: $BOOTSTRAP_PARENT_DIR"
echo "Output aggregation results to: $OUTPUT_DIR"

conda run -n mase_phi_hpc python "$PROCESS_SCRIPT_PATH" "${PATIENT_ID}" \
    --bootstrap-list $bootstrap_list \
    --bootstrap-parent-dir "${BOOTSTRAP_PARENT_DIR}" \
    --output-dir "${OUTPUT_DIR}"

SCRIPT_EXIT_CODE=$?
if [ $SCRIPT_EXIT_CODE -eq 0 ]; then
    echo "Aggregation Python script completed successfully for patient ${PATIENT_ID}."
    echo "Aggregation results are available in: ${OUTPUT_DIR}"
else
    echo "Error: Aggregation Python script failed for patient ${PATIENT_ID} with exit code $SCRIPT_EXIT_CODE."
    exit $SCRIPT_EXIT_CODE
fi

echo "Detailed script execution logs are in: $LOG_DIR_IN_PARENT/ (aggregation_execution.log/err)"
echo "Primary SLURM job log is in the submission directory."
echo "--- Aggregation Script End (redirected output) ---" 