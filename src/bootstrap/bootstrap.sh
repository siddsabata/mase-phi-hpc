#!/bin/bash
# Generates bootstrap samples from SSM data with VAF filtering
# Usage: sbatch bootstrap.sh <input_ssm_file> <output_directory> <code_directory> [num_bootstraps]

# --- Argument Parsing ---
if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
  echo "Error: Incorrect number of arguments."
  echo "Usage: sbatch $0 <input_ssm_file> <output_directory> <code_directory> [num_bootstraps (default: 100)]"
  exit 1
fi

INPUT_SSM_FILE=$1
OUTPUT_BASE_DIR=$2
CODE_DIR=$3
NUM_BOOTSTRAPS=${4:-100} # Default to 100 if not provided

echo "--- Bootstrap Script Start (initial output to SLURM default log) ---"

# --- Setup Log and Output Directories ---
mkdir -p "$OUTPUT_BASE_DIR/logs"

# Define the actual directory where bootstrapN folders will be created
BOOTSTRAP_DATA_DIR="$OUTPUT_BASE_DIR/bootstraps"
mkdir -p "$BOOTSTRAP_DATA_DIR"

echo "--- Bootstrap Script Execution ---"
echo "Input SSM File: $INPUT_SSM_FILE"
echo "Output Base Directory: $OUTPUT_BASE_DIR"
echo "Code Directory: $CODE_DIR"
echo "Number of Bootstraps: $NUM_BOOTSTRAPS"

# Verify the input file exists
if [ ! -f "$INPUT_SSM_FILE" ]; then
  echo "Error: Input SSM file not found at $INPUT_SSM_FILE"
  exit 1
fi

# --- Environment Setup ---
echo "Using conda environment mase_phi_hpc..."
cd "$CODE_DIR"

# --- Run Bootstrap Script ---
echo "Running bootstrap.py..."

# Use the absolute path to the bootstrap.py file based on the provided code directory
BOOTSTRAP_PY_PATH="${CODE_DIR}/src/bootstrap/bootstrap.py"
echo "Using bootstrap.py at: $BOOTSTRAP_PY_PATH"

# Verify the bootstrap.py file exists
if [ ! -f "$BOOTSTRAP_PY_PATH" ]; then
  echo "Error: Could not find bootstrap.py at $BOOTSTRAP_PY_PATH"
  echo "Current working directory: $(pwd)"
  echo "Files in code directory: $(ls -la $CODE_DIR/src/bootstrap/)"
  exit 1
fi

conda run -n mase_phi_hpc python "$BOOTSTRAP_PY_PATH" \
  -i "$INPUT_SSM_FILE" \
  -o "$BOOTSTRAP_DATA_DIR" \
  -n "$NUM_BOOTSTRAPS"

SCRIPT_EXIT_CODE=$?
if [ $SCRIPT_EXIT_CODE -eq 0 ]; then
  echo "bootstrap.py completed successfully."
else
  echo "Error: bootstrap.py script failed with exit code $SCRIPT_EXIT_CODE."
  exit $SCRIPT_EXIT_CODE
fi

echo "Bootstrapped files should be in subdirectories within: $BOOTSTRAP_DATA_DIR"
echo "Primary SLURM job log is in the submission directory."
echo "--- Bootstrap Script End ---"

