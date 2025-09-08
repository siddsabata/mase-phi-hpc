#!/bin/bash
#SBATCH --job-name=longitudinal_simple
#SBATCH --partition=pool1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=longitudinal_%j.out
#SBATCH --error=longitudinal_%j.err

# Simple SLURM script for longitudinal analysis
# Usage: sbatch longitudinal_simple.sh <patient_id> <aggregation_dir> <ssm_file> <longitudinal_data> <output_dir>
# Example: sbatch longitudinal_simple.sh CRUK0044 /path/to/aggregation_results /path/to/ssm.txt /path/to/liquid.csv /path/to/output

set -e

# Arguments
PATIENT_ID=$1
AGGREGATION_DIR=$2
SSM_FILE=$3
LONGITUDINAL_DATA=$4
OUTPUT_DIR=$5

# Validate arguments
if [ $# -ne 5 ]; then
    echo "Error: Expected 5 arguments"
    echo "Usage: sbatch $0 <patient_id> <aggregation_dir> <ssm_file> <longitudinal_data> <output_dir>"
    exit 1
fi

# Validate input files exist
if [ ! -d "$AGGREGATION_DIR" ]; then
    echo "Error: Aggregation directory not found: $AGGREGATION_DIR"
    exit 1
fi

if [ ! -f "$SSM_FILE" ]; then
    echo "Error: SSM file not found: $SSM_FILE"
    exit 1
fi

if [ ! -f "$LONGITUDINAL_DATA" ]; then
    echo "Error: Longitudinal data file not found: $LONGITUDINAL_DATA"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Get script directory (where this script is located)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
CODE_DIR="$(cd "$SCRIPT_DIR/../.." &>/dev/null && pwd)"

echo "=== Longitudinal Analysis Started ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Patient ID: $PATIENT_ID"
echo "Aggregation dir: $AGGREGATION_DIR"
echo "SSM file: $SSM_FILE"
echo "Longitudinal data: $LONGITUDINAL_DATA"
echo "Output dir: $OUTPUT_DIR"
echo "Code dir: $CODE_DIR"
echo "Node: $SLURMD_NODENAME"
echo "=================================="

# Load required modules
echo "Loading Gurobi module..."
module load gurobi1102

# Change to code directory
cd "$CODE_DIR"

# Run the longitudinal analysis with fixed parameters
echo "Running longitudinal analysis..."
conda run -n mase_phi_hpc python src/longitudinal/longitudinal_main.py \
    --patient_id "$PATIENT_ID" \
    --aggregation_dir "$AGGREGATION_DIR" \
    --ssm_file "$SSM_FILE" \
    --longitudinal_data "$LONGITUDINAL_DATA" \
    --output_dir "$OUTPUT_DIR" \
    --n_markers 5 \
    --read_depth 90000 \
    --method phylowgs \
    --track_clone_freq

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo "Longitudinal analysis completed successfully for patient $PATIENT_ID"
else
    echo "Error: Longitudinal analysis failed with exit code $EXIT_CODE"
    exit $EXIT_CODE
fi

echo "=== Longitudinal Analysis Complete ==="