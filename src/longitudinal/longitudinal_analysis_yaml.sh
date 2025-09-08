#!/bin/bash
# YAML Parser and SLURM Submission Script for Unified Longitudinal Analysis
# This script parses YAML configuration and submits job to SLURM
# Usage: bash longitudinal_analysis_yaml.sh <config_yaml_file> [additional_flags]

set -e

# --- Argument Parsing and Validation ---
if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Error: Incorrect number of arguments."
    echo "Usage: $0 <config_yaml_file> [additional_flags]"
    echo ""
    echo "Examples:"
    echo "  bash $0 configs/templates/longitudinal_minimal.yaml"
    echo "  bash $0 configs/templates/longitudinal_minimal.yaml '--debug'"
    echo ""
    echo "Configuration files should be in YAML format with all required parameters."
    echo "Uses unified longitudinal pipeline (v3.0)"
    exit 1
fi

CONFIG_FILE=$1
ADDITIONAL_FLAGS=${2:-""} # Optional additional flags like --debug, --no-plots

# --- Validate Configuration File ---
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file '$CONFIG_FILE' not found."
    exit 1
fi

# Convert to absolute path
if [[ ! "$CONFIG_FILE" = /* ]]; then
    CONFIG_FILE="$(pwd)/${CONFIG_FILE}"
    echo "Converted config file path to absolute: ${CONFIG_FILE}"
fi

echo "=== Longitudinal Analysis YAML Parser ==="
echo "Configuration file: $CONFIG_FILE"
echo "Additional flags: $ADDITIONAL_FLAGS"

# --- Extract Parameters from YAML ---
echo "Parsing YAML configuration..."

# Extract basic parameters
PATIENT_ID=$(grep "^patient_id:" "$CONFIG_FILE" | sed 's/patient_id: *"\?\([^"]*\)"\?/\1/' | tr -d '"' | xargs)
OUTPUT_BASE=$(grep "base_dir:" "$CONFIG_FILE" | sed 's/.*base_dir: *"\?\([^"]*\)"\?/\1/' | tr -d '"' | xargs)
CODE_DIR=$(grep "code_dir:" "$CONFIG_FILE" | sed 's/.*code_dir: *"\?\([^"]*\)"\?/\1/' | tr -d '"' | xargs)

# Extract HPC configuration
PARTITION=$(grep -A10 "^hpc:" "$CONFIG_FILE" | grep "partition:" | head -1 | sed 's/.*partition: *"\?\([^"]*\)"\?/\1/' | tr -d '"' | xargs)
CPUS=$(grep -A10 "^hpc:" "$CONFIG_FILE" | grep "cpus_per_task:" | head -1 | sed 's/.*cpus_per_task: *"\?\([^"]*\)"\?/\1/' | tr -d '"' | xargs)
MEMORY=$(grep -A10 "^hpc:" "$CONFIG_FILE" | grep "memory:" | head -1 | sed 's/.*memory: *"\?\([^"]*\)"\?/\1/' | tr -d '"' | xargs)
JOB_NAME=$(grep -A10 "^hpc:" "$CONFIG_FILE" | grep "job_name:" | head -1 | sed 's/.*job_name: *"\?\([^"]*\)"\?/\1/' | tr -d '"' | xargs)

# Set defaults if not found
PARTITION=${PARTITION:-"pool1"}
CPUS=${CPUS:-"2"}
MEMORY=${MEMORY:-"16G"}
JOB_NAME=${JOB_NAME:-"longitudinal_analysis"}

# Expand tilde paths
OUTPUT_BASE="${OUTPUT_BASE/#\~/$HOME}"
CODE_DIR="${CODE_DIR/#\~/$HOME}"

# --- Validation ---
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

if [ ! -d "$CODE_DIR" ]; then
    echo "Error: Code directory '$CODE_DIR' (from config) not found."
    exit 1
fi

# --- Display Configuration ---
echo "Parsed Configuration:"
echo "  Patient ID: $PATIENT_ID"
echo "  Code directory: $CODE_DIR"
echo "  Output base directory: $OUTPUT_BASE"
echo "  HPC Settings:"
echo "    Partition: $PARTITION"
echo "    CPUs: $CPUS"
echo "    Memory: $MEMORY"
echo "    Job name: ${PATIENT_ID}_${JOB_NAME}"

# --- Locate SLURM Execution Script ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
SLURM_SCRIPT="${SCRIPT_DIR}/longitudinal_slurm.sh"

if [ ! -f "$SLURM_SCRIPT" ]; then
    echo "Error: SLURM execution script not found at $SLURM_SCRIPT"
    exit 1
fi

# --- Setup Output Directory and Logs ---
mkdir -p "${OUTPUT_BASE}"
LOG_DIR="${OUTPUT_BASE}/logs"
mkdir -p "${LOG_DIR}"

# --- Submit to SLURM ---
echo "Submitting job to SLURM..."

JOB_ID=$(sbatch --parsable \
    --job-name="${PATIENT_ID}_${JOB_NAME}" \
    --partition="$PARTITION" \
    --cpus-per-task="$CPUS" \
    --mem="$MEMORY" \
    --time="02:00:00" \
    --output="${LOG_DIR}/longitudinal_${PATIENT_ID}_%j.out" \
    --error="${LOG_DIR}/longitudinal_${PATIENT_ID}_%j.err" \
    "$SLURM_SCRIPT" "$CONFIG_FILE" "$ADDITIONAL_FLAGS")

if [ $? -eq 0 ]; then
    echo "Job submitted successfully!"
    echo "  Job ID: $JOB_ID"
    echo "  Patient: $PATIENT_ID"
    echo "  Log files will be in: $LOG_DIR"
    echo "  Monitor with: squeue -j $JOB_ID"
    echo "  View output: tail -f ${LOG_DIR}/longitudinal_${PATIENT_ID}_${JOB_ID}.out"
else
    echo "Error: Failed to submit job to SLURM"
    exit 1
fi

echo "=== Longitudinal Analysis Submission Complete ==="