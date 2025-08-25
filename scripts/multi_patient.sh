#!/bin/bash
# Multi-Patient Pipeline Script for TracerX Marker Selection Pipeline
# Processes multiple SSM files using the existing single-patient run_pipeline.sh
# Usage: bash multi_patient_pipeline.sh <ssm_directory> <config_template> <output_base_directory> [--delay=N] [--dry-run]
# 
# Example: bash multi_patient_pipeline.sh data/patients/ configs/template_multi_patient.yaml /path/to/results/
# Example: bash multi_patient_pipeline.sh data/patients/ configs/template_multi_patient.yaml /path/to/results/ --delay=30 --dry-run

set -e  # Exit on any error

# --- Function to print usage ---
print_usage() {
    echo "Usage: $0 <ssm_directory> <config_template> <output_base_directory> [options]"
    echo ""
    echo "Arguments:"
    echo "  ssm_directory         Directory containing SSM files (*.txt, *.ssm)"
    echo "  config_template       YAML template file for patient configurations"
    echo "  output_base_directory Base directory where patient folders will be created"
    echo ""
    echo "Options:"
    echo "  --delay=N            Delay N seconds between patient submissions (default: 0)"
    echo "  --dry-run            Test configuration generation without submitting jobs"
    echo "  --help               Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 data/patients/ configs/template_multi_patient.yaml /path/to/results/"
    echo "  $0 data/patients/ configs/template_multi_patient.yaml /path/to/results/ --delay=60"
    echo "  $0 data/patients/ configs/template_multi_patient.yaml /path/to/results/ --dry-run"
}

# --- Input Validation ---
if [ "$#" -lt 3 ]; then
    echo "Error: Insufficient arguments."
    print_usage
    exit 1
fi

# Parse required arguments
SSM_DIRECTORY=$1
CONFIG_TEMPLATE=$2
OUTPUT_BASE_DIRECTORY=$3
shift 3

# Parse optional arguments
DELAY_SECONDS=0
DRY_RUN=""

while [ "$#" -gt 0 ]; do
    case $1 in
        --delay=*)
            DELAY_SECONDS="${1#*=}"
            if ! [[ "$DELAY_SECONDS" =~ ^[0-9]+$ ]]; then
                echo "Error: --delay must be a positive integer"
                exit 1
            fi
            ;;
        --dry-run)
            DRY_RUN="--dry-run"
            ;;
        --help)
            print_usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option '$1'"
            print_usage
            exit 1
            ;;
    esac
    shift
done

# --- Get script directory for absolute paths ---
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
echo "Multi-patient pipeline script directory: ${SCRIPT_DIR}"

# --- Validate Inputs ---
if [ ! -d "$SSM_DIRECTORY" ]; then
    echo "Error: SSM directory not found: $SSM_DIRECTORY"
    exit 1
fi

if [ ! -f "$CONFIG_TEMPLATE" ]; then
    echo "Error: Config template file not found: $CONFIG_TEMPLATE"
    exit 1
fi

# Convert to absolute paths
if [[ ! "$SSM_DIRECTORY" = /* ]]; then
    SSM_DIRECTORY="${SCRIPT_DIR}/${SSM_DIRECTORY}"
fi

if [[ ! "$CONFIG_TEMPLATE" = /* ]]; then
    CONFIG_TEMPLATE="${SCRIPT_DIR}/${CONFIG_TEMPLATE}"
fi

if [[ ! "$OUTPUT_BASE_DIRECTORY" = /* ]]; then
    OUTPUT_BASE_DIRECTORY="${SCRIPT_DIR}/${OUTPUT_BASE_DIRECTORY}"
fi

# Validate run_pipeline.sh exists
MAIN_PIPELINE="${SCRIPT_DIR}/run_pipeline.sh"
if [ ! -f "$MAIN_PIPELINE" ]; then
    echo "Error: run_pipeline.sh not found at: $MAIN_PIPELINE"
    exit 1
fi

# --- Create Output Structure ---
MULTI_CONFIGS_DIR="${OUTPUT_BASE_DIRECTORY}/configs/generated"
MULTI_PATIENTS_DIR="${OUTPUT_BASE_DIRECTORY}/patients"
MULTI_LOGS_DIR="${OUTPUT_BASE_DIRECTORY}/logs"

if [ "$DRY_RUN" != "--dry-run" ]; then
    mkdir -p "${MULTI_CONFIGS_DIR}"
    mkdir -p "${MULTI_PATIENTS_DIR}"
    mkdir -p "${MULTI_LOGS_DIR}"
fi

# --- Setup Master Log ---
MASTER_LOG="${MULTI_LOGS_DIR}/multi_patient_master.log"
if [ "$DRY_RUN" != "--dry-run" ]; then
    exec > >(tee -a "${MASTER_LOG}") 2>&1
fi

echo "=== Multi-Patient Pipeline Start: $(date) ==="
echo "SSM Directory: ${SSM_DIRECTORY}"
echo "Config Template: ${CONFIG_TEMPLATE}"
echo "Output Base Directory: ${OUTPUT_BASE_DIRECTORY}"
echo "Delay Between Submissions: ${DELAY_SECONDS} seconds"
echo "Dry Run Mode: ${DRY_RUN:-false}"
echo "Master Pipeline Script: ${MASTER_PIPELINE}"
echo "Generated Configs Directory: ${MULTI_CONFIGS_DIR}"
echo "Patients Results Directory: ${MULTI_PATIENTS_DIR}"
echo "----------------------------------------"

# --- Find SSM Files ---
echo "Scanning for SSM files in: ${SSM_DIRECTORY}"
SSM_FILES=($(find "${SSM_DIRECTORY}" -maxdepth 1 -name "*.txt" -o -name "*.ssm" | sort))

if [ ${#SSM_FILES[@]} -eq 0 ]; then
    echo "Error: No SSM files (*.txt, *.ssm) found in ${SSM_DIRECTORY}"
    exit 1
fi

echo "Found ${#SSM_FILES[@]} SSM files:"
for ssm_file in "${SSM_FILES[@]}"; do
    echo "  - $(basename "$ssm_file")"
done
echo "----------------------------------------"

# --- Function to Extract Patient ID from filename ---
extract_patient_id() {
    local ssm_file="$1"
    local filename=$(basename "$ssm_file")
    
    # Remove file extension
    local patient_id="${filename%.*}"
    
    # Clean up common prefixes/suffixes if needed
    patient_id="${patient_id%_ssm}"
    patient_id="${patient_id%_subset}"
    
    echo "$patient_id"
}

# --- Function to Generate Patient Config ---
generate_patient_config() {
    local patient_id="$1"
    local ssm_file="$2"
    local patient_output_dir="$3"
    local config_output_path="$4"
    
    echo "Generating config for patient: $patient_id"
    echo "  SSM file: $ssm_file"
    echo "  Output directory: $patient_output_dir"
    echo "  Config path: $config_output_path"
    
    # Read template and substitute placeholders
    # Use project root directory (parent of scripts directory) for CODE_DIR
    PROJECT_ROOT_DIR=$(dirname "${SCRIPT_DIR}")
    sed -e "s|PLACEHOLDER_PATIENT_ID|${patient_id}|g" \
        -e "s|PLACEHOLDER_SSM_FILE|${ssm_file}|g" \
        -e "s|PLACEHOLDER_OUTPUT_DIR|${patient_output_dir}|g" \
        -e "s|PLACEHOLDER_CODE_DIR|${PROJECT_ROOT_DIR}|g" \
        "$CONFIG_TEMPLATE" > "$config_output_path"
    
    echo "  Config generated successfully"
}

# --- Function to Submit Patient Pipeline ---
submit_patient_pipeline() {
    local patient_id="$1"
    local config_path="$2"
    local patient_log_file="$3"
    
    echo "Submitting pipeline for patient: $patient_id"
    echo "  Config: $config_path"
    echo "  Log: $patient_log_file"
    
    if [ "$DRY_RUN" == "--dry-run" ]; then
        echo "  DRY RUN: Would execute: bash ${MAIN_PIPELINE} ${config_path}"
        return 0
    fi
    
    # Submit the master pipeline and capture output
    {
        echo "=== Patient $patient_id Pipeline Submission: $(date) ==="
        echo "Config: $config_path"
        echo "Command: bash ${MAIN_PIPELINE} ${config_path}"
        echo "----------------------------------------"
        
        bash "${MAIN_PIPELINE}" "${config_path}"
        local exit_code=$?
        
        echo "----------------------------------------"
        echo "Pipeline submission completed with exit code: $exit_code"
        echo "=== Patient $patient_id Submission End: $(date) ==="
        
        return $exit_code
        
    } > "$patient_log_file" 2>&1
    
    local submission_status=$?
    if [ $submission_status -eq 0 ]; then
        echo "  ✓ Pipeline submitted successfully"
    else
        echo "  ✗ Pipeline submission failed (exit code: $submission_status)"
        echo "  Check log: $patient_log_file"
    fi
    
    return $submission_status
}

# --- Process Each Patient ---
PATIENTS_SUBMITTED=0
PATIENTS_FAILED=0

for ssm_file in "${SSM_FILES[@]}"; do
    # Extract patient ID from filename
    PATIENT_ID=$(extract_patient_id "$ssm_file")
    
    # Define paths for this patient
    PATIENT_OUTPUT_DIR="${MULTI_PATIENTS_DIR}/${PATIENT_ID}"
    PATIENT_CONFIG="${MULTI_CONFIGS_DIR}/${PATIENT_ID}_config.yaml"
    PATIENT_LOG="${MULTI_LOGS_DIR}/${PATIENT_ID}_submission.log"
    
    echo ""
    echo "Processing patient: $PATIENT_ID"
    echo "----------------------------------------"
    
    # Create patient output directory
    if [ "$DRY_RUN" != "--dry-run" ]; then
        mkdir -p "$PATIENT_OUTPUT_DIR"
    fi
    
    # Generate patient-specific config
    if [ "$DRY_RUN" != "--dry-run" ]; then
        generate_patient_config "$PATIENT_ID" "$ssm_file" "$PATIENT_OUTPUT_DIR" "$PATIENT_CONFIG"
    else
        echo "DRY RUN: Would generate config for patient $PATIENT_ID"
        echo "  Template: $CONFIG_TEMPLATE"
        echo "  Output config: $PATIENT_CONFIG"
    fi
    
    # Submit patient pipeline
    if submit_patient_pipeline "$PATIENT_ID" "$PATIENT_CONFIG" "$PATIENT_LOG"; then
        PATIENTS_SUBMITTED=$((PATIENTS_SUBMITTED + 1))
    else
        PATIENTS_FAILED=$((PATIENTS_FAILED + 1))
    fi
    
    # Apply delay if specified (except for last patient)
    if [ $DELAY_SECONDS -gt 0 ] && [ "$ssm_file" != "${SSM_FILES[-1]}" ]; then
        echo "Waiting $DELAY_SECONDS seconds before next submission..."
        if [ "$DRY_RUN" != "--dry-run" ]; then
            sleep $DELAY_SECONDS
        fi
    fi
done

# --- Final Summary ---
echo ""
echo "========================================"
echo "Multi-Patient Pipeline Summary"
echo "========================================"
echo "Total SSM files found: ${#SSM_FILES[@]}"
echo "Patients successfully submitted: $PATIENTS_SUBMITTED"
echo "Patients failed: $PATIENTS_FAILED"
echo "Output directory: $OUTPUT_BASE_DIRECTORY"
echo "Master log: $MASTER_LOG"

if [ "$DRY_RUN" == "--dry-run" ]; then
    echo ""
    echo "DRY RUN COMPLETED - No jobs were actually submitted"
    echo "All configurations would be generated successfully"
    echo "Use without --dry-run to submit actual jobs"
else
    echo ""
    echo "Monitor all jobs with: squeue -u $USER"
    echo "Individual logs in: $MULTI_LOGS_DIR"
    echo "Patient results in: $MULTI_PATIENTS_DIR"
fi

echo "=== Multi-Patient Pipeline End: $(date) ==="

# Exit with error if any patients failed
if [ $PATIENTS_FAILED -gt 0 ]; then
    exit 1
fi 