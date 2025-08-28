#!/bin/bash
# Multi-Patient Pipeline Script for TracerX Marker Selection Pipeline
# Processes multiple SSM files using the existing single-patient run_pipeline.sh
# Usage: bash multi_patient_pipeline.sh <ssm_directory> <config_template> <output_base_directory> [--delay=N]
# 
# Example: bash multi_patient_pipeline.sh data/patients/ configs/template_multi_patient.yaml /path/to/results/
# Example: bash multi_patient_pipeline.sh data/patients/ configs/template_multi_patient.yaml /path/to/results/ --delay=30

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
    echo "  --help               Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 data/patients/ configs/template_multi_patient.yaml /path/to/results/"
    echo "  $0 data/patients/ configs/template_multi_patient.yaml /path/to/results/ --delay=60"
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

while [ "$#" -gt 0 ]; do
    case $1 in
        --delay=*)
            DELAY_SECONDS="${1#*=}"
            if ! [[ "$DELAY_SECONDS" =~ ^[0-9]+$ ]]; then
                echo "Error: --delay must be a positive integer"
                exit 1
            fi
            ;;
        --help)
            print_usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option '$1'. Only '--delay=N' and '--help' are supported."
            print_usage
            exit 1
            ;;
    esac
    shift
done

# --- Helper Functions ---
make_absolute_path() {
    local path="$1"
    if [[ ! "$path" = /* ]]; then
        echo "${SCRIPT_DIR}/${path}"
    else
        echo "$path"
    fi
}

setup_output_structure() {
    mkdir -p "${MULTI_CONFIGS_DIR}" "${MULTI_PATIENTS_DIR}" "${MULTI_LOGS_DIR}"
}

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
SSM_DIRECTORY=$(make_absolute_path "$SSM_DIRECTORY")
CONFIG_TEMPLATE=$(make_absolute_path "$CONFIG_TEMPLATE")
OUTPUT_BASE_DIRECTORY=$(make_absolute_path "$OUTPUT_BASE_DIRECTORY")

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
MASTER_LOG="${MULTI_LOGS_DIR}/multi_patient_master.log"

setup_output_structure

echo "=== Multi-Patient Pipeline Start: $(date) ==="
echo "SSM Directory: ${SSM_DIRECTORY}"
echo "Config Template: ${CONFIG_TEMPLATE}"
echo "Output Base Directory: ${OUTPUT_BASE_DIRECTORY}"
echo "Delay Between Submissions: ${DELAY_SECONDS} seconds"
echo "Master Pipeline Script: ${MASTER_PIPELINE}"
echo "Generated Configs Directory: ${MULTI_CONFIGS_DIR}"
echo "Patients Results Directory: ${MULTI_PATIENTS_DIR}"
echo "----------------------------------------"

# Setup logging after initial output
exec > >(tee -a "${MASTER_LOG}") 2>&1
echo "Master log: ${MASTER_LOG}"

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
    
    # Read template and substitute placeholders
    PROJECT_ROOT_DIR=$(dirname "${SCRIPT_DIR}")
    sed -e "s|PLACEHOLDER_PATIENT_ID|${patient_id}|g" \
        -e "s|PLACEHOLDER_SSM_FILE|${ssm_file}|g" \
        -e "s|PLACEHOLDER_OUTPUT_DIR|${patient_output_dir}|g" \
        -e "s|PLACEHOLDER_CODE_DIR|${PROJECT_ROOT_DIR}|g" \
        "$CONFIG_TEMPLATE" > "$config_output_path"
}

# --- Function to Submit Patient Pipeline ---
submit_patient_pipeline() {
    local patient_id="$1"
    local config_path="$2"
    local patient_log_file="$3"
    
    echo "Submitting pipeline for patient: $patient_id"
    echo "  Config: $config_path"
    echo "  Log: $patient_log_file"
    
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
    
    # For the config template, use the parent directory so pipeline can append /{patient_id}/initial
    PATIENT_BASE_DIR="${MULTI_PATIENTS_DIR}"
    
    echo ""
    echo "Processing patient: $PATIENT_ID"
    echo "----------------------------------------"
    
    # Create patient output directory and generate config
    mkdir -p "$PATIENT_OUTPUT_DIR"
    generate_patient_config "$PATIENT_ID" "$ssm_file" "$PATIENT_BASE_DIR" "$PATIENT_CONFIG"
    
    # Submit patient pipeline
    if submit_patient_pipeline "$PATIENT_ID" "$PATIENT_CONFIG" "$PATIENT_LOG"; then
        PATIENTS_SUBMITTED=$((PATIENTS_SUBMITTED + 1))
    else
        PATIENTS_FAILED=$((PATIENTS_FAILED + 1))
    fi
    
    # Apply delay if specified (except for last patient)
    if [ $DELAY_SECONDS -gt 0 ] && [ "$ssm_file" != "${SSM_FILES[-1]}" ]; then
        echo "Waiting $DELAY_SECONDS seconds before next submission..."
        sleep $DELAY_SECONDS
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

echo ""
echo "Monitor all jobs with: squeue -u $USER"
echo "Individual logs in: $MULTI_LOGS_DIR"
echo "Patient results in: $MULTI_PATIENTS_DIR"

echo "=== Multi-Patient Pipeline End: $(date) ==="

# Exit with error if any patients failed
if [ $PATIENTS_FAILED -gt 0 ]; then
    exit 1
fi 