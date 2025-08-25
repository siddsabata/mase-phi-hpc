#!/bin/bash
# Validation script to test the new directory structure and path references
# This script verifies that all components are properly organized and accessible

set -e

echo "=== TracerX Pipeline Structure Validation ==="
echo "Validating directory structure and path references..."
echo

# Function to check if a file exists
check_file() {
    local file_path="$1"
    local description="$2"
    
    if [ -f "$file_path" ]; then
        echo "✅ $description: $file_path"
    else
        echo "❌ $description: $file_path (NOT FOUND)"
        return 1
    fi
}

# Function to check if a directory exists
check_directory() {
    local dir_path="$1"
    local description="$2"
    
    if [ -d "$dir_path" ]; then
        echo "✅ $description: $dir_path"
    else
        echo "❌ $description: $dir_path (NOT FOUND)"
        return 1
    fi
}

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"

echo "Root directory: $ROOT_DIR"
echo "Script directory: $SCRIPT_DIR"
echo

# Check main directories
echo "--- Checking Main Directory Structure ---"
check_directory "$ROOT_DIR/src" "Source code directory"
check_directory "$ROOT_DIR/scripts" "Scripts directory"
check_directory "$ROOT_DIR/configs" "Configuration directory"
check_directory "$ROOT_DIR/data" "Data directory"
check_directory "$ROOT_DIR/docs" "Documentation directory"
echo

# Check source subdirectories
echo "--- Checking Source Subdirectories ---"
check_directory "$ROOT_DIR/src/bootstrap" "Bootstrap stage"
check_directory "$ROOT_DIR/src/phylowgs" "PhyloWGS stage"
check_directory "$ROOT_DIR/src/aggregation" "Aggregation stage"
check_directory "$ROOT_DIR/src/markers" "Markers stage"
check_directory "$ROOT_DIR/src/longitudinal" "Longitudinal stage"
check_directory "$ROOT_DIR/src/common" "Common utilities"
echo

# Check main scripts
echo "--- Checking Main Scripts ---"
check_file "$ROOT_DIR/scripts/run_pipeline.sh" "Main pipeline script"
check_file "$ROOT_DIR/scripts/multi_patient.sh" "Multi-patient script"
check_file "$ROOT_DIR/scripts/install.sh" "Installation script"
check_file "$ROOT_DIR/scripts/validate_setup.sh" "Validation script"
echo

# Check stage scripts
echo "--- Checking Stage Scripts ---"
check_file "$ROOT_DIR/src/bootstrap/bootstrap.sh" "Bootstrap SLURM script"
check_file "$ROOT_DIR/src/phylowgs/phylowgs.sh" "PhyloWGS SLURM script"
check_file "$ROOT_DIR/src/aggregation/aggregation.sh" "Aggregation SLURM script"
check_file "$ROOT_DIR/src/markers/marker_selection.sh" "Marker selection SLURM script"
check_file "$ROOT_DIR/src/longitudinal/longitudinal_analysis_yaml.sh" "Longitudinal analysis SLURM script"
echo

# Check Python files
echo "--- Checking Key Python Files ---"
check_file "$ROOT_DIR/src/bootstrap/step1_bootstrap.py" "Bootstrap Python script"
check_file "$ROOT_DIR/src/aggregation/step3_aggregate.py" "Aggregation Python script"
check_file "$ROOT_DIR/src/markers/step4_run_data_multi_sample.py" "Markers Python script"
check_file "$ROOT_DIR/src/longitudinal/longitudinal_main.py" "Longitudinal main Python script"
echo

# Check configuration files
echo "--- Checking Configuration Files ---"
check_directory "$ROOT_DIR/configs/analysis" "Analysis configurations"
check_directory "$ROOT_DIR/configs/templates" "Template configurations"
echo

# Check data directories
echo "--- Checking Data Directories ---"
check_directory "$ROOT_DIR/data/input" "Input data directory"
check_directory "$ROOT_DIR/data/test" "Test data directory"
check_directory "$ROOT_DIR/data/raw" "Raw data directory"
echo

# Check documentation
echo "--- Checking Documentation ---"
check_file "$ROOT_DIR/README.md" "Main README"
check_directory "$ROOT_DIR/docs" "Documentation directory"
echo

# Test path resolution in scripts
echo "--- Testing Path Resolution ---"
echo "Testing path references in main pipeline script..."

# Check if the main pipeline script uses the correct paths
if grep -q "src/bootstrap/bootstrap.sh" "$ROOT_DIR/scripts/run_pipeline.sh"; then
    echo "✅ Bootstrap path correctly updated in main pipeline"
else
    echo "❌ Bootstrap path not found in main pipeline"
fi

if grep -q "src/phylowgs/phylowgs.sh" "$ROOT_DIR/scripts/run_pipeline.sh"; then
    echo "✅ PhyloWGS path correctly updated in main pipeline"
else
    echo "❌ PhyloWGS path not found in main pipeline"
fi

if grep -q "src/aggregation/aggregation.sh" "$ROOT_DIR/scripts/run_pipeline.sh"; then
    echo "✅ Aggregation path correctly updated in main pipeline"
else
    echo "❌ Aggregation path not found in main pipeline"
fi

if grep -q "src/markers/marker_selection.sh" "$ROOT_DIR/scripts/run_pipeline.sh"; then
    echo "✅ Markers path correctly updated in main pipeline"
else
    echo "❌ Markers path not found in main pipeline"
fi

echo

# Test multi-patient script
echo "Testing multi-patient script..."
if grep -q "run_pipeline.sh" "$ROOT_DIR/scripts/multi_patient.sh"; then
    echo "✅ Multi-patient script correctly references run_pipeline.sh"
else
    echo "❌ Multi-patient script still references old master_pipeline.sh"
fi

echo

# Final summary
echo "=== Validation Complete ==="
echo "If all checks show ✅, the directory structure and path references are correctly updated."
echo "The pipeline is ready for testing with the new organization."
echo

# Test dry run capability
echo "--- Testing Dry Run Capability ---"
echo "Testing configuration file access..."

# Check if we can find a test configuration
if [ -f "$ROOT_DIR/configs/analysis/test_analysis.yaml" ]; then
    echo "✅ Test configuration found"
    echo "You can now test the pipeline with:"
    echo "  bash scripts/run_pipeline.sh configs/analysis/test_analysis.yaml --dry-run"
elif [ -f "$ROOT_DIR/configs/analysis/standard_analysis.yaml" ]; then
    echo "✅ Standard configuration found"
    echo "You can now test the pipeline with:"
    echo "  bash scripts/run_pipeline.sh configs/analysis/standard_analysis.yaml --dry-run"
else
    echo "❌ No test configuration found in configs/analysis/"
    echo "Please verify configuration files are properly located"
fi

echo
echo "=== Next Steps ==="
echo "1. Test the pipeline with a dry run:"
echo "   bash scripts/run_pipeline.sh configs/analysis/[config_file] --dry-run"
echo
echo "2. For multi-patient processing:"
echo "   bash scripts/multi_patient.sh data/input/ configs/analysis/template_multi_patient.yaml /path/to/output/ --dry-run"
echo
echo "3. Review the new README.md for complete usage instructions"
echo

exit 0