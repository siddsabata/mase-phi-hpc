# Multi-Patient Pipeline Implementation Guide

## Overview
This guide provides step-by-step instructions to implement a multi-patient processing system for the TracerX Marker Selection Pipeline. The system will process multiple SSM files using the existing single-patient `master_pipeline.sh` without any modifications to existing code.

## Prerequisites
- Existing TracerX-MP pipeline is working (test with `master_pipeline.sh`)
- Access to SLURM cluster with array job support
- Understanding of bash scripting and YAML configuration files
- Familiarity with the existing pipeline structure

## Implementation Steps

### Step 1: Create Main Multi-Patient Script
**File**: `multi_patient_pipeline.sh`
**Estimated time**: 2-3 hours

#### Task 1.1: Create script skeleton
```bash
#!/bin/bash
# Multi-Patient Pipeline Script for TracerX Marker Selection Pipeline
# Add header with usage examples and documentation
```

#### Task 1.2: Implement input validation and argument parsing
- Parse 3 required arguments: `ssm_directory`, `config_template`, `output_base_directory`
- Handle optional arguments: `--delay=N`, `--dry-run`, `--help`
- Validate all input files and directories exist
- Convert relative paths to absolute paths

#### Task 1.3: Implement directory structure creation
```bash
# Create these directories:
MULTI_CONFIGS_DIR="${OUTPUT_BASE_DIRECTORY}/configs/generated"
MULTI_PATIENTS_DIR="${OUTPUT_BASE_DIRECTORY}/patients"
MULTI_LOGS_DIR="${OUTPUT_BASE_DIRECTORY}/logs"
```

#### Task 1.4: Implement SSM file discovery
- Find all `*.txt` and `*.ssm` files in the SSM directory
- Sort files for consistent processing order
- Validate at least one SSM file is found

#### Task 1.5: Implement patient ID extraction function
```bash
extract_patient_id() {
    local ssm_file="$1"
    local filename=$(basename "$ssm_file")
    # Remove extension: .txt, .ssm
    # Remove common suffixes: _ssm, _subset
    # Return clean patient ID
}
```

#### Task 1.6: Implement configuration generation function
```bash
generate_patient_config() {
    # Use sed to replace placeholders:
    # PLACEHOLDER_PATIENT_ID → actual patient ID
    # PLACEHOLDER_SSM_FILE → absolute path to SSM file
    # PLACEHOLDER_OUTPUT_DIR → patient-specific output directory
    # PLACEHOLDER_CODE_DIR → TracerX-MP repository directory
}
```

#### Task 1.7: Implement pipeline submission function
```bash
submit_patient_pipeline() {
    # Call: bash master_pipeline.sh patient_config.yaml
    # Capture output to patient-specific log file
    # Return success/failure status
    # Handle dry-run mode
}
```

#### Task 1.8: Implement main processing loop
- Iterate through all SSM files
- For each patient:
  - Extract patient ID
  - Create patient output directory
  - Generate patient-specific config
  - Submit patient pipeline
  - Apply delay if specified
  - Track success/failure counts

#### Task 1.9: Implement logging and status reporting
- Setup master log file with timestamp
- Log all operations and decisions
- Provide final summary with success/failure counts
- Include monitoring instructions for users

**Validation for Step 1:**
```bash
# Test script creation and basic validation
bash multi_patient_pipeline.sh --help
bash multi_patient_pipeline.sh nonexistent_dir template.yaml output/ --dry-run
# Should show appropriate error messages
```

### Step 2: Create Configuration Templates
**Files**: `configs/template_multi_patient.yaml`, `configs/template_multi_patient_test.yaml`
**Estimated time**: 1 hour

#### Task 2.1: Create standard production template
**File**: `configs/template_multi_patient.yaml`

- Copy structure from existing `configs/standard_analysis.yaml`
- Replace patient-specific values with placeholders:
  ```yaml
  patient_id: "PLACEHOLDER_PATIENT_ID"
  ssm_file: "PLACEHOLDER_SSM_FILE"
  base_dir: "PLACEHOLDER_OUTPUT_DIR"
  code_dir: "PLACEHOLDER_CODE_DIR"
  ```
- Use conservative resource settings for multi-patient workload:
  - `parallel_limit: 20` (instead of 50)
  - `num_bootstraps: 100` (standard)
  - Memory allocations: 8-16GB per step
- Add comments explaining multi-patient specific settings

#### Task 2.2: Create testing template
**File**: `configs/template_multi_patient_test.yaml`

- Copy structure from `configs/test_analysis.yaml`
- Replace patient-specific values with same placeholders
- Use reduced resources for testing:
  - `num_bootstraps: 20`
  - `parallel_limit: 5`
  - Memory allocations: 4-8GB per step
  - Shorter walltimes: 30min-1hour per step
- Enable debug logging: `log_level: "DEBUG"`

**Validation for Step 2:**
```bash
# Test YAML syntax
python -c "import yaml; yaml.safe_load(open('configs/template_multi_patient.yaml'))"
python -c "import yaml; yaml.safe_load(open('configs/template_multi_patient_test.yaml'))"

# Test placeholder substitution
sed -e "s|PLACEHOLDER_PATIENT_ID|TEST|g" \
    -e "s|PLACEHOLDER_SSM_FILE|/test.txt|g" \
    -e "s|PLACEHOLDER_OUTPUT_DIR|/output|g" \
    -e "s|PLACEHOLDER_CODE_DIR|/code|g" \
    configs/template_multi_patient.yaml
```

### Step 3: Create Documentation
**File**: `README_MULTI_PATIENT.md`
**Estimated time**: 2 hours

#### Task 3.1: Write comprehensive user guide
Include sections:
- Overview and key features
- Directory structure explanation
- Usage examples with command-line options
- Input requirements and file formats
- Configuration template documentation
- Resource management and SLURM considerations
- Monitoring and debugging instructions
- Performance considerations and timing estimates
- Troubleshooting guide with common issues

#### Task 3.2: Add integration examples
- How to convert existing single-patient configs to templates
- How to use with existing data
- How to customize for different analysis types

#### Task 3.3: Add testing instructions
- Step-by-step testing procedures
- Dry-run examples
- Small batch testing recommendations

**Validation for Step 3:**
```bash
# Review documentation completeness
# Ensure all command examples are accurate
# Verify all file paths and references are correct
```

### Step 4: Create Test Utilities
**File**: `data/create_test_patients.sh`
**Estimated time**: 30 minutes

#### Task 4.1: Create test data generation script
```bash
#!/bin/bash
# Create multiple SSM files for testing:
# - Copy existing ssm.txt as CRUK0044.txt
# - Copy ssm_subset.txt as multiple patients with different names
# - Test different file extensions (.txt, .ssm)
# - Test suffix removal (patient_subset.txt → patient)
```

#### Task 4.2: Make script executable and test
```bash
chmod +x data/create_test_patients.sh
bash data/create_test_patients.sh
ls -la data/test_patients/
```

**Validation for Step 4:**
```bash
# Run script and verify it creates multiple patient files
bash data/create_test_patients.sh
# Should create 5+ test patient files with different naming patterns
```

### Step 5: Integration Testing
**Estimated time**: 2-3 hours

#### Task 5.1: Test configuration generation
```bash
# Create test patients
bash data/create_test_patients.sh

# Test dry-run mode
bash multi_patient_pipeline.sh data/test_patients/ configs/template_multi_patient_test.yaml test_output/ --dry-run

# Verify:
# - All SSM files are discovered
# - Patient IDs are extracted correctly
# - Configurations would be generated properly
# - No actual jobs submitted
```

#### Task 5.2: Test single patient submission
```bash
# Test with one patient using testing template
mkdir -p data/single_test/
cp data/test_patients/CRUK0044.txt data/single_test/

bash multi_patient_pipeline.sh data/single_test/ configs/template_multi_patient_test.yaml test_single/

# Verify:
# - Configuration is generated correctly
# - master_pipeline.sh is called with correct arguments
# - Jobs are submitted to SLURM
# - Directory structure is created properly
```

#### Task 5.3: Test small batch
```bash
# Test with 2-3 patients
mkdir -p data/small_batch/
cp data/test_patients/CRUK004{4,5}.txt data/small_batch/

bash multi_patient_pipeline.sh data/small_batch/ configs/template_multi_patient_test.yaml test_batch/ --delay=30

# Verify:
# - Multiple patients processed sequentially
# - Delay is applied between submissions
# - Each patient gets independent directory structure
# - Logs are properly separated
```

#### Task 5.4: Test error handling
```bash
# Test with invalid inputs
bash multi_patient_pipeline.sh nonexistent/ template.yaml output/
bash multi_patient_pipeline.sh data/ nonexistent.yaml output/
bash multi_patient_pipeline.sh data/ template.yaml /root/forbidden/

# Verify appropriate error messages
```

#### Task 5.5: Test monitoring and logs
```bash
# During a test run, verify:
tail -f test_output/logs/multi_patient_master.log
tail -f test_output/logs/CRUK0044_submission.log
ls -la test_output/patients/CRUK0044/initial/logs/

# Check SLURM job status
squeue -u $USER
```

**Validation for Step 5:**
- [ ] Dry-run mode works correctly
- [ ] Single patient submission completes successfully
- [ ] Multi-patient batch processes correctly
- [ ] Error handling works as expected
- [ ] Logging is comprehensive and useful
- [ ] Directory structure matches design
- [ ] SLURM jobs are submitted with correct parameters

### Step 6: Performance Testing
**Estimated time**: 1-2 hours

#### Task 6.1: Test resource management
```bash
# Test with testing template (reduced resources)
bash multi_patient_pipeline.sh data/test_patients/ configs/template_multi_patient_test.yaml perf_test/

# Monitor cluster resource usage:
squeue -u $USER
sinfo
htop  # if available
```

#### Task 6.2: Test delay functionality
```bash
# Test different delay settings
bash multi_patient_pipeline.sh data/test_patients/ configs/template_multi_patient_test.yaml delay_test/ --delay=60

# Verify jobs are submitted with appropriate spacing
```

#### Task 6.3: Validate array job independence
```bash
# Submit multiple patients and verify:
# - Each gets separate array job IDs
# - Array jobs don't conflict
# - Total array elements don't exceed limits
squeue -u $USER | grep ARRAY
```

**Validation for Step 6:**
- [ ] Resource usage is reasonable for cluster
- [ ] Delay functionality works correctly
- [ ] Array jobs are properly isolated
- [ ] No cluster limit violations

### Step 7: Documentation Updates
**Estimated time**: 1 hour

#### Task 7.1: Update main project documentation
Add multi-patient pipeline information to:
- `project.md` - Add section describing multi-patient capabilities
- `README.md` - Add usage examples for multi-patient processing

#### Task 7.2: Create quick-start guide
Create a simple example workflow:
```bash
# Quick-start example
bash data/create_test_patients.sh
bash multi_patient_pipeline.sh data/test_patients/ configs/template_multi_patient_test.yaml my_analysis/ --dry-run
bash multi_patient_pipeline.sh data/test_patients/ configs/template_multi_patient_test.yaml my_analysis/
```

#### Task 7.3: Add troubleshooting section
Document common issues and solutions:
- "No SSM files found" → check file extensions
- "Pipeline submission failed" → check individual patient logs
- "Array jobs exceed limits" → use testing template or delays

**Validation for Step 7:**
- [ ] Documentation is complete and accurate
- [ ] Examples can be copy-pasted and work
- [ ] Troubleshooting guide addresses likely issues

### Step 8: Final Validation and Cleanup
**Estimated time**: 1 hour

#### Task 8.1: End-to-end testing
```bash
# Complete workflow test:
1. Create test data
2. Run multi-patient pipeline with testing template
3. Monitor progress through completion
4. Verify all expected outputs are generated
5. Check that results match single-patient runs
```

#### Task 8.2: Code review checklist
- [ ] All scripts have proper error handling
- [ ] All file paths use absolute paths
- [ ] All user inputs are validated
- [ ] All functions have clear purposes
- [ ] Code follows existing project conventions
- [ ] Comments explain complex logic

#### Task 8.3: File permissions and executable bits
```bash
chmod +x multi_patient_pipeline.sh
chmod +x data/create_test_patients.sh
# Ensure all scripts are executable
```

#### Task 8.4: Clean up test files
```bash
# Remove temporary test directories if desired
rm -rf test_output/ test_single/ test_batch/ perf_test/ delay_test/
# Keep data/test_patients/ for future testing
```

## Deliverables Checklist

### Code Files
- [ ] `multi_patient_pipeline.sh` - Main orchestration script
- [ ] `configs/template_multi_patient.yaml` - Production template
- [ ] `configs/template_multi_patient_test.yaml` - Testing template
- [ ] `data/create_test_patients.sh` - Test data generation utility

### Documentation
- [ ] `README_MULTI_PATIENT.md` - Comprehensive user guide
- [ ] Updates to `project.md` and main `README.md`
- [ ] This implementation guide (`todo.md`)

### Testing
- [ ] Test data files in `data/test_patients/`
- [ ] Validation of dry-run functionality
- [ ] Validation of single and multi-patient workflows
- [ ] Performance and resource usage testing

## Success Criteria

### Functional Requirements
- [ ] Process multiple SSM files from a directory
- [ ] Generate patient-specific configurations from templates
- [ ] Submit independent pipelines for each patient
- [ ] Maintain complete isolation between patients
- [ ] Support dry-run testing mode
- [ ] Provide comprehensive logging

### Non-Functional Requirements
- [ ] Zero modifications to existing `master_pipeline.sh`
- [ ] Cluster-friendly resource management
- [ ] Clear error messages and troubleshooting
- [ ] Performance suitable for 10+ patients
- [ ] Documentation sufficient for independent use

### Integration Requirements
- [ ] Works with existing YAML configuration system
- [ ] Compatible with existing SLURM job submission
- [ ] Maintains existing directory structures
- [ ] Supports existing conda environments

## Timeline Estimate

**Total estimated time: 10-15 hours**

- **Day 1 (4-5 hours)**: Steps 1-2 (Main script + templates)
- **Day 2 (3-4 hours)**: Steps 3-4 (Documentation + test utilities)
- **Day 3 (3-4 hours)**: Steps 5-6 (Integration + performance testing)
- **Day 4 (1-2 hours)**: Steps 7-8 (Final documentation + validation)

## Risk Mitigation

### Potential Issues
1. **SLURM array job limits**: Mitigated by conservative template settings
2. **Cluster overload**: Mitigated by sequential submission and delay options
3. **Path resolution**: Mitigated by absolute path conversion
4. **Configuration errors**: Mitigated by dry-run mode and validation

### Rollback Plan
- Keep all existing code unchanged
- Multi-patient functionality is additive only
- Can disable/remove without affecting existing workflows

## Post-Implementation

### User Training
- Demonstrate multi-patient workflow to team
- Provide examples with real patient data
- Document any site-specific customizations needed

### Monitoring
- Monitor cluster resource usage with multi-patient workloads
- Collect user feedback on performance and usability
- Iterate on template configurations based on experience

### Future Enhancements
- Web interface for multi-patient submission
- Integration with workflow management systems
- Advanced patient filtering and selection options
- Automated result aggregation across patients
