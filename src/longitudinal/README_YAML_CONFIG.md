# Longitudinal Analysis with YAML Configuration

This directory contains the updated longitudinal cancer evolution analysis pipeline that uses YAML configuration files for clean parameter management, eliminating the complexity of shell argument parsing.

## Overview

The YAML-based approach provides several advantages:
- **Clean Configuration**: All parameters in readable YAML format
- **Version Control**: Configuration files can be tracked and versioned
- **No Shell Escaping**: Eliminates issues with special characters like `>` in mutation names
- **Reproducibility**: Exact configurations can be easily shared and rerun
- **Validation**: Built-in configuration validation and error reporting

## Files

### Core Scripts
- `longitudinal_update.py` - Main Python analysis script (updated for YAML)
- `longitudinal_analysis_yaml.sh` - Simplified SLURM submission script
- `config_schema.yaml` - Configuration schema and documentation

### Configuration Files
- `configs/cruk0044_fixed_markers.yaml` - Fixed marker analysis for CRUK0044
- `configs/cruk0044_dynamic.yaml` - Dynamic marker analysis for CRUK0044  
- `configs/cruk0044_both.yaml` - Comparative analysis (both approaches)

## Usage

### 1. Prepare Configuration File

Copy and modify one of the example configuration files:

```bash
# Copy template
cp configs/cruk0044_fixed_markers.yaml configs/my_patient_analysis.yaml

# Edit paths and parameters
nano configs/my_patient_analysis.yaml
```

Update the following sections in your config file:
- `patient_id`: Your patient identifier
- `input_files`: Paths to your aggregation results, SSM file, longitudinal data, and code directory
- `output.base_dir`: Where to save results
- `fixed_markers`: Gene identifiers (for fixed mode)
- `parameters`: Analysis parameters as needed

### 2. Run Analysis

#### Local Execution
```bash
# Run with YAML config
python longitudinal_update.py --config configs/my_patient_analysis.yaml

# Run with debug mode
python longitudinal_update.py --config configs/my_patient_analysis.yaml --debug

# Run without plots
python longitudinal_update.py --config configs/my_patient_analysis.yaml --no-plots
```

#### SLURM Execution
```bash
# Basic submission
sbatch longitudinal_analysis_yaml.sh configs/my_patient_analysis.yaml

# With debug mode
sbatch longitudinal_analysis_yaml.sh configs/my_patient_analysis.yaml "--debug"

# With multiple flags and custom log suffix
sbatch longitudinal_analysis_yaml.sh configs/my_patient_analysis.yaml "--debug --no-plots" "test_run"
```

## Configuration Structure

### Required Sections

```yaml
# Patient identification
patient_id: "PATIENT_ID"
analysis_mode: "fixed"  # or "dynamic" or "both"

# Input file paths
input_files:
  aggregation_dir: "/path/to/aggregation_results"
  ssm_file: "/path/to/ssm.txt"
  longitudinal_data: "/path/to/longitudinal_data.csv"
  code_dir: "/path/to/tracerx-mp"

# Output configuration
output:
  base_dir: "/path/to/output"
  create_subdirs: true
```

### Analysis Parameters

```yaml
parameters:
  n_markers: 2          # Number of markers (dynamic mode)
  read_depth: 90000     # Expected ddPCR read depth
  method: "phylowgs"    # Phylogenetic method
  lambda1: 0.0         # Weight for fraction-based objective (0.0 = ignore fractions)
  lambda2: 1.0         # Weight for structure-based objective (1.0 = focus on tree structure)
  focus_sample: 0      # Sample focus index
```

**Lambda Parameter Guide:**
- `lambda1=1.0, lambda2=0.0` → Pure fraction-based optimization (focus on VAF differences)
- `lambda1=0.0, lambda2=1.0` → Pure structure-based optimization (focus on tree topology)
- `lambda1=0.5, lambda2=0.5` → Balanced approach (consider both fractions and structure)
- `lambda1=0.0, lambda2=0.0` → Invalid (no objective function)

### Fixed Markers (for fixed/both modes)

```yaml
fixed_markers:
  - "DLG2_11_83544685_T>A"
  - "GPR65_14_88477948_G>T"
  - "C12orf74_12_93100715_G>T"
  - "CSMD1_8_2975974_G>T"
  - "OR51D1_11_4661941_G>T"
```

### Optional Sections

```yaml
# Timepoint filtering
filtering:
  timepoints: []  # Empty = all timepoints
  # timepoints: ["2021-01-15", "2021-03-15"]  # Specific dates

# Visualization options
visualization:
  generate_plots: true
  plot_format: "png"
  save_intermediate: false

# Validation and debugging
validation:
  validate_inputs: true
  debug_mode: false
```

## Analysis Modes

### Fixed Marker Analysis
- Uses clinically specified markers consistently across all timepoints
- Ideal for standardized clinical monitoring
- Configuration: `analysis_mode: "fixed"`

### Dynamic Marker Analysis  
- Selects optimal markers at each timepoint based on tree discrimination
- Maximizes information gain but requires different assays per timepoint
- Configuration: `analysis_mode: "dynamic"`

### Comparative Analysis
- Runs both approaches and generates comparison metrics
- Useful for method evaluation and clinical decision making
- Configuration: `analysis_mode: "both"`

## Output Structure

Results are organized in the specified output directory:

```
output_directory/
├── logs/                           # Execution logs
├── fixed_marker_analysis/          # Fixed approach results (if applicable)
│   ├── updated_trees/
│   ├── marker_tracking/
│   └── fixed_marker_complete_results.json
├── dynamic_marker_analysis/        # Dynamic approach results (if applicable)
│   ├── updated_trees/
│   ├── marker_selections/
│   └── dynamic_marker_complete_results.json
├── comparative_analysis/           # Comparison results (if both modes)
└── analysis_summary.json          # Overall summary
```

## Troubleshooting

### Configuration Errors
- Check YAML syntax with `python -c "import yaml; yaml.safe_load(open('config.yaml'))"`
- Ensure all required sections are present
- Verify file paths exist and are accessible

### Missing Dependencies
- Ensure `pyyaml` is installed: `pip install pyyaml`
- Verify all other dependencies are available in your conda environment

### Path Issues
- Use absolute paths in configuration files
- Ensure output directories exist or can be created
- Check file permissions for input and output locations

## Migration from Command Line

If you have existing command line calls, convert them to YAML:

**Old:**
```bash
python longitudinal_update.py CRUK0044 \
    --aggregation-dir /data/aggregation \
    --ssm-file /data/ssm.txt \
    --longitudinal-data /data/liquid.csv \
    --output-dir /data/output \
    --analysis-mode fixed \
    --fixed-markers "DLG2_11_83544685_T>A" "GPR65_14_88477948_G>T" \
    --read-depth 90000
```

**New:**
```yaml
# config.yaml
patient_id: "CRUK0044"
analysis_mode: "fixed"
input_files:
  aggregation_dir: "/data/aggregation"
  ssm_file: "/data/ssm.txt"
  longitudinal_data: "/data/liquid.csv"
  code_dir: "/path/to/tracerx-mp"
output:
  base_dir: "/data/output"
parameters:
  read_depth: 90000
  lambda1: 0.0  # Structure-focused optimization
  lambda2: 1.0
fixed_markers:
  - "DLG2_11_83544685_T>A"
  - "GPR65_14_88477948_G>T"
```

```bash
python longitudinal_update.py --config config.yaml
```

This approach is much cleaner and eliminates the shell parsing issues you encountered! 