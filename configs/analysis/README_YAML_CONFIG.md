# TracerX Marker Selection Pipeline - YAML Configuration Guide

## Overview

The TracerX marker selection pipeline (steps 1-4) now supports unified YAML configuration, replacing complex command-line arguments with structured, readable configuration files. This system provides:

- **Single configuration file** for entire pipeline execution
- **Environment-specific** configurations for testing, standard analysis, and high-performance computing
- **Comprehensive validation** and error handling
- **Dry-run capability** for testing configurations
- **Resource optimization** with step-specific SLURM settings

## Quick Start

### 1. Choose or Create Configuration

```bash
# Use existing template
cp configs/standard_analysis.yaml configs/my_analysis.yaml

# Edit paths and patient information
vim configs/my_analysis.yaml
```

### 2. Execute Pipeline

```bash
# Submit to SLURM
sbatch master_pipeline.sh configs/my_analysis.yaml

# Test configuration (dry run)
bash master_pipeline.sh configs/my_analysis.yaml --dry-run
```

### 3. Monitor Progress

```bash
# Check job status
squeue -u $USER

# View logs
tail -f /path/to/patient_data/PATIENT_ID/initial/logs/pipeline_master.log
```

## Configuration Templates

### Standard Analysis (`standard_analysis.yaml`)

Complete production analysis with default settings:
- **100 bootstraps** for robust statistics
- **Standard resource allocation** (2-16GB memory per step)
- **Full marker selection** with multiple optimization strategies
- **Comprehensive visualizations**

**Use case**: Regular production analysis of patient samples

### Test Analysis (`test_analysis.yaml`)

Rapid testing configuration for development:
- **20 bootstraps** for faster execution
- **Reduced resources** (1-8GB memory per step)
- **Limited marker iterations** for quick validation
- **Debug logging enabled**

**Use case**: Development, debugging, and rapid validation

### High-Depth Analysis (`high_depth_analysis.yaml`)

Intensive computational analysis:
- **500 bootstraps** for maximum robustness
- **High resources** (8-32GB memory per step)
- **Multiple optimization strategies**
- **Extended walltime** for thorough analysis

**Use case**: Research studies requiring maximum statistical power

## Configuration Structure

### Required Sections

```yaml
# Patient identification
patient_id: "CRUK0044"
analysis_type: "standard"  # standard, test, high_depth

# Input files (REQUIRED - update paths for your environment)
input:
  ssm_file: "/path/to/data/ssm.txt"
  code_dir: "/path/to/tracerx-mp/"

# Output configuration
output:
  base_dir: "/path/to/patient_data/CRUK0044/"
  create_subdirs: true
  log_level: "INFO"
```

### Step-Specific Configuration

#### Bootstrap (Step 1)
```yaml
bootstrap:
  num_bootstraps: 100    # Number of bootstrap samples
  random_seed: null      # Set integer for reproducible results
```

#### PhyloWGS (Step 2)
```yaml
phylowgs:
  num_chains: 5          # MCMC chains
  parallel_limit: 10     # SLURM array job throttling
```

#### Aggregation (Step 3)
```yaml
aggregation:
  method: "phylowgs"
  generate_visualizations: true
  sample_prefix: "Region"
  custom_sample_names: []  # Optional custom names
```

#### Marker Selection (Step 4)
```yaml
marker_selection:
  read_depth: 1500
  filter_strategy: "any_high"    # any_high, all_high, majority_high
  filter_threshold: 0.9
  filter_samples: []             # For specific_samples strategy
  
  optimization:
    lambda1_values: [1.0, 0.0]    # Fraction weights
    lambda2_values: [0.0, 1.0]    # Structure weights
    max_iterations: null          # null = all markers
```

### HPC/SLURM Configuration

Each step can have custom resource allocation:

```yaml
hpc:
  bootstrap:
    partition: "pool1"
    cpus_per_task: 2
    memory: "8G"
    walltime: "02:00:00"
    conda_env: "preprocess_env"
  
  marker_selection:
    partition: "pool1"
    cpus_per_task: 1
    memory: "16G"
    walltime: "06:00:00"
    conda_env: "markers_env"
    modules: ["gurobi1102"]  # Required modules
```

## Configuration Validation

The pipeline automatically validates configurations and provides clear error messages:

### Path Validation
- Checks existence of input files and directories
- Converts relative paths to absolute paths
- Validates code directory structure

### Parameter Validation
- Ensures required fields are present
- Validates parameter ranges and types
- Checks HPC resource specifications

### Dry Run Testing
```bash
# Test configuration without submitting jobs
bash master_pipeline.sh configs/my_analysis.yaml --dry-run
```

## Advanced Configuration

### Custom Sample Names
```yaml
aggregation:
  sample_prefix: "Tumor"
  custom_sample_names: ["Primary", "Metastasis_1", "Metastasis_2"]
```

### Filtering Strategies
```yaml
marker_selection:
  filter_strategy: "specific_samples"
  filter_samples: [0, 1]  # Only use first two samples
  filter_threshold: 0.9
```

### Resource Optimization
```yaml
# High-memory configuration for large datasets
hpc:
  aggregation:
    memory: "64G"
    cpus_per_task: 4
  
  marker_selection:
    memory: "128G"
    cpus_per_task: 8
    walltime: "24:00:00"
```

## Troubleshooting

### Common Configuration Errors

1. **Missing required fields**
   ```
   Error: Missing required configuration parameters:
     PATIENT_ID: ''
   ```
   **Solution**: Ensure all required fields are filled in configuration

2. **File not found errors**
   ```
   Error: Input SSM file not found: /path/to/data/ssm.txt
   ```
   **Solution**: Update file paths in configuration to match your environment

3. **Invalid YAML syntax**
   ```
   Error parsing configuration: while parsing a block mapping...
   ```
   **Solution**: Check YAML syntax, ensure proper indentation and no tab characters

### Debugging Tips

1. **Use dry-run mode**
   ```bash
   bash master_pipeline.sh configs/my_analysis.yaml --dry-run
   ```

2. **Enable debug logging**
   ```yaml
   output:
     log_level: "DEBUG"
   validation:
     debug_mode: true
   ```

3. **Start with test configuration**
   ```bash
   # Test with smaller dataset first
   cp configs/test_analysis.yaml configs/my_test.yaml
   # Edit paths and run
   sbatch master_pipeline.sh configs/my_test.yaml
   ```

## Migration from Command Line

### Old Command Line Approach
```bash
bash main_init.sh CRUK0044 data/ssm.txt /path/output/ 100 1500
```

### New YAML Approach
```yaml
# configs/cruk0044.yaml
patient_id: "CRUK0044"
input:
  ssm_file: "data/ssm.txt"
  code_dir: "/path/to/tracerx-mp/"
output:
  base_dir: "/path/output/"
bootstrap:
  num_bootstraps: 100
marker_selection:
  read_depth: 1500
```

```bash
sbatch master_pipeline.sh configs/cruk0044.yaml
```

## Performance Optimization

### Resource Scaling Guidelines

| Dataset Size | Bootstrap | PhyloWGS | Aggregation | Marker Selection |
|--------------|-----------|----------|-------------|------------------|
| Small (<10 mutations) | 2G, 1 CPU | 4G, 2 CPU | 8G, 1 CPU | 8G, 1 CPU |
| Medium (10-50 mutations) | 8G, 2 CPU | 8G, 4 CPU | 16G, 1 CPU | 16G, 2 CPU |
| Large (>50 mutations) | 16G, 4 CPU | 16G, 8 CPU | 32G, 2 CPU | 32G, 4 CPU |

### Walltime Recommendations

| Analysis Type | Total Runtime | Critical Path |
|---------------|---------------|---------------|
| Test | 2-4 hours | PhyloWGS arrays |
| Standard | 8-12 hours | PhyloWGS arrays |
| High-depth | 24-48 hours | PhyloWGS arrays |

## Best Practices

1. **Version control configurations**: Track configuration files in git
2. **Environment-specific configs**: Separate configs for dev/test/prod
3. **Resource monitoring**: Check `seff` after completion to optimize resources
4. **Incremental testing**: Start with test config, then scale up
5. **Documentation**: Comment complex configurations for team use

## Examples

See the `configs/` directory for complete working examples:
- `standard_analysis.yaml`: Production analysis template
- `test_analysis.yaml`: Development and testing template  
- `high_depth_analysis.yaml`: Research-grade analysis template

Each template includes detailed comments explaining all configuration options.