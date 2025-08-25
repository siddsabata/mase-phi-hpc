# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Common Development Commands

### Pipeline Execution
```bash
# Run main pipeline with configuration
bash scripts/run_pipeline.sh configs/analysis/standard_analysis.yaml

# Test configuration (dry run)
bash scripts/run_pipeline.sh configs/analysis/test_analysis.yaml --dry-run

# Multi-patient batch processing
bash scripts/multi_patient.sh data/input/ configs/analysis/template_multi_patient.yaml /path/to/output/

# Longitudinal analysis
cd src/longitudinal
python longitudinal_main.py --config configs/templates/cruk0044_fixed_markers.yaml
```

### Environment Setup
```bash
# Install PhyloWGS and dependencies
bash scripts/install.sh

# Validate setup
bash scripts/validate_setup.sh

# Conda environments are defined in environment.yml files in each src/ subdirectory
# Each pipeline stage has its own environment (bootstrap_env, phylowgs_env, aggregation_env, markers_env)
```

### Job Monitoring (SLURM)
```bash
# Check job status
squeue -u $USER

# Monitor main pipeline log
tail -f /path/to/patient/initial/logs/pipeline_master.log

# Check stage-specific logs
tail -f /path/to/patient/initial/logs/slurm_bootstrap_*.out
```

## Code Architecture

### Pipeline Structure
The codebase implements a 5-stage phylogenetic analysis pipeline for cancer genomics:

1. **Bootstrap** (`src/bootstrap/`): Generates multiple bootstrap samples from mutation data (SSM files)
2. **PhyloWGS** (`src/phylowgs/`): Phylogenetic tree inference using external PhyloWGS software
3. **Aggregation** (`src/aggregation/`): Aggregates results from multiple bootstraps and creates visualizations
4. **Markers** (`src/markers/`): Selects optimal genetic markers based on tree structures
5. **Longitudinal** (`src/longitudinal/`): Temporal cancer evolution tracking with tissue and blood samples

### Key Architecture Patterns

- **SLURM Integration**: Each stage has a corresponding `.sh` script that handles SLURM job submission with stage-specific resource requirements
- **YAML Configuration**: All pipeline parameters are defined in YAML files with schema validation
- **Modular Design**: Each stage can run independently with clear input/output interfaces
- **Parallel Processing**: Bootstrap stages run as SLURM array jobs; dependency management ensures proper execution order

### Configuration System

- `configs/analysis/`: Analysis-specific configurations (standard_analysis.yaml, test_analysis.yaml, high_depth_analysis.yaml)
- `configs/templates/`: Template configurations for multi-patient and longitudinal analysis
- Configuration parser in `scripts/run_pipeline.sh` handles YAML to shell variable conversion
- Each stage has HPC-specific settings (partition, memory, walltime, conda environment)

### Data Flow

1. **Input**: SSM (Somatic Single Mutation) files containing mutation data
2. **Bootstrap**: Creates multiple resampled datasets → `bootstraps/` directories
3. **PhyloWGS**: Generates phylogenetic trees → `phylowgs_output/` in each bootstrap
4. **Aggregation**: Combines trees across bootstraps → `aggregation_results/`
5. **Markers**: Selects optimal markers → `markers/` directory
6. **Longitudinal**: Updates trees with temporal data → analysis-specific output directories

### Longitudinal Analysis Module

The longitudinal analysis (`src/longitudinal/`) implements a modular design:

- **config_handler.py**: YAML configuration parsing and validation
- **data_loader.py**: Standardized data loading for tree distributions, SSM files, and longitudinal data
- **fixed_analysis.py**: Analysis with predefined marker sets
- **dynamic_analysis.py**: Dynamic marker selection and optimization
- **longitudinal_main.py**: Main orchestration script supporting both analysis modes
- **clone_frequency.py**: Clone frequency tracking and visualization when enabled

### Important Implementation Notes

- All paths in SLURM jobs must be absolute due to execution context
- Each stage activates its own conda environment defined in `environment.yml`
- Pipeline uses job dependencies (`afterok`/`afterany`) to ensure proper execution order
- Error handling allows partial bootstrap failures while continuing aggregation
- Memory and CPU allocation is configurable per stage via YAML configuration
- The codebase is designed for HPC clusters with SLURM workload manager

### Testing and Validation

- Use `--dry-run` flag to validate configurations without job submission
- `test_analysis.yaml` provides faster testing with reduced bootstrap counts (20 vs 100)
- `scripts/validate_setup.sh` verifies directory structure and path references
- Each stage logs to both SLURM output files and centralized pipeline logs