# TracerX Phylogenetic Analysis Pipeline

A comprehensive computational pipeline for phylogenetic analysis and marker selection in cancer genomics research.

## 📁 Project Structure

```
├── src/                    # Source code organized by analysis stage
│   ├── bootstrap/          # Data bootstrapping and sampling
│   ├── phylowgs/           # Phylogenetic tree inference
│   ├── aggregation/        # Result aggregation and visualization
│   ├── markers/            # Marker selection and optimization
│   ├── longitudinal/       # Longitudinal analysis and tracking
│   └── common/             # Shared utilities and libraries
├── scripts/                # Main execution scripts
│   ├── run_pipeline.sh     # Main pipeline orchestration
│   ├── multi_patient.sh    # Multi-patient batch processing
│   └── install.sh          # Installation script
├── configs/                # Configuration files
│   ├── analysis/           # Analysis-specific configurations
│   └── templates/          # Template configurations
├── data/                   # Data files
│   ├── input/              # Input data files
│   ├── test/               # Test datasets
│   └── raw/                # Raw data batches
└── docs/                   # Documentation
```

## 🚀 Quick Start

### Prerequisites
- High-Performance Computing (HPC) cluster with SLURM
- Conda package manager
- Python 3.x
- PhyloWGS software (installed via `scripts/install.sh`)

### Installation
```bash
# Clone the repository
git clone <repository-url>
cd tracerx-mp

# Install PhyloWGS
bash scripts/install.sh

# Set up conda environments (one for each stage)
# See individual environment.yml files in each src/ subdirectory
```

### Running the Pipeline

#### Single Patient Analysis
```bash
# Using configuration file
bash scripts/run_pipeline.sh configs/analysis/standard_analysis.yaml

# Dry run to test configuration
bash scripts/run_pipeline.sh configs/analysis/test_analysis.yaml --dry-run
```

#### Multi-Patient Batch Processing
```bash
# Process multiple patients
bash scripts/multi_patient.sh data/input/ configs/analysis/template_multi_patient.yaml /path/to/output/

# With delays between submissions
bash scripts/multi_patient.sh data/input/ configs/analysis/template_multi_patient.yaml /path/to/output/ --delay=60
```

#### Longitudinal Analysis
```bash
cd src/longitudinal
python longitudinal_main.py --config configs/templates/cruk0044_fixed_markers.yaml
```

## 🔧 Configuration

### Analysis Configurations (`configs/analysis/`)
- `standard_analysis.yaml` - Production analysis (100 bootstraps)
- `test_analysis.yaml` - Development/testing (20 bootstraps)
- `high_depth_analysis.yaml` - Research-grade analysis (500 bootstraps)

### Template Configurations (`configs/templates/`)
- `template_multi_patient.yaml` - Multi-patient processing template
- `cruk0044_fixed_markers.yaml` - Fixed marker longitudinal analysis
- `cruk0044_dynamic.yaml` - Dynamic marker longitudinal analysis

## 🔬 Pipeline Stages

### 1. Bootstrap (`src/bootstrap/`)
- **Purpose**: Generate multiple bootstrap samples from mutation data
- **Input**: SSM (Somatic Single Mutation) files
- **Output**: Bootstrap samples for phylogenetic analysis

### 2. PhyloWGS (`src/phylowgs/`)
- **Purpose**: Phylogenetic tree inference using PhyloWGS
- **Input**: Bootstrap samples
- **Output**: Phylogenetic trees and mutation assignments

### 3. Aggregation (`src/aggregation/`)
- **Purpose**: Aggregate results and create visualizations
- **Input**: Phylogenetic trees from multiple bootstraps
- **Output**: Best tree structures and visualization plots

### 4. Markers (`src/markers/`)
- **Purpose**: Select optimal genetic markers
- **Input**: Aggregated results and original mutation data
- **Output**: Selected marker lists and optimization results

### 5. Longitudinal (`src/longitudinal/`)
- **Purpose**: Temporal cancer evolution tracking
- **Input**: Tissue data and longitudinal blood samples
- **Output**: Updated tree distributions and temporal analyses

## 📊 Data Organization

### Input Data (`data/input/`)
- SSM files containing mutation data
- Patient-specific datasets
- Configuration-specific input files

### Test Data (`data/test/`)
- Small test datasets for development
- Validation datasets
- Example input files

### Raw Data (`data/raw/`)
- Unprocessed batch data
- Original sequencing results
- Archive of historical datasets

## 🔍 Monitoring and Logging

### Job Monitoring
```bash
# Check job status
squeue -u $USER

# Monitor logs
tail -f /path/to/patient/initial/logs/pipeline_master.log

# Check specific stage logs
tail -f /path/to/patient/initial/logs/slurm_bootstrap_*.out
```

### Log Files
- `pipeline_master.log` - Main pipeline orchestration
- `slurm_*.out` - SLURM standard output
- `slurm_*.err` - SLURM error output
- Stage-specific logs in respective output directories

## 🛠️ Development

### Code Organization
- Each stage is modular and can be run independently
- Consistent naming: `stepN_purpose.py` for Python scripts
- Comprehensive error handling and logging
- YAML-based configuration system

### Testing
```bash
# Run test configuration
bash scripts/run_pipeline.sh configs/analysis/test_analysis.yaml --dry-run

# Test multi-patient processing
bash scripts/multi_patient.sh data/test/ configs/analysis/template_multi_patient_test.yaml /tmp/test_output/ --dry-run
```

## 📈 Performance Optimization

### Resource Allocation
- Stage-specific SLURM resource configurations
- Automatic scaling based on dataset size
- Efficient job dependency management

### Parallel Processing
- Bootstrap stages run in parallel
- Optimized memory usage per stage
- Configurable CPU and memory allocation

## 🆘 Troubleshooting

### Common Issues
1. **Path Resolution**: Ensure all paths are absolute when running on SLURM
2. **Environment Issues**: Verify conda environments are properly activated
3. **Permission Errors**: Check file permissions in output directories
4. **Memory Issues**: Adjust memory allocation in configuration files

### Support
- Check logs in the respective `logs/` directories
- Review configuration files for parameter issues
- Verify input data format and completeness

## 📝 Citation

If you use this pipeline in your research, please cite:

```
[Add citation information here]
```

## 🔄 Changelog

### Version 2.0.0 (Current)
- **Major restructuring**: Organized code into logical directories
- **Enhanced configuration**: YAML-based configuration system
- **Improved documentation**: Comprehensive README and inline documentation
- **Multi-patient support**: Batch processing capabilities
- **Longitudinal analysis**: Temporal evolution tracking

### Previous Versions
- See `docs/changelog.md` for detailed version history

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests and documentation
5. Submit a pull request

## 📄 License

[Add license information here]