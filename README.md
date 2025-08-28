# Mase-phi HPC

A phylogenetic analysis pipeline for cancer genomics that processes WGS mutation data (SSM files) through bootstrap sampling, phylogenetic tree inference, aggregation, and marker selection.

## Installation

### Requirements
- **Conda is required** - Install [Miniforge](https://github.com/conda-forge/miniforge) (recommended) or Miniconda
- SLURM workload manager
- Gurobi optimization software (for marker selection)

### Quick Setup
```bash
# Clone repository
git clone <repository-url>
cd mase-phi-hpc

# Install pipeline (creates environments and dependencies)
bash scripts/install.sh
```

The installer automatically uses `mamba` for faster environment creation (10-100x speedup) if available, with fallback to `conda`.

## Pipeline Overview

The pipeline consists of 4 main steps:
1. **Bootstrap**: Generate multiple bootstrap samples from mutation data with VAF filtering
2. **PhyloWGS**: Phylogenetic tree inference using MCMC
3. **Aggregation**: Combine bootstrap results and create visualizations  
4. **Marker Selection**: Select optimal genetic markers using two optimization strategies

## Single Patient Processing

### Configuration
Edit `configs/config_single.yaml` with your paths and settings.

### Execution
```bash
bash scripts/run_pipeline.sh configs/config_single.yaml
```

### Output Structure
```
{base_dir}/{patient_id}/initial/
├── ssm_filtered.txt        # VAF-filtered mutations
├── bootstraps/             # Bootstrap samples with PhyloWGS results
├── aggregation_results/    # Combined tree structures  
├── markers/               # Marker selection results
└── logs/                  # Pipeline logs
```

## Multi-Patient Processing

### Configuration  
Edit `configs/config_multi.yaml` with your settings.

### Execution
```bash
bash scripts/multi_patient.sh <ssm_directory> configs/config_multi.yaml <output_base_directory>

# With delay between submissions (recommended for large batches)
bash scripts/multi_patient.sh data/patients/ configs/config_multi.yaml /results/ --delay=60
```

### Input Structure
```
data/patients/
├── patient_001.txt
├── patient_002.txt  
└── patient_003.txt
```

### Output Structure
```
/results/
├── patients/              
│   ├── patient_001/initial/    # Individual pipeline results
│   ├── patient_002/initial/
│   └── patient_003/initial/
├── configs/generated/          # Auto-generated patient configs
└── logs/                       # Multi-patient orchestration logs
```

## Marker Selection Output

- **λ1=1, λ2=0 (Pure Fraction Optimization)**: Optimizes for VAF-based marker selection (best for clone frequency tracking)
- **λ1=0, λ2=1 (Pure Structure Optimization)**: Optimizes for phylogenetic tree structure differences (best for phylogenetic reconstruction)

## Job Monitoring

```bash
# Check SLURM job status
squeue -u $USER
```

## Troubleshooting

- Ensure Gurobi modules are loaded correctly in HPC configuration
- Check that SSM files have required columns: `id`, `gene`, `a`, `d`, `mu_r`, `mu_v`  
- Verify sufficient disk space for bootstrap samples and tree results
- Review stage-specific logs for detailed error information
