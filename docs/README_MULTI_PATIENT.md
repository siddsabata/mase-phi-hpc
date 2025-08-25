# Multi-Patient Pipeline for TracerX Marker Selection

## Overview

The multi-patient pipeline extends the existing single-patient TracerX marker selection pipeline to process multiple patients in a coordinated manner. This system maintains **100% compatibility** with the existing `master_pipeline.sh` while adding batch processing capabilities.

## Key Features

- **Zero Changes to Existing Code**: The multi-patient system uses the existing `master_pipeline.sh` without any modifications
- **Template-Based Configuration**: Automatically generates patient-specific YAML configurations from templates
- **Resource Management**: Optional delays between submissions to manage cluster load
- **Comprehensive Logging**: Individual logs for each patient plus master orchestration log
- **Dry-Run Support**: Test configuration generation without submitting actual jobs
- **Flexible Input**: Supports various SSM file naming conventions

## Directory Structure

When you run the multi-patient pipeline, it creates the following structure:

```
your_output_directory/
├── configs/
│   └── generated/                    # Auto-generated patient configs
│       ├── patient1_config.yaml
│       ├── patient2_config.yaml
│       └── ...
├── patients/                         # Patient analysis results  
│   ├── patient1/
│   │   └── initial/                  # Standard pipeline structure
│   │       ├── bootstraps/
│   │       ├── aggregation_results/
│   │       ├── markers/
│   │       └── logs/
│   ├── patient2/
│   │   └── initial/
│   └── ...
└── logs/
    ├── multi_patient_master.log      # Overall orchestration log
    ├── patient1_submission.log       # Per-patient submission logs
    ├── patient2_submission.log
    └── ...
```

## Usage

### Basic Usage

```bash
bash multi_patient_pipeline.sh <ssm_directory> <config_template> <output_base_directory>
```

### Examples

```bash
# Process all SSM files using standard template
bash multi_patient_pipeline.sh data/patients/ configs/template_multi_patient.yaml /path/to/results/

# Test with reduced resources and 60-second delays between submissions
bash multi_patient_pipeline.sh data/patients/ configs/template_multi_patient_test.yaml /path/to/test_results/ --delay=60

# Dry run to test configuration generation without submitting jobs
bash multi_patient_pipeline.sh data/patients/ configs/template_multi_patient.yaml /path/to/results/ --dry-run
```

### Command Line Options

- `--delay=N`: Add N seconds delay between patient submissions (default: 0)
- `--dry-run`: Test configuration generation without submitting actual jobs
- `--help`: Show detailed usage information

## Input Requirements

### SSM Directory Structure

Your SSM files should be organized in a single directory:

```
data/patients/
├── CRUK0044.txt
├── CRUK0045.txt
├── CRUK0046_subset.txt
├── patient_001.ssm
└── patient_002.ssm
```

### SSM File Naming

The pipeline automatically extracts patient IDs from filenames:

- `CRUK0044.txt` → Patient ID: `CRUK0044`
- `CRUK0045_subset.txt` → Patient ID: `CRUK0045`
- `patient_001.ssm` → Patient ID: `patient_001`

Supported file extensions: `.txt`, `.ssm`

### SSM File Format

Each SSM file must follow the standard TracerX format:

```
id  gene                      a               d               mu_r  mu_v
s0  GPR65_14_88477948_G>T     1287,1116,1262  2135,1233,1325  0.999 0.499
s1  GLIPR1L1_12_75728663_C>A  238,351,385     359,367,406     0.999 0.499
...
```

## Configuration Templates

### Standard Template (`configs/template_multi_patient.yaml`)

Use this for production analyses:
- 100 bootstraps for robust statistics
- Standard resource allocation (8-16GB per step)
- Conservative SLURM array limits for multi-patient workloads

### Testing Template (`configs/template_multi_patient_test.yaml`)

Use this for development and validation:
- 20 bootstraps for faster testing
- Reduced resources (4-8GB per step)
- Debug logging enabled
- Shorter walltimes

### Template Placeholders

Templates use placeholders that are automatically replaced:

- `PLACEHOLDER_PATIENT_ID`: Replaced with extracted patient ID
- `PLACEHOLDER_SSM_FILE`: Replaced with absolute path to SSM file
- `PLACEHOLDER_OUTPUT_DIR`: Replaced with patient-specific output directory
- `PLACEHOLDER_CODE_DIR`: Replaced with TracerX-MP repository directory

## Resource Management and SLURM Array Considerations

### Array Job Limits

The pipeline handles SLURM array limits intelligently:

- Each patient gets separate array jobs for PhyloWGS stage
- Array jobs are independent (Patient A: `--array=0-99%20`, Patient B: `--array=0-99%20`)
- The 1001 array limit applies per job, not across all patients
- Conservative `parallel_limit` settings in templates prevent overwhelming the cluster

### Managing Cluster Load

**Option 1: Sequential Submission (Recommended)**
```bash
# Submit patients with 60-second delays
bash multi_patient_pipeline.sh data/ template.yaml output/ --delay=60
```

**Option 2: Batch Processing**
```bash
# Process first 5 patients
mkdir -p data/batch1/
cp data/patient_{1..5}.txt data/batch1/
bash multi_patient_pipeline.sh data/batch1/ template.yaml output_batch1/

# Process next 5 patients after first batch completes
mkdir -p data/batch2/
cp data/patient_{6..10}.txt data/batch2/
bash multi_patient_pipeline.sh data/batch2/ template.yaml output_batch2/
```

## Monitoring and Debugging

### Monitor All Jobs

```bash
# Check overall job status
squeue -u $USER

# Check specific patient jobs
squeue -u $USER | grep CRUK0044

# Monitor logs in real-time
tail -f /path/to/results/logs/multi_patient_master.log
tail -f /path/to/results/logs/patient1_submission.log
```

### Common Issues and Solutions

**Issue**: "No SSM files found"
```bash
# Check file extensions and permissions
ls -la data/patients/
# Ensure files have .txt or .ssm extension
```

**Issue**: "Pipeline submission failed"
```bash
# Check individual patient logs
cat /path/to/results/logs/patient1_submission.log
# Check master_pipeline.sh works for single patient
bash master_pipeline.sh /path/to/results/configs/generated/patient1_config.yaml --dry-run
```

**Issue**: Array jobs exceed cluster limits
```bash
# Use testing template with lower limits
bash multi_patient_pipeline.sh data/ configs/template_multi_patient_test.yaml output/
# Or add delays between submissions
bash multi_patient_pipeline.sh data/ template.yaml output/ --delay=120
```

## Testing Your Setup

### 1. Test with Dry Run

```bash
# Test configuration generation
bash multi_patient_pipeline.sh data/test/ configs/template_multi_patient_test.yaml test_output/ --dry-run
```

### 2. Test with Single Patient

```bash
# Create test directory with one SSM file
mkdir -p data/single_test/
cp data/patients/CRUK0044.txt data/single_test/
bash multi_patient_pipeline.sh data/single_test/ configs/template_multi_patient_test.yaml test_single/ --delay=0
```

### 3. Test with Small Batch

```bash
# Test with 2-3 patients using testing template
mkdir -p data/small_batch/
cp data/patients/CRUK004{4,5,6}.txt data/small_batch/
bash multi_patient_pipeline.sh data/small_batch/ configs/template_multi_patient_test.yaml test_batch/ --delay=30
```

## Performance Considerations

### Resource Planning

For N patients with standard template (100 bootstraps each):

- **Bootstrap stage**: N concurrent jobs, ~2 CPU + 8GB each
- **PhyloWGS stage**: N array jobs (100 elements each), ~2 CPU + 8GB per element
- **Aggregation stage**: N concurrent jobs, ~1 CPU + 16GB each  
- **Marker selection**: N concurrent jobs, ~1 CPU + 16GB each

### Timing Estimates

With standard template:
- **Bootstrap**: 30-60 minutes per patient
- **PhyloWGS**: 2-4 hours per patient (depends on parallel_limit)
- **Aggregation**: 15-30 minutes per patient
- **Marker selection**: 30-60 minutes per patient

**Total per patient**: ~4-6 hours
**Total for 10 patients**: ~40-60 hours (with proper parallelization)

### Optimization Tips

1. **Use testing template first**: Validate with reduced resources
2. **Batch processing**: Group patients to control resource usage
3. **Time submissions**: Use `--delay` during peak cluster hours
4. **Monitor resources**: Adjust template parameters based on cluster load

## Integration with Existing Workflows

### Converting Single-Patient Configs

If you have existing single-patient configs, you can create a template:

```bash
# Take an existing config and replace patient-specific values with placeholders
cp configs/existing_patient_config.yaml configs/my_template.yaml
# Edit my_template.yaml to replace:
# patient_id: "CRUK0044" → patient_id: "PLACEHOLDER_PATIENT_ID"
# ssm_file: "/path/to/CRUK0044.txt" → ssm_file: "PLACEHOLDER_SSM_FILE"
# base_dir: "/path/to/CRUK0044/" → base_dir: "PLACEHOLDER_OUTPUT_DIR"
```

### Using with Existing Data

```bash
# If you have existing SSM files in different locations
mkdir -p data/all_patients/
find /path/to/existing/data/ -name "*.txt" -exec cp {} data/all_patients/ \;
bash multi_patient_pipeline.sh data/all_patients/ configs/template_multi_patient.yaml new_analysis/
```

## Troubleshooting

### Configuration Issues

```bash
# Test template validity
python -c "import yaml; yaml.safe_load(open('configs/template_multi_patient.yaml'))"

# Test placeholder substitution
sed -e "s|PLACEHOLDER_PATIENT_ID|TEST_PATIENT|g" \
    -e "s|PLACEHOLDER_SSM_FILE|/test/path.txt|g" \
    -e "s|PLACEHOLDER_OUTPUT_DIR|/test/output|g" \
    -e "s|PLACEHOLDER_CODE_DIR|/test/code|g" \
    configs/template_multi_patient.yaml
```

### SLURM Issues

```bash
# Check queue limits
sinfo
sacctmgr show qos
scontrol show partition pool1

# Check job history
sacct -u $USER --starttime=today

# Cancel all jobs if needed
scancel -u $USER
```

## Next Steps

After running the multi-patient pipeline:

1. **Validate Results**: Check that each patient has complete output in their `initial/` directory
2. **Aggregate Analysis**: Use the individual results for comparative studies
3. **Longitudinal Analysis**: For patients with temporal data, proceed to Stage 5 (longitudinal analysis)

## Support

For issues or questions:

1. Check the individual patient logs in `/path/to/results/logs/`
2. Validate single-patient pipeline works with `master_pipeline.sh`
3. Test with dry-run mode and testing template
4. Review this documentation for common solutions 