# Clone Frequency Tracking Implementation

This document describes the implementation of clone frequency tracking and visualization in the longitudinal cancer evolution analysis pipeline.

## Overview

The clone frequency tracking system computes clone-level frequencies from tree structures and VAF data, implementing the root-correction approach without copy-number adjustments. This extends the existing pipeline to provide clone-level insights alongside mutation-level analysis.

## Features

### Core Functionality
- **Clone Frequency Computation**: Aggregates mutation VAFs by clone membership to compute clone frequencies
- **Root Correction**: Implements subtree remainder calculation to handle unrooted tree structures
- **Time Series Analysis**: Tracks clone frequencies across multiple timepoints
- **Visualization**: Generates trajectory plots and heatmaps for clone frequency evolution

### Output Files
- `clone_frequencies.csv`: Tidy DataFrame with clone frequencies over time
- `clone_frequency_summary.json`: Statistical summary of clone frequency data
- `*_clone_trajectories.png`: Line plots showing clone frequency evolution
- `*_clone_heatmap.png`: Heatmap visualization of clone frequencies

## Usage

### Configuration

Add the clone frequency tracking flag to your YAML configuration:

```yaml
parameters:
  n_markers: 2
  read_depth: 90000
  method: "phylowgs"
  lambda1: 0.0
  lambda2: 1.0
  focus_sample: 0
  track_clone_freq: true  # Enable clone frequency tracking
```

### Command Line

Enable clone frequency tracking via command line:

```bash
# Enable via configuration file
python longitudinal_main.py --config config.yaml

# Override configuration file setting
python longitudinal_main.py --config config.yaml --plot-clonal-freq
```

### Integration Points

The clone frequency tracking is integrated into both analysis modes:

1. **Fixed Analysis**: Computes clone frequencies using fixed marker VAFs across timepoints
2. **Dynamic Analysis**: Computes clone frequencies using selected marker VAFs from optimization

## Implementation Details

### Core Modules

#### `clone_frequency.py`
- `compute_clone_frequencies()`: Main function to compute clone frequencies from tree and VAF data
- `compute_subtree_remainder()`: Implements root correction with remainder pseudo-clone
- `validate_clone_frequency_data()`: Data validation and quality checks

#### `clone_visualizer.py`
- `plot_clone_trajectories()`: Creates line plots of clone frequency evolution
- `plot_clone_heatmap()`: Generates heatmap visualizations
- `create_clone_frequency_summary()`: Computes statistical summaries

#### `output_manager.py`
- `write_clone_frequencies()`: Standardized CSV output for clone frequency data

### Algorithm

1. **Clone Frequency Calculation**:
   - For each clone C and timepoint t:
   - F_raw(C,t) = mean(VAF_m,t) for all mutations m in clone C

2. **Root Correction**:
   - F_subtree(t) = sum of all clone frequencies at timepoint t
   - F_leafsum(t) = sum of leaf clone frequencies at timepoint t
   - F_remainder(t) = F_subtree(t) - F_leafsum(t)
   - Add remainder as pseudo-clone "unrooted_remainder"

3. **Visualization**:
   - Regular clones: colored lines with markers
   - Remainder pseudo-clone: grey dashed line with reduced opacity

## Data Format

### Input
- Tree distribution with converged/final tree structure
- VAF data from ddPCR measurements across timepoints

### Output
Clone frequencies CSV format:
```csv
sample,time,clone_id,freq
sample_1,2021-01-15,0,0.245
sample_1,2021-01-15,1,0.123
sample_1,2021-01-15,unrooted_remainder,0.089
sample_1,2021-03-15,0,0.267
...
```

## File Structure

```
output_directory/
├── fixed_marker_analysis/          # Fixed analysis results
│   ├── clone_frequencies.csv       # Clone frequency data
│   ├── clone_frequency_summary.json # Statistical summary
│   └── visualizations/
│       ├── *_clone_trajectories.png
│       └── *_clone_heatmap.png
├── dynamic_marker_analysis/        # Dynamic analysis results
│   ├── clone_frequencies.csv
│   ├── clone_frequency_summary.json
│   └── visualizations/
│       ├── *_clone_trajectories.png
│       └── *_clone_heatmap.png
└── final_report.json              # Updated to include clone frequency files
```

## Testing

The implementation includes validation functions and error handling:

- Data format validation
- Tree structure consistency checks
- VAF data quality validation
- Graceful handling of missing data

## Integration with Existing Pipeline

The clone frequency tracking is:
- **Optional**: Controlled by configuration flag
- **Non-intrusive**: Doesn't modify existing analysis logic
- **Consistent**: Follows existing code patterns and styles
- **Modular**: Contained in separate modules for easy maintenance

## Performance Considerations

- Computation is performed after tree convergence (minimal overhead)
- VAF data is aggregated efficiently using pandas operations
- Memory usage scales with number of clones × timepoints
- Visualization generation is optional and can be disabled

## Future Enhancements

Potential extensions:
- Multi-sample support (currently single sample per analysis)
- Statistical testing for clone frequency changes
- Integration with phylogenetic uncertainty quantification
- Support for different aggregation methods beyond mean VAF

## Troubleshooting

Common issues and solutions:

1. **No VAF data**: Check that mutations in tree are present in ddPCR measurements
2. **Empty clone frequencies**: Verify tree structure and mutation assignments
3. **Visualization errors**: Ensure required packages (matplotlib, seaborn) are installed
4. **Performance issues**: Consider reducing number of timepoints or clones for large datasets

## References

- Implementation follows patterns from existing `longitudinal_visualizer.py`
- Tree structure handling based on `analyze.py` and `tree_updater.py`
- Output management consistent with `output_manager.py`
- Configuration handling follows `config_handler.py` patterns