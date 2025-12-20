#!/bin/bash
#SBATCH --job-name=longitudinal
#SBATCH --output=longitudinal_%j.log
#SBATCH --error=longitudinal_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=pool1

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="$(dirname "$SCRIPT_DIR")"

# Run longitudinal analysis
conda run -n mase_phi_hpc python "${CODE_DIR}/src/longitudinal/longitudinal_main.py" --config "$1"