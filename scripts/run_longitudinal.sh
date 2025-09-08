#!/bin/bash
#SBATCH --job-name=longitudinal
#SBATCH --output=longitudinal_%j.log
#SBATCH --error=longitudinal_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

# Run longitudinal analysis
conda run -n mase_phi_hpc python ../src/longitudinal/longitudinal_main.py --config "$1"