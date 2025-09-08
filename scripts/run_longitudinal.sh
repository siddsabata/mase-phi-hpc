#!/bin/bash
#SBATCH --job-name=longitudinal
#SBATCH --output=/home/ssabata/masephi-test/cruk0044/cruk0044/longitudinal/longitudinal_%j.log
#SBATCH --error=/home/ssabata/masephi-test/cruk0044/cruk0044/longitudinal/longitudinal_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=pool1

# Run longitudinal analysis
conda run -n mase_phi_hpc python /home/ssabata/mase-phi-hpc/src/longitudinal/longitudinal_main.py --config "$1"