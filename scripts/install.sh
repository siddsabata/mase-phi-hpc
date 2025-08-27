#!/bin/bash

echo "=== MASE-PHI-HPC Setup Script ==="

# make logs directory
mkdir -p logs

# Function to create conda environment with mamba fallback
create_conda_env() {
    local env_file=$1
    local env_description=$2
    
    echo "Setting up $env_description..."
    
    # Try mamba first (much faster), fallback to conda
    if command -v mamba >/dev/null 2>&1; then
        echo "Using mamba for faster environment creation..."
        mamba env create -f "$env_file"
    else
        echo "Mamba not found, using conda..."
        conda env create -f "$env_file"
    fi
    
    if [ $? -eq 0 ]; then
        echo "✓ $env_description created successfully"
    else
        echo "✗ $env_description creation failed"
        exit 1
    fi
}

# Setup main conda environment (Python 3.11+ stages)
create_conda_env "environment.yml" "Main conda environment for Python 3.11+ stages"

# Create PhyloWGS conda environment (Python 2.7)
create_conda_env "src/phylowgs/environment.yml" "PhyloWGS conda environment"

# install phylowgs 
echo "Installing PhyloWGS software..."
cd src/phylowgs
git clone https://github.com/morrislab/phylowgs.git
cd phylowgs
echo "repo cloned; compiling phylowgs"
g++ -o mh.o -O3 mh.cpp util.cpp `gsl-config --cflags --libs`
cd ../../

echo "=== Setup Complete ==="
echo "✓ conda mase_phi_hpc environment ready for bootstrap, aggregation, markers, longitudinal"
echo "✓ conda phylowgs_env ready for PhyloWGS"
echo "✓ PhyloWGS software compiled"