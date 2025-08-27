#!/bin/bash

echo "=== MASE-PHI-HPC Setup Script ==="

# make logs directory
mkdir -p logs

# Setup uv environment (Python 3.11+ stages)
echo "Setting up uv environment for Python 3.11+ stages..."
uv sync
if [ $? -eq 0 ]; then
    echo "✓ uv environment created successfully"
else
    echo "✗ uv environment setup failed"
    exit 1
fi

# Install packages not available in uv registry
echo "Installing additional packages via pip (gurobipy==11.0.2, python-graphviz)..."
uv pip install gurobipy==11.0.2 python-graphviz
if [ $? -eq 0 ]; then
    echo "✓ Additional packages installed successfully"
else
    echo "✗ Package installation failed"
    echo "  Note: Gurobi requires a license for optimization. Academic licenses available at gurobi.com"
fi

# Create PhyloWGS conda environment (Python 2.7)
echo "Creating PhyloWGS conda environment..."
conda env create -f src/phylowgs/environment.yml
if [ $? -eq 0 ]; then
    echo "✓ PhyloWGS conda environment created successfully"
else
    echo "✗ PhyloWGS conda environment creation failed"
    exit 1
fi

# install phylowgs 
echo "Installing PhyloWGS software..."
cd src/phylowgs
git clone https://github.com/morrislab/phylowgs.git
cd phylowgs
echo "repo cloned; compiling phylowgs"
g++ -o mh.o -O3 mh.cpp util.cpp `gsl-config --cflags --libs`
cd ../../

echo "=== Setup Complete ==="
echo "✓ uv environment ready for bootstrap, aggregation, markers, longitudinal"
echo "✓ conda phylowgs_env ready for PhyloWGS"
echo "✓ PhyloWGS software compiled"