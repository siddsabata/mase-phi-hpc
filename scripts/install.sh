#!/bin/bash

# make logs directory
mkdir -p logs

# install phylowgs 
cd src/phylowgs
git clone https://github.com/morrislab/phylowgs.git
cd phylowgs
echo "repo cloned; compiling phylowgs"
g++ -o mh.o -O3 mh.cpp util.cpp `gsl-config --cflags --libs`
cd ../../

echo "phylowgs installed"