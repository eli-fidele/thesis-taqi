#!/bin/bash

# Remove the existing files
rm -rf R/eigen.R
rm -rf R/ensemble.R
rm -rf R/matrices.R
rm -rf R/matrices-diagnostics.R
rm -rf R/analyze.R
rm -rf R/evolve.R
rm -rf R/simulate.R
rm -rf R/vis.R

# Copy them from developing environment
cp ../thesis-taqi/R/eigen.R R/eigen.R
cp ../thesis-taqi/R/ensemble.R R/ensemble.R
cp ../thesis-taqi/R/matrices.R R/matrices.R
cp ../thesis-taqi/R/matrices-diagnostics.R R/matrices-diagnostics.R
cp ../thesis-taqi/R/analyze.R R/analyze.R
cp ../thesis-taqi/R/evolve.R R/evolve.R
cp ../thesis-taqi/R/simulate.R R/simulate.R
cp ../thesis-taqi/R/vis.R R/vis.R

