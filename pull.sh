#!/bin/bash

# Remove the existing files
rm -rf R/eigen.R
rm -rf R/ensemble.R
rm -rf R/matrices.R
rm -rf R/matrices-diagnostics.R

# Copy them from developing environment
cp ../RMAT/R/eigen.R R/eigen.R
cp ../RMAT/R/ensemble.R R/ensemble.R
cp ../RMAT/R/matrices.R R/matrices.R
cp ../RMAT/R/matrices-diagnostics.R R/matrices-diagnostics.R

