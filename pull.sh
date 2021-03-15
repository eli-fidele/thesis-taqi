#!/bin/bash

# Remove the existing files
rm -rf R/spectrum.R
rm -rf R/dispersion.R
rm -rf R/matrices.R
rm -rf R/matrices-diagnostics.R

# Copy them from developing environment
cp ../RMAT/R/spectrum.R R/spectrum.R
cp ../RMAT/R/dispersion.R R/dispersion.R
cp ../RMAT/R/matrices.R R/matrices.R
cp ../RMAT/R/mat-diag.R R/mat-diag.R


