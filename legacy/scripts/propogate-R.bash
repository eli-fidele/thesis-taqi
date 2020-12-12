#!/bin/bash

rm -R code/R/
mkdir code/R/
cp R/matrices.R code/R/
cp R/eigenvectors.R code/R/
cp R/eigenmetrics.R code/R/
cp R/simulate.R code/R/
cp R/animate.R code/R/
rm -R sim_EV/R/
mkdir sim_EV/R/
cp R/matrices.R sim_EV/R/
cp R/eigenvectors.R sim_EV/R/
cp R/eigenmetrics.R sim_EV/R/
cp R/simulate.R sim_EV/R/
cp R/animate.R sim_EV/R/
