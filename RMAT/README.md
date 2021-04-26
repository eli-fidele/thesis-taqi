# RMAT: Random Matrix Analysis Toolkit

RMAT is an R package for simulating random matrices and ensembles as well as analyzing their eigenvalue spectra and dispersions. 

## Installation

NOTE: This package has been submitted to CRAN and is pending release. 

You can install the package in one of two ways:

``` r
# Through CRAN
install.packages("RMAT")
# Through Github
devtools::install_github(repo = "ataqi23/RMAT")
```

## Usage

The package can be thought to contain two priamry modules: the *matrix module* and *spectral statistics module*.

### Matrices

Included in the package are various wrapper functions for generating random matrices very efficiently. They follow the ``RM_xxx`` format in allusion to the `stats` sampling functions like `rnorm`. Every matrix sampling function has an ensemble counterpart in the ``RME_xxx`` format. 

Consider the following example:

``` r
library(RMAT)
# Generate an ensemble of 50x50 Complex Hermitian Standard Normal matrices
ens <- RME_norm(N = 20, cplx = T, herm = T, size = 50)
```
The matrices included are:

- Uniform Matrices
- Normal Matrices
- Tridiagonal Matrices
- Beta Ensemble for b = 1,2,4
- Generalized Beta Ensemble for b > 0
- Stochastic Matrices
- Erdos-Renyi Graph Transition Matrices

### Spectral Statistics 

The two primary functions in the spectral statistics module are ``spectrum`` and ``dispersion``. These two functions take either a matrix or an ensemble, and return the respective eigenvalue spectrum or dispersion. Consider the following example:

``` r
library(RMAT)
# Generate a random matrix ensemble
ens <- RME_beta(N = 20, beta = 4, size = 50)
# Compute the spectral statistics of the ensemble
ens_spectrum <- spectrum(ens, sort_norms = FALSE)
ens_dispersion <- dispersion(ens, pairs = "consecutive", sort_norms = FALSE)
```

## About

This package was developed as a part of my Reed senior thesis in the Spring of 2021.



