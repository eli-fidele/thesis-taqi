
# Import dependencies
require(tidyverse)
require(patchwork)
require(matrixcalc)

# Imports all source code for a variable depth of directory.
.src <- function(depth = 1){
  # Get directory path based on directory depth for use in directories with given variable depth
  dir_path <- paste(c(rep("../",depth-1),"../R/"), sep = "", collapse = "")
  # Get script file names and get their path
  script_files <- dir(path = dir_path, pattern = "\\.R$")
  # Helper function to parse string name
  .ADDpref <- function(filename, dirpath){paste(dirpath, filename, sep = "", collapse = "")}
  script_paths <- as.list(purrr::map_chr(script_files, .ADDpref, dir_path))
  sapply(script_paths, source, .GlobalEnv)
}

# Imports a preset of source code for a variable depth of directory.
.srcp <- function(depth = 1, preset){
  # Get directory path based on directory depth for use in directories with given variable depth
  dir_path <- paste(c(rep("../",depth-1),"../R/"), sep = "", collapse = "")
  # Script file presets
  preset1 <- c("eigen.R", "ensemble.R", "matrices.R", "mat-diag.R")
  preset2 <- c(preset1, "analyze.R", "evolve.R", "simulate.R", "vis.R")
  # Set preset
  if(preset == 1){script_files <- preset1}
  if(preset == 2){script_files <- preset2}
  # Helper function to parse string name
  .ADDpref <- function(filename, dirpath){paste(dirpath, filename, sep = "", "")}
  script_paths <- as.list(purrr::map_chr(script_files, .ADDpref, dir_path))
  sapply(script_paths, source, .GlobalEnv)
}