
# Import dependencies
require(tidyverse)
require(patchwork)
require(matrixcalc)

# Imports all source code for a variable depth of directory.
.src <- function(depth = 1){
  # Get prefix based on directory depth for use in directories with variable depth
  path_prefix <- paste(c(rep("../",depth-1)), sep = "", collapse = "")
  dir_path <- paste(path_prefix,"../R/", sep = "", collapse = "")
  # Get script file names and get their path
  script_files <- dir(path = dir_path, pattern = "\\.R$")
  # Helper function to parse string name
  .ADDpref <- function(filename, dirpath){paste(dirpath, filename, sep = "", "")}
  script_paths <- as.list(purrr::map_chr(script_files, .ADDpref, dir_path))
  #print(script_paths)
  sapply(script_paths, source, .GlobalEnv)
}
