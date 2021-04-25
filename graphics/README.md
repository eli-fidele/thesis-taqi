# Graphics: A reproducible workflow
Contained within this directory is a set of `.Rmd` files in which lies the code for which every plot in this thesis can be reproduced. For graphics with more involved code, they may have their own independent `.Rmd` files such as `wigner.Rmd`.

## Setup:
The setup is as follows.
- In any `.Rmd` file, there will be one plot per code chunk. 
- Running the code chunk should automatically generate the graphic and place it in the appopriate directory.
- After generating the plot, running `pdflatex` should update the graphics! The file reading in the filenames is in `teX/graphics.tex`.

## Code Chunk Structure:
- In each code chunk, there is a `set.seed()` argument to ensure reproducibility.
- There are also manual `run` and `save` booleans that must manually be set to `TRUE` to run the code and save the plot respectively. They are implemented for safety.
- After running the code, the last plot saved will be captured by `ggsave()` and be appopriately written into the correct directory to be read by `LaTeX`.
