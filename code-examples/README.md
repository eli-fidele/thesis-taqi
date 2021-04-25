# Code Examples: A reproducible workflow
Contained within this directory is a set of `.Rmd` files in which lies the code for which every code chunk in this thesis can be reproduced.

## Setup:
The setup is as follows.
- In any `.Rmd` file, there will be one code example per code chunk. See the code chunk breakdown in the subsection below.

## Code Chunk Structure:
- In each code chunk, there is a `set.seed()` argument to ensure reproducibility.
- After setting the seed, the relevant code chunk may be ran. Then, its output is copied manually and it is adjusted to be inserted in a `lstliting` LaTeX environment to prepare the code chunk for output.
- After editing the output till satisfaction, directly copy and paste the code chunks into the `LaTeX` file.
