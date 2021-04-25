# The LaTeX document

The actual thesis document was written in LaTeX using the standard Reed Thesis template. This was a large project, involving many moving parts like constantly changing code examples, graphics, etc. As such, proper refactoring paid off as it streamlined the process of stitching this document together. As such, here are a few pointers.

## Structure:
- Master document: the master `.tex` file is the aptly named `thesis.tex` file. It sources and sets up the Reed template and inputs all the other `.tex` files. By imposing this structure of sourcing smaller documents, we are able to have a tidier thesis.
- Chapter preview: Plus, if we wanted to compile only a single chapter, we could do so by compiling the stripped down `chapter.tex` file, reducing overall developement time! This is one payoff of writing each chapter in its own `chapter#.tex` file.
- Preamble: the preamble of any `LaTeX` document is very critical; it sources packages, sets up environments and macros, and so much more. So, having a well-factored preamble was necessary. So, the master preamble file is called `preamble.tex` and can be found in this directory. The bulk of the preamble, however, was factored as a directory in which there are smaller source files that comprise it. A breakdown is provided below.

## Preamble:
- Macros: all macroes and commands can be found in `commands.tex`
- Graphics: all graphics wrapper commands can be found in `graphics.tex`
- Proof: the proof in the appendix can be found in `proofs.tex`
- Tables: the table wrapper commands can be found in `tables.tex`
