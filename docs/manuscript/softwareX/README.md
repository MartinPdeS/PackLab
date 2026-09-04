# SoftwareX manuscript

This directory contains the PackLab Original Software Publication prepared for
SoftwareX. The source uses Elsevier's `elsarticle` class in single-column
preprint mode and includes line numbers, highlights, the required software
metadata tables, reproducible validation, and the standard SoftwareX section
structure.

Build the article from this directory:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error packlab.tex
```

The quantitative inputs are stored in `data/`, while `figures/` contains both
the rendered figures and their generation scripts. The checked-in mctpy CSV is
a versioned validation reference; building the manuscript does not install or
execute mctpy.

Before submission, complete every open item in `SUBMISSION_CHECKLIST.md`, then
build and inspect the PDF from a clean checkout. Generated manuscript build
files are intentionally ignored by Git.
