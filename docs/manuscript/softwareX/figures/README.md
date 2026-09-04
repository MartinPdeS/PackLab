# Manuscript figures

`figure1_architecture.svg` is the editable source for Figure 1. The matching
PDF is committed because `pdflatex` includes PDF figures directly.

`create_figure2.py` creates Figure 2 from a fixed-seed PackLab simulation. It
outputs the editable SVG and the PDF included by the manuscript:

```bash
python create_figure2.py
```

`create_figure3.py` creates the rendered Figure 4: matched RSA--PY and
Metropolis--PY comparisons in the top row, followed by the wavenumber-grid
convergence panel. The Metropolis comparison uses 12 independently initialized
chains, a 600-sweep burn-in, and a 300--600-sweep stability check:

```bash
python create_figure3.py
```

`create_figure4_independent_validation.py` creates the real- and
reciprocal-space independent-validation figure. Its lower panel uses the same
PackLab--mctpy comparison code as `development/compare_mctpy.py`. The mctpy
results are read from the versioned reference CSV, so generating the manuscript
figure does not require mctpy:

```bash
python create_figure4_independent_validation.py
```

Only maintainers updating the reference dataset need mctpy; its provenance and
regeneration command are documented in `../data/README.md`.

After editing the SVG, regenerate the PDF before compiling the manuscript:

```bash
rsvg-convert --format=pdf \
  --output figure1_architecture.pdf \
  figure1_architecture.svg
```
