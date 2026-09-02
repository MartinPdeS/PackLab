# Manuscript figures

`figure1_architecture.svg` is the editable source for Figure 1. The matching
PDF is committed because `pdflatex` includes PDF figures directly.

`create_figure2.py` creates Figure 2 from a fixed-seed PackLab simulation. It
outputs the editable SVG and the PDF included by the manuscript:

```bash
python create_figure2.py
```

`create_figure3.py` creates the Monte-Carlo/analytical validation and
wavenumber-grid convergence figure using the same interface:

```bash
python create_figure3.py
```

After editing the SVG, regenerate the PDF before compiling the manuscript:

```bash
rsvg-convert --format=pdf \
  --output figure1_architecture.pdf \
  figure1_architecture.svg
```
