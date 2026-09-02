# PackLab manuscript draft

`packlab.tex` is the scientific manuscript draft and `references.bib` contains
its bibliography. Build it locally with:

```bash
cd docs/manuscript
latexmk -pdf packlab.tex
```

JOSS submissions are made as `paper.md` and `paper.bib` in the repository
root, using the [JOSS paper template](https://joss.readthedocs.io/en/latest/submitting.html).
This LaTeX draft is intentionally separate from the final JOSS submission
format so the fuller manuscript can evolve without constraining that workflow.
