# PackLab manuscripts

The manuscript sources are separated by target journal:

- `softwareX/` contains the full Original Software Publication prepared with
  Elsevier's single-column `elsarticle` format and SoftwareX section structure.
- `joss/` is reserved for the shorter JOSS `paper.md` and `paper.bib`
  submission, which will be prepared separately if that venue is selected.

Build the SoftwareX draft locally with:

```bash
cd docs/manuscript/softwareX
latexmk -pdf packlab.tex
```

See `softwareX/SUBMISSION_CHECKLIST.md` for the journal-format checks and the
author declarations that still require human confirmation.
