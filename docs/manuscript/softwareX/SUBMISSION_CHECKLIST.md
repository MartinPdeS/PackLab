# SoftwareX submission checklist

This draft follows the SoftwareX Original Software Publication structure and
Elsevier's single-column article format. Confirm the journal's online Guide for
Authors again immediately before submission because editorial requirements can
change.

## Implemented in the manuscript

- [x] Software name followed by a descriptive title.
- [x] Approximately 100-word abstract and no more than six keywords.
- [x] Single-column `elsarticle` layout with line numbers.
- [x] Code and executable-software metadata tables.
- [x] Sections for Motivation and significance, Software description,
  Illustrative examples, Impact, and Conclusions.
- [x] At most six figures.
- [x] Reproducible examples, empirical validation, limitations, and references.
- [x] Data-availability statement.
- [x] Generative-AI disclosure immediately before the references.
- [x] Separate highlights file with four concise highlights.
- [x] MIT-licensed public repository and version-specific source link.

## Author confirmation required before submission

- [ ] Confirm the corresponding author and preferred institutional email.
- [ ] Confirm the full postal affiliation and spelling for every author.
- [ ] Confirm the drafted CRediT roles with all authors.
- [ ] Confirm the drafted no-funding statement with all authors.
- [ ] Confirm the drafted no-competing-interest declaration with all authors.
- [ ] Confirm that the AI disclosure accurately describes the authors' use and
  oversight of OpenAI Codex.
- [ ] Confirm that DOI `10.5281/zenodo.20207810` is the correct archive for the
  submitted version, or create a new version-specific Zenodo deposit.
- [ ] Confirm that the v0.7.0 GitHub tag contains exactly the software evaluated
  in the article.
- [ ] Add any research projects, external users, citations, or adoption evidence
  that can substantiate the Impact section.

## Final packaging checks

- [ ] Keep the main text at or below 3000 words, including the abstract,
  captions, and footnotes but excluding metadata tables and references.
- [ ] Build from a clean checkout with `latexmk -pdf packlab.tex`.
- [ ] Inspect every figure at final size and supply editable/high-resolution
  source files.
- [ ] Verify every URL, DOI, software version, and numerical claim.
- [ ] Upload the manuscript PDF and complete LaTeX source package.
- [ ] Upload `highlights.txt` separately if requested by the submission system.
- [ ] Consider using `figures/architecture.png` as an optional graphical
  abstract after checking the current size requirements.
