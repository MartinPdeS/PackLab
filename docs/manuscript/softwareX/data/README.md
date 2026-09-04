# Reference data

`tsang_2001_figure_8_3_4_digitized.csv` contains 254 manually digitised marker
locations from Figure 8.3.4 of Tsang, Kong, Ding, and Ao (2001), cited in the
manuscript as `Tsang2001`.

`tsang_2001_figure_8_1_3_digitized.csv` contains 211 manually digitised
monodisperse \gls{py} curve locations from Figure 8.1.3 at volume fractions
0.2 and 0.3. Its columns are `reduced_separation` (the book coordinate
\(r/b\)), `g_r`, and `volume_fraction`.

Columns are:

- `reduced_separation`: the book coordinate \(r/(2a_1)\);
- `g_ij`: the reported partial pair correlation;
- `component_pair`: one of `g11`, `g12`, or `g22`.

The values are digitised from a printed figure, not the authors' original
numerical dataset. They are intended for visual comparison in the manuscript's
independent-validation figure, so their plotting uncertainty and digitisation
precision should not support a high-precision numerical error claim.

## mctpy structure-factor reference

`mctpy_0_0_2_structure_factors.csv` contains partial structure factors produced
with `mctpy` 0.0.2.post1. It uses the density-scaled convention

\[
S_{ij}(k)=\delta_{ij}+\sqrt{\rho_i\rho_j}\,\widetilde h_{ij}(k),
\]

so the structure-factor matrix approaches the identity at high wavenumber. The
columns are:

- `mixture`: `monodisperse`, `binary`, or `ternary`;
- `wavenumber_um_inverse`: wavenumber in $\mathrm{\mu m}^{-1}$;
- `row`, `column`: zero-based species indices for a unique matrix component;
- `S_ij`: the dimensionless partial structure factor.

All cases use 601 points over $0\leq k\leq30~\mathrm{\mu m}^{-1}$. Their
physical inputs are defined in `development/compare_mctpy.py`. The data were
generated from `mctpy.structurefactors.hsmPY`, with number densities expressed
in $\mathrm{\mu m}^{-3}$ and diameters in $\mathrm{\mu m}$. The source release
is archived at <https://doi.org/10.5281/zenodo.19799936>.

Normal PackLab validation reads the checked-in CSV and does not require mctpy.
To deliberately regenerate the reference:

```bash
python -m pip install mctpy==0.0.2.post1 numba scipy h5py
python development/extract_mctpy_reference.py
```

Review changes to the CSV before committing them. Updating the reference
implementation or mixture definitions requires updating this provenance note
and the manuscript's reported errors.
