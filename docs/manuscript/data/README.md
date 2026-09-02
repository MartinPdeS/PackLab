# Digitised reference data

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
numerical dataset. They are intended for visual comparison in manuscript
Figure 4, so their plotting uncertainty and digitisation precision should not
support a high-precision numerical error claim.
