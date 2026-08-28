# PackLab contributor guide

## Project purpose

PackLab is a Python package for three-dimensional hard-sphere packings. It
offers two intentionally distinct scientific workflows:

- `PackLab.monte_carlo` generates explicit random sequential adsorption (RSA)
  configurations and their statistics.
- `PackLab.analytical` evaluates a Percus--Yevick hard-sphere-mixture reference.

Do not describe these methods as equivalent. RSA is an irreversible,
history-dependent packing process; Percus--Yevick is an analytical reference
model. Validation should match physical inputs and explain expected differences.

## Repository map

| Path | Purpose |
| --- | --- |
| `PackLab/monte_carlo/` | Python-facing RSA workflow and result helpers. |
| `PackLab/analytical/` | Python-facing Percus--Yevick workflow and grid helper. |
| `PackLab/samplers.*` | Native radius-distribution bindings. |
| `PackLab/scattering/` | Optional PyMieSim-based scattering workflow. |
| `PackLab/cpp/` | C++20 implementations and pybind11 interfaces. |
| `tests/` | Public-API, packing, sampler, analytical, and scattering tests. |
| `docs/source/` | Sphinx documentation. |
| `docs/examples/` | Executable gallery examples: `monte_carlo`, `analytical`, and `validation`. |
| `docs/manuscript/` | Longer LaTeX scientific manuscript draft and bibliography. |

Native extension modules are installed alongside their Python APIs within
`PackLab/` and its public subpackages. Keep this layout; do not introduce a
separate `binary` package or an `interface_*` naming layer.

## Local setup and verification

Use Python 3.10 or newer. Native builds require CMake, a C++20 compiler,
pybind11, and OpenMP. On macOS, install `libomp` first if CMake cannot find
OpenMP.

```bash
python -m pip install -e ".[testing]"
python -m pytest
```

Install documentation dependencies before building docs:

```bash
python -m pip install -e ".[documentation,testing]"
python -m sphinx -E -b html docs/source docs/build/html
```

Useful focused checks:

```bash
python -m pytest tests/test_monte_carlo.py
python -m pytest tests/test_analytical.py
python -m pytest tests/test_samplers.py
python -m pytest tests/test_scattering.py
```

After modifying C++ sources, rebuild the editable install before testing. Do
not hand-edit `PackLab/_version.py`; the build backend generates it from Git.

## API and numerical conventions

- Treat `PackLab` imports and documented names as the public contract.
- Make breaking renames cleanly: remove obsolete names rather than adding
  compatibility aliases or legacy wrappers.
- Keep Monte-Carlo and analytical APIs clearly separated. Shared utilities are
  acceptable only when they do not blur that distinction.
- Use `wavenumber`, never the historical short name `p`, in new public Python,
  C++, type-stub, test, and documentation code.
- Prefer `wavenumber="auto"` in beginner-facing analytical examples. If a
  manual grid is used, explain the resolution choice and preserve the warning
  for inadequate Fourier sampling.
- Inputs carrying physical dimensions should use `TypedUnit`; convert units at
  the boundary and keep native kernels operating on explicit SI values.
- Add NumPy-style docstrings to Python APIs and pybind11 docstrings plus useful
  `__repr__` output to exposed C++ classes.

## C++ and bindings

- CMake uses C++20 and builds extensions with pybind11.
- Put native implementation beside its binding interface under `PackLab/cpp/`.
- Register a new native target in its local `CMakeLists.txt`; make its install
  destination match the corresponding Python package.
- Keep binding validation explicit: check shapes, units, finite values, and
  physical ranges before entering a numerical kernel.
- Tests must exercise the public Python API, not only a private native module.

## Tests, examples, and documentation

- Add or update tests for every behaviour change. Numerical features require a
  convergence or accuracy test, not only a smoke test.
- Use deterministic seeds for Monte-Carlo tests. Avoid assertions that depend
  on an exact stochastic trajectory unless that trajectory is the behaviour
  under test.
- Keep gallery examples small enough to execute during a documentation build.
  Put each example in the appropriate `monte_carlo`, `analytical`, or
  `validation` category.
- Explain assumptions and limitations in docs, especially the distinction
  between RSA and the analytical model.
- Do not commit generated documentation, coverage reports, build directories,
  or Python cache files unless a release workflow explicitly requires them.

## Packaging and CI

- Runtime dependencies and build configuration live in `pyproject.toml`.
- cibuildwheel builds CPython 3.11--3.13 wheels for Linux, macOS, and Windows.
- Linux uses `manylinux_2_28` so current scientific Python dependencies can be
  installed as binary wheels during isolated wheel tests. Do not downgrade this
  baseline without checking the complete dependency and wheel policy.
- Wheel tests must run outside the source checkout, as configured in
  `pyproject.toml`, so they test the installed wheel rather than imported local
  files.

## Citations and manuscript

Keep `CITATION.cff` and `.zenodo.json` aligned with release metadata. The
scientific manuscript draft lives in `docs/manuscript/packlab.tex`; its
bibliography is `docs/manuscript/references.bib`.

Compile the draft locally with:

```bash
cd docs/manuscript
latexmk -pdf packlab.tex
```

JOSS submissions use `paper.md` and `paper.bib`, not LaTeX. Treat the LaTeX
file as the fuller scientific source. When replacing a figure placeholder,
include parameters, units, uncertainty or convergence information, and a
caption that explains the scientific claim rather than merely naming the plot.

## Change discipline

Keep changes focused and avoid unrelated formatting churn. Before handing off,
run the narrowest relevant tests and then the full suite when native or public
API behaviour changed. Report any dependency, compiler, platform, or numerical
assumption that prevents verification.
