# Contributing to PackLab

Thank you for helping improve PackLab. Bug reports, documentation corrections,
numerical validation cases, and focused code contributions are welcome.

## Before opening a change

For substantial features or changes to scientific conventions, open an issue
first so the scope and validation strategy can be discussed. Please describe
the physical model, expected behavior, units, and relevant references. Security
issues should not be reported publicly; contact the maintainer through the
email listed in `pyproject.toml`.

PackLab deliberately separates two kinds of workflow:

- `PackLab.monte_carlo` generates explicit RSA configurations and samples
  fixed-volume hard-sphere configurations with Metropolis Monte Carlo.
- `PackLab.analytical` evaluates the Percus--Yevick equilibrium reference.

Contributions must preserve this distinction. Comparisons should use matched
physical inputs and explain expected model differences.

## Development setup

PackLab requires Python 3.10 or newer, CMake, a C++20 compiler, pybind11, and
OpenMP. On macOS, install `libomp` if CMake cannot find OpenMP.

```bash
git clone https://github.com/MartinPdeS/PackLab.git
cd PackLab
make setup
make doctor
```

After changing C++ code, rebuild the editable installation before testing.
Do not edit `PackLab/_version.py`; it is generated from Git metadata.

## Contribution expectations

- Add or update tests for behavior changes. Numerical features need an accuracy
  or convergence test, not only a smoke test.
- Test public Python APIs and use deterministic seeds for stochastic tests.
- Use `TypedUnit` for dimensioned public inputs and `wavenumber` for the public
  Fourier-space variable name.
- Add NumPy-style Python docstrings and useful pybind11 docstrings for native
  APIs. Keep gallery examples small enough to run during a documentation build.
- Update documentation when behavior, assumptions, limitations, or public APIs
  change. Do not commit generated documentation, build products, or caches.
- Keep pull requests focused and avoid unrelated formatting changes.

Run `make help` to see the supported workflows. Run the relevant focused test
first, followed by the full local check when changing native code or public
behavior:

```bash
python -m pytest tests/test_analytical.py
python -m pytest tests/test_monte_carlo.py
python -m pytest tests/test_samplers.py
python -m pytest tests/test_scattering.py
make check
```

Documentation and the SoftwareX manuscript can be checked independently with:

```bash
make docs
make manuscript-check
```

`make reproduce-paper` regenerates the quantitative manuscript figures, checks
the stored independent validation data, and rebuilds the article. It is slower
than the normal pre-push workflow and is intended for changes affecting
scientific results or manuscript figures.

## Pull requests

Explain what changed, why it is needed, and how it was verified. Include the
platform and compiler for native-code changes, and report numerical assumptions
or known limitations. By submitting a contribution, you agree that it may be
distributed under PackLab's MIT License.
