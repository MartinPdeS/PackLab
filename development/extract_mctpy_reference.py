"""Regenerate the versioned mctpy structure-factor reference CSV.

This provenance utility is only needed when deliberately updating the stored
reference. Normal comparisons and manuscript figure generation do not import
or install mctpy.

Install the reference implementation and its currently undeclared runtime
dependencies before executing this script::

    python -m pip install mctpy==0.0.2.post1 numba scipy h5py
    python development/extract_mctpy_reference.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import mctpy
import numpy as np
from mctpy.structurefactors import hsmPY

from compare_mctpy import MIXTURES, REFERENCE_DATA_PATH, partial_number_densities


EXPECTED_MCTPY_VERSION = "0.0.2.post1"


def write_reference(output_path: Path = REFERENCE_DATA_PATH) -> None:
    """Evaluate mctpy and write all unique partial structure factors."""
    if mctpy.__version__ != EXPECTED_MCTPY_VERSION:
        raise RuntimeError(
            f"Expected mctpy {EXPECTED_MCTPY_VERSION}, found {mctpy.__version__}. "
            "Review and update the provenance before changing reference versions."
        )
    wavenumber = np.linspace(0.0, 30.0, 601)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.writer(output_file, lineterminator="\n")
        writer.writerow(("mixture", "wavenumber_um_inverse", "row", "column", "S_ij"))
        for mixture in MIXTURES:
            densities = partial_number_densities(mixture)
            model = hsmPY(densities, 2.0 * np.asarray(mixture.radii_um))
            structure_factor, _ = model.Sq(wavenumber)
            for row in range(len(mixture.radii_um)):
                for column in range(row, len(mixture.radii_um)):
                    for wavenumber_value, structure_factor_value in zip(
                        wavenumber, structure_factor[:, row, column]
                    ):
                        writer.writerow(
                            (
                                mixture.name,
                                f"{wavenumber_value:.17g}",
                                row,
                                column,
                                f"{structure_factor_value:.17g}",
                            )
                        )

    print(f"Wrote mctpy reference data to {output_path}")


if __name__ == "__main__":
    write_reference()
