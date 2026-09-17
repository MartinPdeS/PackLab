"""Portable, validated persistence for Monte-Carlo packing configurations."""


import json
from pathlib import Path
from typing import Any

import numpy as np

from PackLab.monte_carlo.domain import PackingDomain
from PackLab.monte_carlo.simulator import PackingConfiguration
from PackLab.units import ureg

_SCHEMA = "packlab.monte_carlo.packing"
_VERSION = 1
_SI_UNITS = {"positions_m": "meter", "radii_m": "meter", "box_lengths_m": "meter"}
_REQUIRED_ARCHIVE_KEYS = {
    "positions_m",
    "radii_m",
    "classes_index",
    "box_lengths_m",
    "periodic",
    "metadata_json",
}


class LoadedPackingConfiguration:
    """
    A validated packing configuration reconstructed from a portable archive.

    Parameters
    ----------
    sphere_configuration : PackingConfiguration
        Native configuration reconstructed from SI arrays.
    domain : PackingDomain
        Native domain reconstructed from SI box dimensions.
    metadata : dict
        Archive metadata, including format version and source when known.

    Notes
    -----
    This object deliberately does not manufacture simulation statistics. Its
    ``sphere_configuration`` can be passed directly to
    :class:`MetropolisSimulator`.
    """

    def __init__(
        self,
        sphere_configuration: PackingConfiguration,
        domain: PackingDomain,
        metadata: dict[str, Any],
    ):
        self.sphere_configuration = sphere_configuration
        self.domain = domain
        self.metadata = metadata
        self.source = metadata["source"]
        self.run_metadata = dict(metadata.get("run_metadata", {}))

    @property
    def positions(self) -> Any:
        """Sphere centres as a quantity with metre units."""
        return self.sphere_configuration.positions

    @property
    def radii(self) -> Any:
        """Sphere radii as a quantity with metre units."""
        return self.sphere_configuration.radii

    @property
    def classes_index(self) -> np.ndarray:
        """Non-negative particle class labels."""
        return np.asarray(self.sphere_configuration.classes_index)

    def compute_partial_pair_correlation_function(
        self, **kwargs: Any
    ) -> tuple[Any, np.ndarray]:
        """
        Compute partial pair correlations for this loaded configuration.

        Parameters
        ----------
        **kwargs : dict
            ``n_bins`` and optional ``maximum_pairs`` forwarded to the native
            configuration method.
        """
        centers, values = self.sphere_configuration.compute_partial_pair_correlation_function(
            self.domain, **kwargs
        )
        return np.asarray(centers) * ureg.meter, np.asarray(values)

    def save(self, path: str | Path, **kwargs: Any) -> None:
        """Save this loaded configuration to another portable archive."""
        save_packing(self, path, **kwargs)


def _quantity_values(value: Any, unit: str, name: str) -> np.ndarray:
    try:
        values = np.asarray(value.to(unit).magnitude, dtype=float)
    except AttributeError as error:
        raise TypeError(
            f"{name} must be a unit-bearing quantity convertible to {unit}."
        ) from error
    return values


def _extract_packing(
    packing: Any, domain: PackingDomain | None
) -> tuple[Any, PackingDomain, str, dict[str, Any]]:
    if isinstance(packing, LoadedPackingConfiguration):
        return packing.sphere_configuration, packing.domain, packing.source, packing.run_metadata
    if hasattr(packing, "sphere_configuration") and hasattr(packing, "domain"):
        return (
            packing.sphere_configuration,
            packing.domain,
            getattr(packing, "source", "unknown"),
            dict(getattr(packing, "run_metadata", {})),
        )
    if isinstance(packing, PackingConfiguration):
        if domain is None:
            raise ValueError("domain is required when saving a PackingConfiguration directly.")
        return packing, domain, "unknown", {}
    raise TypeError(
        "packing must be a PackingResult, LoadedPackingConfiguration, or PackingConfiguration."
    )


def _validate_arrays(
    positions: np.ndarray,
    radii: np.ndarray,
    classes: np.ndarray,
    lengths: np.ndarray,
    periodic: bool,
) -> None:
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("archive positions_m must have shape (N, 3).")
    if radii.ndim != 1 or radii.shape != (positions.shape[0],):
        raise ValueError("archive radii_m must have shape (N,) matching positions_m.")
    if classes.ndim != 1 or classes.shape != (positions.shape[0],):
        raise ValueError("archive classes_index must have shape (N,) matching positions_m.")
    if lengths.shape != (3,):
        raise ValueError("archive box_lengths_m must have shape (3,).")
    if not np.all(np.isfinite(positions)):
        raise ValueError("archive positions_m must contain only finite values.")
    if not np.all(np.isfinite(radii)) or np.any(radii <= 0.0):
        raise ValueError("archive radii_m must be finite and positive.")
    if not np.all(np.isfinite(lengths)) or np.any(lengths <= 0.0):
        raise ValueError("archive box_lengths_m must be finite and positive.")
    if not np.issubdtype(classes.dtype, np.integer) or np.any(classes < 0):
        raise ValueError("archive classes_index must contain non-negative integers.")
    if periodic:
        if np.any(positions < 0.0) or np.any(positions >= lengths):
            raise ValueError("archive periodic positions must lie in [0, box length).")
    elif np.any(positions < radii[:, None]) or np.any(positions > lengths - radii[:, None]):
        raise ValueError("archive contains a sphere outside the non-periodic domain.")

    for first in range(positions.shape[0]):
        deltas = positions[first + 1 :] - positions[first]
        if periodic and deltas.size:
            deltas -= lengths * np.round(deltas / lengths)
        distances_squared = np.einsum("ij,ij->i", deltas, deltas)
        minimum_squared = (radii[first] + radii[first + 1 :]) ** 2
        if np.any(distances_squared < minimum_squared):
            raise ValueError("archive contains overlapping spheres.")


def save_packing(
    packing: Any,
    path: str | Path,
    *,
    domain: PackingDomain | None = None,
    source: str | None = None,
    metadata: dict[str, Any] | None = None,
) -> None:
    """
    Save a packing configuration as a compressed, portable ``.npz`` archive.

    Parameters
    ----------
    packing : PackingResult, LoadedPackingConfiguration, or PackingConfiguration
        Configuration to persist. Supplying a bare configuration requires
        ``domain``.
    path : str or pathlib.Path
        Destination path, which must end in ``.npz``.
    domain : PackingDomain, optional
        Domain for a bare :class:`PackingConfiguration`.
    source : {"rsa", "metropolis", "unknown"}, optional
        Override the workflow source stored in metadata.
    metadata : dict, optional
        JSON-serializable caller metadata. It is stored alongside PackLab's
        schema, version, source, and available run options.
    """
    target = Path(path)
    if target.suffix != ".npz":
        raise ValueError("path must end in '.npz'.")
    configuration, packing_domain, inferred_source, run_metadata = _extract_packing(packing, domain)
    source = inferred_source if source is None else source
    if source not in {"rsa", "metropolis", "unknown"}:
        raise ValueError("source must be 'rsa', 'metropolis', or 'unknown'.")

    positions = _quantity_values(configuration.positions, "meter", "positions")
    radii = _quantity_values(configuration.radii, "meter", "radii")
    classes = np.asarray(configuration.classes_index)
    lengths = np.asarray(
        [
            packing_domain.length_x.to("meter").magnitude,
            packing_domain.length_y.to("meter").magnitude,
            packing_domain.length_z.to("meter").magnitude,
        ],
        dtype=float,
    )
    periodic = bool(packing_domain.use_periodic_boundaries)
    _validate_arrays(positions, radii, classes, lengths, periodic)
    archive_metadata = {
        "schema": _SCHEMA,
        "version": _VERSION,
        "source": source,
        "run_metadata": run_metadata,
        "user_metadata": metadata or {},
        "units": _SI_UNITS,
    }
    try:
        metadata_json = json.dumps(archive_metadata, sort_keys=True)
    except (TypeError, ValueError) as error:
        raise TypeError("metadata and run metadata must be JSON-serializable.") from error
    np.savez_compressed(
        target,
        positions_m=positions,
        radii_m=radii,
        classes_index=classes,
        box_lengths_m=lengths,
        periodic=np.asarray(periodic),
        metadata_json=np.asarray(metadata_json),
    )


def load_packing(path: str | Path) -> LoadedPackingConfiguration:
    """
    Load and validate a compressed PackLab packing archive.

    Parameters
    ----------
    path : str or pathlib.Path
        Existing ``.npz`` archive produced by :func:`save_packing`.

    Returns
    -------
    LoadedPackingConfiguration
        Validated native configuration and domain. No artificial run
        statistics are reconstructed.

    Raises
    ------
    ValueError
        If schema metadata, array shapes, physical ranges, domain containment,
        or hard-sphere non-overlap constraints are invalid.
    """
    target = Path(path)
    if target.suffix != ".npz":
        raise ValueError("path must end in '.npz'.")
    try:
        with np.load(target, allow_pickle=False) as archive:
            if set(archive.files) != _REQUIRED_ARCHIVE_KEYS:
                missing = _REQUIRED_ARCHIVE_KEYS - set(archive.files)
                raise ValueError(f"archive keys are invalid; missing {sorted(missing)}.")
            positions = np.asarray(archive["positions_m"], dtype=float)
            radii = np.asarray(archive["radii_m"], dtype=float)
            classes = np.asarray(archive["classes_index"])
            lengths = np.asarray(archive["box_lengths_m"], dtype=float)
            periodic_value = np.asarray(archive["periodic"])
            if periodic_value.shape != () or not np.issubdtype(periodic_value.dtype, np.bool_):
                raise ValueError("archive periodic flag must be a scalar boolean.")
            periodic = bool(periodic_value.item())
            try:
                archive_metadata = json.loads(str(np.asarray(archive["metadata_json"]).item()))
            except (json.JSONDecodeError, ValueError) as error:
                raise ValueError("archive metadata_json must contain valid JSON.") from error
    except OSError as error:
        raise ValueError(f"could not read packing archive '{target}'.") from error

    if not isinstance(archive_metadata, dict) or archive_metadata.get("schema") != _SCHEMA:
        raise ValueError("archive has an unsupported packing schema.")
    if archive_metadata.get("version") != _VERSION:
        raise ValueError("archive has an unsupported packing schema version.")
    if archive_metadata.get("source") not in {"rsa", "metropolis", "unknown"}:
        raise ValueError("archive metadata source is invalid.")
    if archive_metadata.get("units") != _SI_UNITS:
        raise ValueError("archive metadata must declare SI metre units.")
    if not isinstance(archive_metadata.get("run_metadata"), dict):
        raise ValueError("archive run_metadata must be a JSON object.")
    _validate_arrays(positions, radii, classes, lengths, periodic)
    domain = PackingDomain(*(lengths * ureg.meter), use_periodic_boundaries=periodic)
    configuration = PackingConfiguration.from_arrays(
        positions * ureg.meter, radii * ureg.meter, classes
    )
    return LoadedPackingConfiguration(configuration, domain, archive_metadata)
