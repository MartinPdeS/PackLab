"""Low-memory diagnostics for fixed-volume Metropolis hard-sphere sampling."""

from __future__ import annotations

from dataclasses import dataclass
from statistics import NormalDist
from typing import Any, Callable

import numpy as np

from PackLab.units import ureg


@dataclass(frozen=True, slots=True)
class MetropolisDiagnosticReport:
    """
    Retained diagnostics from a sampled Metropolis trajectory.

    Notes
    -----
    These summaries quantify the recorded observable only. They do not prove
    or automatically assert equilibration.
    """

    sample_sweeps: np.ndarray
    acceptance_rates: np.ndarray
    accepted_move_deltas: np.ndarray
    rejected_move_deltas: np.ndarray
    mean_squared_displacements: Any
    observable: np.ndarray
    observable_name: str
    burn_in_sweeps: int
    autocorrelation_time: float
    effective_sample_size: float
    block_means: np.ndarray
    standard_error: float
    confidence_interval: tuple[float, float]
    confidence_level: float
    final_result: Any


def _validate_integer(name: str, value: int, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer greater than or equal to {minimum}.")
    return int(value)


def _autocorrelation_time(values: np.ndarray) -> float:
    centered = values - values.mean()
    variance = np.dot(centered, centered) / values.size
    if variance == 0.0:
        return 0.5
    correlations = np.correlate(centered, centered, mode="full")[values.size - 1 :] / (
        variance * np.arange(values.size, 0, -1)
    )
    tau = 0.5
    for correlation in correlations[1:]:
        if correlation <= 0.0:
            break
        tau += float(correlation)
    return tau


def run_metropolis_diagnostics(
    simulator: Any,
    number_of_sweeps: int,
    *,
    sample_interval: int = 10,
    burn_in_sweeps: int = 0,
    observable: Callable[[Any], float] | None = None,
    observable_name: str = "mean_squared_displacement_m2",
    block_size: int | None = None,
    confidence_level: float = 0.95,
) -> MetropolisDiagnosticReport:
    """
    Run sampled chunks of a Metropolis simulation without retaining positions.

    Parameters
    ----------
    simulator : MetropolisSimulator
        Native public simulator to advance. It is not reset by this function.
    number_of_sweeps : int
        Number of additional sweeps to execute.
    sample_interval : int, default=10
        Number of sweeps between retained observations.
    burn_in_sweeps : int, default=0
        Initial sweeps excluded from autocorrelation and block statistics.
    observable : callable, optional
        Function taking the current ``PackingConfiguration`` and returning one
        finite scalar. The default records mean squared displacement (MSD).
    observable_name : str, default="mean_squared_displacement_m2"
        Label for the scalar observable in the report.
    block_size : int, optional
        Number of post-burn-in samples per block. By default it is selected
        from the estimated autocorrelation time.
    confidence_level : float, default=0.95
        Open-interval normal confidence level for the block-mean estimate.

    Returns
    -------
    MetropolisDiagnosticReport
        Move deltas, per-sample rates, MSD relative to the starting state, and
        autocorrelation/block summaries. The report does not establish
        equilibration.
    """
    total = _validate_integer("number_of_sweeps", number_of_sweeps, 1)
    interval = _validate_integer("sample_interval", sample_interval, 1)
    burn_in = _validate_integer("burn_in_sweeps", burn_in_sweeps, 0)
    if burn_in >= total:
        raise ValueError("burn_in_sweeps must be smaller than number_of_sweeps.")
    if not 0.0 < confidence_level < 1.0:
        raise ValueError("confidence_level must lie strictly between 0 and 1.")
    if not callable(getattr(simulator, "run_sweeps", None)):
        raise TypeError("simulator must be a MetropolisSimulator exposing run_sweeps().")

    start = np.asarray(
        simulator.sphere_configuration.positions.to("meter").magnitude, dtype=float
    ).copy()
    lengths = np.asarray(
        [
            simulator.domain.length_x.to("meter").magnitude,
            simulator.domain.length_y.to("meter").magnitude,
            simulator.domain.length_z.to("meter").magnitude,
        ],
        dtype=float,
    )
    periodic = bool(simulator.domain.use_periodic_boundaries)
    previous_accepted = int(simulator.statistics.accepted_moves)
    previous_rejected = int(simulator.statistics.rejected_moves)
    sweeps: list[int] = []
    rates: list[float] = []
    accepted: list[int] = []
    rejected: list[int] = []
    msd: list[float] = []
    observations: list[float] = []
    completed = 0
    final_result = None
    while completed < total:
        chunk = min(interval, total - completed)
        final_result = simulator.run_sweeps(chunk)
        completed += chunk
        current_statistics = simulator.statistics
        accepted_delta = current_statistics.accepted_moves - previous_accepted
        rejected_delta = current_statistics.rejected_moves - previous_rejected
        attempted_delta = accepted_delta + rejected_delta
        positions = np.asarray(
            simulator.sphere_configuration.positions.to("meter").magnitude, dtype=float
        )
        displacement = positions - start
        if periodic:
            displacement -= lengths * np.round(displacement / lengths)
        current_msd = float(
            np.mean(np.einsum("ij,ij->i", displacement, displacement))
        )
        value = (
            current_msd if observable is None else float(observable(simulator.sphere_configuration))
        )
        if not np.isfinite(value):
            raise ValueError("observable must return a finite scalar.")
        sweeps.append(completed)
        rates.append(accepted_delta / attempted_delta if attempted_delta else 0.0)
        accepted.append(accepted_delta)
        rejected.append(rejected_delta)
        msd.append(current_msd)
        observations.append(value)
        previous_accepted = int(current_statistics.accepted_moves)
        previous_rejected = int(current_statistics.rejected_moves)

    sample_sweeps = np.asarray(sweeps, dtype=int)
    selected = sample_sweeps > burn_in
    post_burn_values = np.asarray(observations, dtype=float)[selected]
    if post_burn_values.size < 2:
        raise ValueError("burn-in leaves fewer than two retained samples for diagnostics.")
    tau = _autocorrelation_time(post_burn_values)
    ess = post_burn_values.size / (2.0 * tau)
    selected_block_size = (
        int(np.ceil(2.0 * tau))
        if block_size is None
        else _validate_integer("block_size", block_size, 1)
    )
    blocks = [
        post_burn_values[index : index + selected_block_size].mean()
        for index in range(0, post_burn_values.size, selected_block_size)
        if post_burn_values[index : index + selected_block_size].size == selected_block_size
    ]
    if len(blocks) < 2:
        raise ValueError(
            "fewer than two complete blocks remain; reduce block_size or collect more samples."
        )
    block_means = np.asarray(blocks)
    standard_error = float(block_means.std(ddof=1) / np.sqrt(block_means.size))
    z_value = NormalDist().inv_cdf((1.0 + confidence_level) / 2.0)
    mean = float(block_means.mean())
    return MetropolisDiagnosticReport(
        sample_sweeps=sample_sweeps,
        acceptance_rates=np.asarray(rates),
        accepted_move_deltas=np.asarray(accepted, dtype=int),
        rejected_move_deltas=np.asarray(rejected, dtype=int),
        mean_squared_displacements=np.asarray(msd) * ureg.meter**2,
        observable=np.asarray(observations),
        observable_name=observable_name,
        burn_in_sweeps=burn_in,
        autocorrelation_time=tau,
        effective_sample_size=ess,
        block_means=block_means,
        standard_error=standard_error,
        confidence_interval=(mean - z_value * standard_error, mean + z_value * standard_error),
        confidence_level=confidence_level,
        final_result=final_result,
    )
