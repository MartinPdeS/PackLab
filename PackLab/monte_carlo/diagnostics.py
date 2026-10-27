"""Low-memory diagnostics for fixed-volume Metropolis hard-sphere sampling."""


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
    observable: np.ndarray | None
    observable_name: str | None
    burn_in_sweeps: int
    autocorrelation_time: float | None
    effective_sample_size: float | None
    block_means: np.ndarray | None
    standard_error: float | None
    confidence_interval: tuple[float, float] | None
    confidence_level: float | None
    final_result: Any


def _validate_integer(name: str, value: int, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer greater than or equal to {minimum}.")
    return int(value)


def _autocorrelation_time(values: np.ndarray) -> float | None:
    centered = values - values.mean()
    variance = np.dot(centered, centered) / values.size
    if variance == 0.0:
        return None
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
    observable_name: str = "observable",
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
        finite scalar. Supply an approximately stationary observable to obtain
        autocorrelation and uncertainty estimates. If omitted, only mobility
        and move diagnostics are reported.
    observable_name : str, default="observable"
        Label for the supplied scalar observable.
    block_size : int, optional
        Number of post-burn-in samples per block. By default it is selected
        from the estimated autocorrelation time.
    confidence_level : float, default=0.95
        Open-interval normal confidence level for the block-mean estimate.

    Returns
    -------
    MetropolisDiagnosticReport
        Move deltas, per-sample rates, MSD relative to the starting state, and
        autocorrelation/block summaries when an observable is supplied; those
        fields are otherwise ``None``. MSD uses accumulated accepted moves,
        including periodic crossings, relative to the start of this call.
        A final sample at a shorter interval is retained but excluded from
        autocorrelation and block statistics. The report does not establish
        stationarity or equilibration; uncertainty estimates rely on both.
        If the retained observable is constant, its autocorrelation and
        uncertainty are indeterminate and the corresponding fields are ``None``.
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
    if observable is not None and not callable(observable):
        raise TypeError("observable must be callable.")
    if block_size is not None:
        block_size = _validate_integer("block_size", block_size, 1)
        if observable is None:
            raise ValueError("block_size requires an observable.")
    if observable is not None:
        retained_count = total // interval - burn_in // interval
        if retained_count < 2:
            raise ValueError("burn-in leaves fewer than two regularly spaced samples.")
        if block_size is not None and retained_count // block_size < 2:
            raise ValueError(
                "fewer than two complete blocks remain; reduce block_size or collect more samples."
            )

    start_displacement = np.asarray(
        simulator.cumulative_displacements.to("meter").magnitude, dtype=float
    ).copy()
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
        displacement = np.asarray(
            simulator.cumulative_displacements.to("meter").magnitude, dtype=float
        ) - start_displacement
        current_msd = float(
            np.mean(np.einsum("ij,ij->i", displacement, displacement))
        )
        if observable is not None:
            value = float(observable(simulator.sphere_configuration))
            if not np.isfinite(value):
                raise ValueError("observable must return a finite scalar.")
            observations.append(value)
        sweeps.append(completed)
        rates.append(accepted_delta / attempted_delta if attempted_delta else 0.0)
        accepted.append(accepted_delta)
        rejected.append(rejected_delta)
        msd.append(current_msd)
        previous_accepted = int(current_statistics.accepted_moves)
        previous_rejected = int(current_statistics.rejected_moves)

    sample_sweeps = np.asarray(sweeps, dtype=int)
    if observable is None:
        return MetropolisDiagnosticReport(
            sample_sweeps=sample_sweeps,
            acceptance_rates=np.asarray(rates),
            accepted_move_deltas=np.asarray(accepted, dtype=int),
            rejected_move_deltas=np.asarray(rejected, dtype=int),
            mean_squared_displacements=np.asarray(msd) * ureg.meter**2,
            observable=None,
            observable_name=None,
            burn_in_sweeps=burn_in,
            autocorrelation_time=None,
            effective_sample_size=None,
            block_means=None,
            standard_error=None,
            confidence_interval=None,
            confidence_level=None,
            final_result=final_result,
        )
    selected = (sample_sweeps > burn_in) & (sample_sweeps % interval == 0)
    post_burn_values = np.asarray(observations, dtype=float)[selected]
    if post_burn_values.size < 2:
        raise ValueError("burn-in leaves fewer than two retained samples for diagnostics.")
    tau = _autocorrelation_time(post_burn_values)
    ess = post_burn_values.size / (2.0 * tau) if tau is not None else None
    selected_block_size = (
        (int(np.ceil(2.0 * tau)) if tau is not None else 1)
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
    standard_error = None
    confidence_interval = None
    if tau is not None:
        standard_error = float(block_means.std(ddof=1) / np.sqrt(block_means.size))
        z_value = NormalDist().inv_cdf((1.0 + confidence_level) / 2.0)
        mean = float(block_means.mean())
        confidence_interval = (mean - z_value * standard_error, mean + z_value * standard_error)
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
        confidence_interval=confidence_interval,
        confidence_level=confidence_level if confidence_interval is not None else None,
        final_result=final_result,
    )
