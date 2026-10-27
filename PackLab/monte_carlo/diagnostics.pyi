from typing import Any, Callable

class MetropolisDiagnosticReport:
    sample_sweeps: Any
    acceptance_rates: Any
    accepted_move_deltas: Any
    rejected_move_deltas: Any
    mean_squared_displacements: Any
    observable: Any | None
    observable_name: str | None
    burn_in_sweeps: int
    autocorrelation_time: float | None
    effective_sample_size: float | None
    block_means: Any | None
    standard_error: float | None
    confidence_interval: tuple[float, float] | None
    confidence_level: float | None
    final_result: Any

def run_metropolis_diagnostics(
    simulator: Any,
    number_of_sweeps: int,
    *,
    sample_interval: int = ...,
    burn_in_sweeps: int = ...,
    observable: Callable[[Any], float] | None = ...,
    observable_name: str = ...,
    block_size: int | None = ...,
    confidence_level: float = ...,
) -> MetropolisDiagnosticReport: ...
