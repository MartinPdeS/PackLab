from typing import Any

from PackLab.monte_carlo.domain import PackingDomain
from PackLab.monte_carlo.simulator import PackingConfiguration


class MetropolisOptions:
    random_seed: int
    number_of_sweeps: int
    maximum_displacement: Any


class MetropolisStatistics:
    attempted_moves: int
    accepted_moves: int
    rejected_moves: int
    completed_sweeps: int
    total_runtime_seconds: float
    acceptance_rate: float


class MetropolisSimulator:
    def __init__(
        self,
        domain: PackingDomain,
        initial_configuration: PackingConfiguration,
        options: MetropolisOptions,
    ) -> None: ...
    sphere_configuration: PackingConfiguration
    statistics: MetropolisStatistics
    def reset(self) -> None: ...
    def run(self) -> Any: ...
