from typing import Any
from PackLab.monte_carlo.domain import PackingDomain
from PackLab.monte_carlo.simulator import RSAOptions
from PackLab.samplers import RadiusSampler

class PackingEstimate:
    centers: Any
    mean_g: Any
    std_g: Any
    number_of_species: int
    number_of_bins: int

class PackingEstimator:
    def __init__(self, domain: PackingDomain, radius_sampler: RadiusSampler, options: RSAOptions, number_of_bins: int) -> None: ...
    def estimate(self, number_of_samples: int, maximum_pairs: int = 0) -> PackingEstimate: ...
