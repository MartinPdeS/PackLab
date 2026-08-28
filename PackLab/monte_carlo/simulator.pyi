from typing import Any
from PackLab.monte_carlo.domain import PackingDomain
from PackLab.samplers import RadiusSampler

class PackingConfiguration:
    @property
    def count(self) -> int: ...
    @property
    def positions(self) -> Any: ...
    @property
    def radii(self) -> Any: ...
    @property
    def classes_index(self) -> Any: ...
    @property
    def number_of_classes(self) -> int: ...
    def total_sphere_volume(self) -> Any: ...

class RSAOptions:
    random_seed: int
    maximum_attempts: int
    maximum_spheres: int
    maximum_consecutive_rejections: int
    target_packing_fraction: float
    minimum_center_separation_addition: float
    containment_padding: float
    spatial_grid_cell_size: float
    enforce_radii_distribution: bool

class RSASimulator:
    def __init__(self, domain: PackingDomain, radius_sampler: RadiusSampler, options: RSAOptions) -> None: ...
    sphere_configuration: PackingConfiguration
    def reset(self) -> None: ...
    def run(self) -> Any: ...
