from typing import Any

class PackingResult:
    radii: Any
    statistics: Any
    sphere_configuration: Any
    partial_volume_fractions: Any
    partial_volumes: Any
    classes: Any
    number_of_classes: int
    domain: Any
    def compute_partial_pair_distances(self, maximum_pairs: int) -> tuple[Any, Any, Any]: ...
    def compute_partial_pair_correlation_function(self, n_bins: int, maximum_pairs: int = 1_000_000) -> tuple[Any, Any]: ...
