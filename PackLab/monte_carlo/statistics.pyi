class PackingStatistics:
    attempted_insertions: int
    accepted_insertions: int
    rejected_insertions: int
    consecutive_rejections: int
    sphere_count: int
    packing_fraction_geometry: float
    packing_fraction_simulator: float
    radius_min: float
    radius_max: float
    radius_mean: float
    radius_median: float
    radius_std: float
    total_runtime_seconds: float
    def print(self) -> None: ...
