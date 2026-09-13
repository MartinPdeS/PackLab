from pathlib import Path
from typing import Any

class PackingResult:
    binding: Any
    source: str
    run_metadata: dict[str, Any]
    positions: Any
    radii: Any
    statistics: Any
    sphere_configuration: Any
    domain: Any
    partial_volume_fractions: Any
    partial_volumes: Any
    def __init__(
        self, binding: Any, source: str | None = ..., run_metadata: dict[str, Any] | None = ...
    ) -> None: ...
    def save(self, path: str | Path, **kwargs: Any) -> None: ...
    def compute_partial_pair_correlation_function(self, **kwargs: Any) -> tuple[Any, Any]: ...
