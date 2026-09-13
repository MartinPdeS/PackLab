from pathlib import Path
from typing import Any

from PackLab.monte_carlo.domain import PackingDomain
from PackLab.monte_carlo.simulator import PackingConfiguration

class LoadedPackingConfiguration:
    sphere_configuration: PackingConfiguration
    domain: PackingDomain
    metadata: dict[str, Any]
    source: str
    run_metadata: dict[str, Any]
    positions: Any
    radii: Any
    classes_index: Any
    def compute_partial_pair_correlation_function(self, **kwargs: Any) -> tuple[Any, Any]: ...
    def save(self, path: str | Path, **kwargs: Any) -> None: ...

def save_packing(
    packing: Any,
    path: str | Path,
    *,
    domain: PackingDomain | None = ...,
    source: str | None = ...,
    metadata: dict[str, Any] | None = ...,
) -> None: ...
def load_packing(path: str | Path) -> LoadedPackingConfiguration: ...
