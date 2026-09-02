from typing import Any, Literal

class PercusYevickResult:
    radii: Any
    densities: Any
    wavenumber: Any
    distances: Any
    g: Any
    h: Any
    H: Any

class PercusYevickSolver:
    def __init__(
        self,
        densities: Any,
        radii: Any,
        wavenumber: Any | Literal["auto"] = "auto",
        radial_resolution: Any | None = None,
        samples_per_oscillation: int = 12,
    ) -> None: ...
    def compute(self, distances: Any) -> PercusYevickResult: ...
