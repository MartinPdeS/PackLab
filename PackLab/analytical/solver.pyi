from typing import Any

class PercusYevickResult:
    radii: Any
    densities: Any
    p: Any
    distances: Any
    g: Any
    h: Any
    H: Any

class PercusYevickSolver:
    def __init__(self, densities: Any, radii: Any, p: Any) -> None: ...
    def compute(self, distances: Any) -> PercusYevickResult: ...
