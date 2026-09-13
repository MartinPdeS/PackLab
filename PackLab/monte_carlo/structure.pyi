from typing import Any, Literal

class EmpiricalStructure:
    source: str
    densities: Any
    radii: Any
    classes: Any
    distances: Any
    g: Any
    h: Any
    wavenumber: Any
    H: Any
    S: Any
    def scattering_inputs(self) -> tuple[Any, Any, Any]: ...

def empirical_structure(
    packing: Any,
    *,
    domain: Any | None = ...,
    n_bins: int = ...,
    maximum_pairs: int = ...,
    wavenumber: Any | Literal["auto"] = ...,
    samples_per_oscillation: int = ...,
) -> EmpiricalStructure: ...
