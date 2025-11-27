from typing import Tuple

from dataclasses import dataclass

@dataclass(frozen=True)
class DomainInfo:
    length_x: float
    length_y: float
    length_z: float
    use_periodic_boundaries: bool

    @property
    def box_lengths(self) -> Tuple[float, float, float]:
        return (self.length_x, self.length_y, self.length_z)

    @property
    def volume(self) -> float:
        return self.length_x * self.length_y * self.length_z
