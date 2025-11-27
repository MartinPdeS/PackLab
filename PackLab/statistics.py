from __future__ import annotations

from typing import Any, Dict, Optional, Iterable, Tuple
from dataclasses import dataclass

from tabulate import tabulate


@dataclass(frozen=True)
class Statistics:
    sphere_count: int
    packing_fraction_geometry: float
    packing_fraction_simulator: Optional[float]
    attempted_insertions: Optional[int]
    accepted_insertions: Optional[int]
    rejected_insertions: Optional[int]
    consecutive_rejections: Optional[int]
    radius_min: float
    radius_max: float
    radius_mean: float
    radius_median: float
    radius_std: float

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sphere_count": self.sphere_count,
            "packing_fraction_geometry": self.packing_fraction_geometry,
            "packing_fraction_simulator": self.packing_fraction_simulator,
            "attempted_insertions": self.attempted_insertions,
            "accepted_insertions": self.accepted_insertions,
            "rejected_insertions": self.rejected_insertions,
            "consecutive_rejections": self.consecutive_rejections,
            "radius_min": self.radius_min,
            "radius_max": self.radius_max,
            "radius_mean": self.radius_mean,
            "radius_median": self.radius_median,
            "radius_std": self.radius_std,
        }

    def format_table(
        self,
        float_format: str = ".6g",
        table_format: str = "github",
        sort_keys: bool = False,
        include_none_rows: bool = True,
    ) -> str:
        def format_value(value: Any) -> str:
            if value is None:
                return "None"
            if isinstance(value, float):
                return format(value, float_format)
            return str(value)

        items: Iterable[Tuple[str, Any]] = self.to_dict().items()
        if sort_keys:
            items = sorted(items, key=lambda pair: pair[0])

        rows = []
        for key, value in items:
            if (value is None) and (not include_none_rows):
                continue
            rows.append((key, format_value(value)))

        return tabulate(rows, headers=("Metric", "Value"), tablefmt=table_format)

    def print(
        self,
        float_format: str = ".6g",
        table_format: str = "github",
        sort_keys: bool = False,
        include_none_rows: bool = True,
    ) -> None:
        print(self.format_table(
            float_format=float_format,
            table_format=table_format,
            sort_keys=sort_keys,
            include_none_rows=include_none_rows,
        ))
