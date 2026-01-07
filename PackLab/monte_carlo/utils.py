#!/usr/bin/env python
# -*- coding: utf-8 -*-
from typing import Tuple
import numpy as np

from pydantic import ConfigDict

# Configuration dictionary for the Pydantic dataclass
config_dict = ConfigDict(
    kw_only=True, slots=True, extra="forbid", arbitrary_types_allowed=True
)



def _minimum_image_displacement(delta: np.ndarray, box_length: float) -> np.ndarray:
    return delta - box_length * np.round(delta / box_length)

def _pair_distance_pdf(
    positions: np.ndarray,
    box_lengths: Tuple[float, float, float],
    use_periodic_boundaries: bool,
    maximum_pairs: int,
    bins: int,
    random_seed: int,
) -> Tuple[np.ndarray, np.ndarray]:
    positions = np.asarray(positions, dtype=float)
    n = positions.shape[0]
    if n < 2:
        return np.array([]), np.array([])

    random_generator = np.random.default_rng(random_seed)

    total_pairs = n * (n - 1) // 2
    compute_all = total_pairs <= maximum_pairs

    if compute_all:
        i_indices, j_indices = np.triu_indices(n, k=1)
    else:
        i_indices = random_generator.integers(0, n, size=maximum_pairs)
        j_indices = random_generator.integers(0, n, size=maximum_pairs)
        valid = i_indices != j_indices
        i_indices = i_indices[valid]
        j_indices = j_indices[valid]

    delta = positions[i_indices] - positions[j_indices]

    if use_periodic_boundaries:
        delta[:, 0] = _minimum_image_displacement(delta[:, 0], box_lengths[0])
        delta[:, 1] = _minimum_image_displacement(delta[:, 1], box_lengths[1])
        delta[:, 2] = _minimum_image_displacement(delta[:, 2], box_lengths[2])

    distances = np.sqrt(np.sum(delta * delta, axis=1))
    hist, edges = np.histogram(distances, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:] )
    return centers, hist