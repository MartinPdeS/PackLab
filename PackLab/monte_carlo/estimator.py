from typing import Optional, Tuple
import numpy as np

from PackLab.monte_carlo.simulator import Simulator
from TypedUnit.units import ureg

class PartialGEstimator:
    """
    Monte Carlo estimator for partial g_ij(r) by repeatedly running RSA and aggregating g_ij.

    This class does not cache RSA simulators or results.
    """

    def __init__(self, *, domain: object, radius_sampler: object, options: object, n_bins: int) -> None:
        self.domain = domain
        self.radius_sampler = radius_sampler
        self.options = options
        self.n_bins = int(n_bins)

    def estimate(self, n_sample: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Run n_sample independent simulations and compute mean and standard deviation of g_ij(r).

        Returns
        - centers: (B,)
        - mean_g: (K, K, B)
        - std_g: (K, K, B)  (sample standard deviation, ddof=1 when n_sample > 1)
        """
        if n_sample <= 0:
            raise ValueError("n_sample must be > 0")

        centers: Optional[np.ndarray] = None

        # Welford accumulators
        sample_count: int = 0
        mean_g: Optional[np.ndarray] = None
        m2_g: Optional[np.ndarray] = None

        rsa_simulator = Simulator(
            domain=self.domain,
            radius_sampler=self.radius_sampler,
            options=self.options,
        )

        for _ in range(n_sample):
            rsa_simulator._cpp_reset()
            result = rsa_simulator.run()

            centers_i, g_matrix_i = result.binding.compute_partial_pair_correlation_function(n_bins=self.n_bins, maximum_pairs=0)

            sample_count += 1

            if mean_g is None:
                centers = np.asarray(centers_i)
                mean_g = g_matrix_i.copy()
                m2_g = np.zeros_like(mean_g, dtype=float)
                continue

            # Welford update
            delta = g_matrix_i - mean_g
            mean_g += delta / float(sample_count)
            delta2 = g_matrix_i - mean_g
            m2_g += delta * delta2

        if n_sample > 1:
            variance_g = m2_g / float(n_sample - 1)  # sample variance
        else:
            variance_g = np.zeros_like(mean_g, dtype=float)

        std_g = np.sqrt(variance_g)
        return centers * ureg.meter, mean_g, std_g