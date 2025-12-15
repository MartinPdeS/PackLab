from dataclasses import dataclass
from typing import Any, List, Tuple
from TypedUnit import ureg
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, sin, cos

from MPSPlots import helper



@dataclass(slots=True)
class PercusYevickResult:
    """Container for Percus Yevick mixture outputs.

    This object is responsible for plotting, since it owns the computed arrays.

    Attributes:
        epsilons: List [e_0, e_1, e_2, e_3].
        R_ij: Pairwise sum matrix from radii.
        S_ij: Pairwise difference matrix from radii.
        A_i: Species coefficients A_i.
        B_i: Species coefficients B_i.
        D_ij: Pairwise coefficients D_ij.
        Cpy: C(p) tensor, shape (N, N, P).
        H: H(p) tensor, shape (N, N, P).
        h: h(r) tensor, shape (N, N, Nr).
        g: g(r) tensor, shape (N, N, Nr).
        distances: r grid.
        p: p grid.
        densities: densities used.
        radii: radii used.
    """

    epsilons: Any
    R_ij: Any
    S_ij: Any
    A_i: Any
    B_i: Any
    D_ij: Any
    Cpy: Any
    H: Any
    h: Any
    g: Any
    distances: Any
    p: Any
    densities: Any
    radii: Any

    @helper.post_mpl_plot
    def plot_pair_correlation(self) -> None:
        """Plot all g[i, j](r) curves on a single axis.

        Args:
            show_legend: Whether to display the legend.
            legend_ncol: Number of columns in the legend.
        """
        r = self.distances.magnitude
        n = self.g.shape[0]

        figure, ax = plt.subplots(1, 1)

        for i in range(n):
            for j in range(n):
                ax.plot(r, self.g[i, j, :].magnitude, label=f"g[{i},{j}]")

        ax.set_xlabel("r")
        ax.set_ylabel("g(r)")
        ax.set_title("Radial distribution functions")

        ax.legend()

        return figure


class PercusYevickSolver:
    """Percus Yevick mixture helper for computing C(p), H(p), h(r), and g(r).

    Notes:
        The radii provided at construction time are always used.
        Debug printing is only produced if you call debug_Cpy_units.

    Attributes:
        densities: Species densities.
        radii: Species radii.
        p: Fourier grid.
    """

    def __init__(self, densities: Any, radii: Any, p: Any):
        self.densities = densities
        self.radii = radii
        self.p = p

    def get_epsilons(self) -> List[Any]:
        """Compute epsilons e_0..e_3 using self.densities and self.radii."""
        epsilons = []
        for alpha in [0, 1, 2, 3]:
            term0 = pi / 6
            term1 = self.densities * (2 * self.radii) ** alpha
            epsilon_i = term0 * np.sum(term1)
            epsilons.append(epsilon_i)
        return epsilons

    def get_parameters(self, epsilons: List[Any]) -> Tuple[Any, Any, Any, Any, Any]:
        """Compute (R_ij, S_ij, A_i, B_i, D_ij) from epsilons and self.radii."""
        radii = self.radii
        radii_ = np.broadcast_to(radii.magnitude, (radii.size, radii.size)) * radii.units

        R_ij = radii_ + radii_.T
        S_ij = radii_.T - radii_

        e_0, e_1, e_2, e_3 = epsilons
        denominator = (1 - e_3)

        A_i = (1 - e_3 + 6 * radii * e_2) / denominator ** 2
        B_i = -6 * radii ** 2 * e_2 / denominator ** 2

        D_ij = -(A_i[:, np.newaxis] * R_ij ** 2) / 2 - (B_i[:, np.newaxis] * R_ij)

        return R_ij, S_ij, A_i, B_i, D_ij

    @staticmethod
    def _safe_sin_over_x(x, small: float = 1e-12) -> np.ndarray:
        """Compute sin(x)/x for Pint Quantity x, using the limit 1 at x -> 0."""
        x_ = x.magnitude
        out = np.empty_like(x_, dtype=float)
        mask = np.abs(x_) > small
        out[mask] = np.sin(x_[mask]) / x_[mask]
        out[~mask] = 1.0
        return out / x.units

    @staticmethod
    def _safe_sin_over_x3_minus_cos_over_x2(x, small: float = 1e-12) -> np.ndarray:
        """Compute sin(x)/x^3 - cos(x)/x^2 for Pint Quantity x, limit 1/3 at x -> 0."""
        x_ = x.magnitude
        out = np.empty_like(x_, dtype=float)
        mask = np.abs(x_) > small
        out[mask] = np.sin(x_[mask]) / (x_[mask] ** 3) - np.cos(x_[mask]) / (x_[mask] ** 2)
        out[~mask] = 1.0 / 3.0
        return out / x.units ** 2

    def get_Cpy_(
        self,
        index_i: int,
        index_j: int,
        epsilons: List[Any],
        debug_units: bool = False,
    ) -> Any:
        """Compute one curve Cpy[i, j, :] using your implementation."""
        densities = self.densities
        radii = self.radii
        p = self.p

        shape = (radii.size, p.size)

        R = 2 * np.broadcast_to(radii[:, np.newaxis], shape=shape)
        p_grid = np.broadcast_to(p[np.newaxis, :], shape=shape)

        X = R * p_grid / 2
        X = X.to("dimensionless")

        sin_over_X = self._safe_sin_over_x(X)
        sin_over_X3_minus_cos_over_X2 = self._safe_sin_over_x3_minus_cos_over_x2(X)

        N = R ** 2 * sin_over_X
        M = 3 * R ** 3 * sin_over_X3_minus_cos_over_X2

        X_i = X[index_i]
        X_j = X[index_j]

        n_i = densities[index_i]
        n_j = densities[index_j]

        R_i = R[index_i]
        R_j = R[index_j]

        e_0, e_1, e_2, e_3 = epsilons
        denominator = 1 - e_3

        N_i = N[index_i]
        N_j = N[index_j]

        M_i = M[index_i]
        M_j = M[index_j]

        if debug_units:
            def _u(name, value):
                units = getattr(value, "units", "no units")
                print(f"{name}: {units}")

            print("---- get_Cpy_ units ----")
            _u("radii", radii)
            _u("p", p)
            _u("densities", densities)
            _u("R", R)
            _u("p_grid", p_grid)
            _u("X", X)
            print("safe kernels: sin_over_X, sin_over_X3_minus_cos_over_X2 are dimensionless float arrays")
            _u("N", N)
            _u("M", M)
            _u("X_i", X_i)
            _u("X_j", X_j)
            _u("n_i", n_i)
            _u("n_j", n_j)
            _u("R_i", R_i)
            _u("R_j", R_j)
            print("epsilons:")
            for k, e in enumerate(epsilons):
                _u(f"e_{k}", e)
            print(f"denominator: {getattr(denominator, 'units', 'no units')}")
            _u("N_i", N_i)
            _u("N_j", N_j)
            _u("M_i", M_i)
            _u("M_j", M_j)
            print("------------------------")

        term0 = -(pi / 6) * np.sqrt(n_i * n_j) / denominator

        term2 = M_j * (
            cos(X_i)
            + X_i * sin(X_i)
            + 3 * e_2 * R_i * cos(X_i) / denominator
            + 3 * e_1 * N_i / denominator
            + 9 * e_2 ** 2 * N_i / denominator ** 2
        )

        term3 = M_i * (
            cos(X_j)
            + X_j * sin(X_j)
            + 3 * e_2 * R_j * cos(X_j) / denominator
            + 3 * e_1 * N_j / denominator
            + 9 * e_2 ** 2 * N_j / denominator ** 2
        )

        term4 = M_i * M_j * (
            e_0 / denominator
            + p ** 2 * e_2 / (4 * denominator)
            + 6 * e_1 * e_2 / denominator ** 2
            + 9 * e_2 ** 3 / denominator ** 3
        )

        term5 = 3 * N_i * R_j * cos(X_j)
        term6 = 3 * N_j * R_i * cos(X_i)

        term7 = 9 * e_2 * N_i * N_j / denominator

        output = term0 * (term2 + term3 + term4 + term5 + term6 + term7)

        return output.to("dimensionless")

    def get_Cpy(self, epsilons: List[Any]) -> Any:
        """Compute the full Cpy tensor with shape (N, N, P)."""
        radii = self.radii
        p = self.p

        Cpy = np.zeros([radii.size, radii.size, p.size])

        for index_i in range(radii.size):
            for index_j in range(radii.size):
                Cpy[index_i, index_j, :] = self.get_Cpy_(
                    index_i=index_i,
                    index_j=index_j,
                    epsilons=epsilons,
                    debug_units=False,
                ).magnitude

        return Cpy * ureg.dimensionless

    @staticmethod
    def solve_H_from_C_batch(C: Any) -> Any:
        """Solve H = C + C H for each p slice independently."""
        I = np.eye(C.shape[0])[:, :, None]
        A = I - C

        output = np.zeros(C.shape)

        for i in range(C.shape[-1]):
            res = np.linalg.solve(A[:, :, i], C[:, :, i])
            output[:, :, i] = res.to("dimensionless").magnitude

        return output

    @staticmethod
    def get_radial_fourier_of_H(H: Any, distances: Any, p: Any, densities: Any) -> Any:
        """Compute h(r) from H(p) using your radial transform implementation."""
        if H.shape[-1] != p.size:
            raise ValueError(f"H last axis ({H.shape[-1]}) must match p.size ({p.size}).")

        r_p = (distances[:, None] * p[None, :]).to("").magnitude

        kernel = np.ones_like(r_p, dtype=float)
        mask = r_p != 0.0
        kernel[mask] = np.sin(r_p[mask]) / r_p[mask]

        integrand = (
            H[:, :, None, :] *
            (p ** 2)[None, None, None, :] *
            kernel[None, None, :, :]
        )

        integral = np.trapezoid(y=integrand, x=p, axis=-1)

        output = integral
        output = output / (4 * pi)

        densities_factor = np.sqrt(np.outer(densities.magnitude, densities.magnitude)) * densities.units
        densities_factor = np.broadcast_to(densities_factor[:, :, np.newaxis], output.shape)

        return (output / densities_factor).to("dimensionless")

    def compute(self, distances: Any) -> PercusYevickResult:
        """Compute epsilons, parameters, Cpy, H, h, g and return a result object."""
        epsilons = self.get_epsilons()
        R_ij, S_ij, A_i, B_i, D_ij = self.get_parameters(epsilons=epsilons)

        Cpy = self.get_Cpy(epsilons=epsilons)
        H = self.solve_H_from_C_batch(Cpy)

        h = self.get_radial_fourier_of_H(H=H, distances=distances, p=self.p, densities=self.densities)
        g = 1 + h

        return PercusYevickResult(
            epsilons=epsilons,
            R_ij=R_ij,
            S_ij=S_ij,
            A_i=A_i,
            B_i=B_i,
            D_ij=D_ij,
            Cpy=Cpy,
            H=H,
            h=h,
            g=g,
            distances=distances,
            p=self.p,
            densities=self.densities,
            radii=self.radii,
        )

    def debug_Cpy_units(self, index_i: int = 0, index_j: int = 0) -> None:
        """Print intermediate Pint units for a selected Cpy[i, j, :] computation."""
        epsilons = self.get_epsilons()
        _ = self.get_Cpy_(
            index_i=index_i,
            index_j=index_j,
            epsilons=epsilons,
            debug_units=True,
        )
