from dataclasses import dataclass
from typing import Any

from TypedUnit import ureg
from numpy.typing import NDArray
import numpy as np


@dataclass(slots=True)
class ScatteringData:
    """Far-field scattering data for a single particle diameter.

    The quantity-bearing fields intentionally use ``Any`` because TypedUnit/Pint
    quantities are generic at runtime and may wrap scalar or array values.
    """

    S1: Any
    S2: Any
    k: Any
    Csca: Any
    phi: Any

    def plot(self, *, tight_layout: bool = True):
        """Plot the magnitudes of the two scattering amplitudes."""
        import matplotlib.pyplot as plt

        phi = self.phi.to("radian").magnitude if hasattr(self.phi, "to") else self.phi
        s1 = np.abs(self.S1.magnitude if hasattr(self.S1, "magnitude") else self.S1)
        s2 = np.abs(self.S2.magnitude if hasattr(self.S2, "magnitude") else self.S2)

        figure, axes = plt.subplots()
        axes.plot(phi, s1, label="|S1|")
        axes.plot(phi, s2, label="|S2|")
        axes.set_xlabel("polar angle (rad)")
        axes.set_ylabel("amplitude")
        axes.legend()
        if tight_layout:
            figure.tight_layout()
        return figure


# Convention note:
# In typical scattering convention, the polar scattering angle is in [0, pi]
# and the azimuthal rotation around the incident axis is in [0, 2*pi).
# In this class, `self.phi` is used with sin(self.phi) weights, so it behaves like a polar angle.
class ScatteringDataset(list):
    """
    Container for multi size scattering data and mixture level post processing.

    This class stores :class:`ScatteringData` instances, one for each diameter. It also provides mixture formulas
    for a number density distribution over diameters and an inter particle correlation
    term H(wavenumber).

    Attributes
    ----------
    e_theta, e_phi : NDArray, shape (2,)
        Unit basis vectors for the polarization basis used in `get_F_matrix`.
        This code treats the scattering amplitude as a 2 component vector in a
        (theta, phi) polarization basis.

    k
        Optical wavenumber in the surrounding medium (1/length).
    phi
        Polar scattering angle grid in radians. Expected range is [0, pi].
    S1, S2
        Arrays assembled by `process()` from each per diameter element.
    Csca
        Per diameter scattering cross section assembled by `process()`.

    Notes
    -----
    Derived arrays are refreshed automatically by the mixture methods. Calling
    :meth:`process` explicitly remains useful when inspecting those arrays.
    """

    e_theta = np.array([1.0, 0.0])
    e_phi   = np.array([0.0, 1.0])

    def process(self):
        """
        Stack per diameter fields into dense arrays.

        After calling this method, the object provides:

        * `self.S1`: ndarray, shape (N, Pphi)
        * `self.S2`: ndarray, shape (N, Pphi)
        * `self.Csca`: Pint Quantity, shape (N,), units meter**2

        where N is the number of diameters (len(self)) and Pphi is the number of
        angular samples in each per diameter result.

        Returns
        -------
        None
        """
        if not self:
            raise ValueError("ScatteringDataset must contain at least one ScatteringData item.")

        if not hasattr(self, "phi"):
            raise ValueError("ScatteringDataset.phi must be set before processing.")

        phi = self._validated_phi()
        expected_shape = phi.shape
        for index, data in enumerate(self):
            if not isinstance(data, ScatteringData):
                raise TypeError("ScatteringDataset items must be ScatteringData instances.")
            if np.asarray(data.S1.magnitude).shape != expected_shape:
                raise ValueError(f"S1 for species {index} must match the phi grid shape {expected_shape}.")
            if np.asarray(data.S2.magnitude).shape != expected_shape:
                raise ValueError(f"S2 for species {index} must match the phi grid shape {expected_shape}.")

        self.S1 = np.asarray([d.S1.magnitude for d in self])
        self.S2 = np.asarray([d.S2.magnitude for d in self])

        self.Csca = np.asarray([d.Csca.to("meter**2").magnitude for d in self]) * ureg.meter ** 2

    def _validated_phi(self) -> np.ndarray:
        """Return the polar-angle grid in radians after validating it."""
        try:
            phi = np.asarray(self.phi.to("radian").magnitude, dtype=float)
        except AttributeError as error:
            raise TypeError("phi must be a unit-bearing angle quantity.") from error
        if phi.ndim != 1 or phi.size < 2:
            raise ValueError("phi must be a one-dimensional angle grid with at least two points.")
        if not np.all(np.isfinite(phi)) or np.any(np.diff(phi) <= 0):
            raise ValueError("phi must be finite and strictly increasing.")
        if phi[0] < 0.0 or phi[-1] > np.pi:
            raise ValueError("phi must lie within the polar-angle interval [0, pi].")
        return phi

    def _validate_mixture_inputs(
        self, densities: Any, H: Any, wavenumber: Any
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Validate and normalize inputs shared by mixture calculations."""
        self.process()
        if not hasattr(self, "k"):
            raise ValueError("ScatteringDataset.k must be set before mixture calculations.")
        try:
            k = self.k.to("1 / meter")
            density_values = np.asarray(densities.to("1 / meter**3").magnitude, dtype=float)
            wavenumber_values = np.asarray(wavenumber.to("1 / meter").magnitude, dtype=float)
        except AttributeError as error:
            raise TypeError("densities and wavenumber must be unit-bearing quantities.") from error
        if not np.isfinite(k.magnitude) or k.magnitude <= 0:
            raise ValueError("k must be finite and positive.")
        if density_values.ndim != 1 or density_values.size != len(self):
            raise ValueError("densities must be one-dimensional with one value per species.")
        if not np.all(np.isfinite(density_values)) or np.any(density_values < 0):
            raise ValueError("densities must be finite and non-negative.")
        if wavenumber_values.ndim != 1 or wavenumber_values.size < 2:
            raise ValueError("wavenumber must be a one-dimensional grid with at least two points.")
        if not np.all(np.isfinite(wavenumber_values)) or np.any(np.diff(wavenumber_values) <= 0):
            raise ValueError("wavenumber must be finite and strictly increasing.")
        H_values = np.asarray(H)
        expected_shape = (len(self), len(self), wavenumber_values.size)
        if H_values.shape != expected_shape:
            raise ValueError(f"H must have shape {expected_shape}; got {H_values.shape}.")
        if not np.all(np.isfinite(H_values)):
            raise ValueError("H must contain only finite values.")
        return density_values, H_values, wavenumber_values

    def get_alpha_beta_factor(self, densities: NDArray) -> tuple[NDArray, NDArray]:
        """
        Build number density factors for mixture sums.

        Parameters
        ----------
        densities : NDArray, shape (N,)
            Species number densities n_alpha. Units should be 1/length**3.
            N must match the number of size classes stored in this object.

        Returns
        -------
        n_alpha : NDArray, shape (N,)
            Same array as input (returned for convenience).
        sqrt_alpha_beta : NDArray, shape (N, N)
            Matrix with entries sqrt(n_alpha * n_beta). This is commonly used in
            symmetric mixture formulations to weight cross correlations.

        Notes
        -----
        This uses the outer product and then takes a square root:

        sqrt_alpha_beta[a, b] = sqrt(n_alpha[a] * n_alpha[b])
        """
        n_alpha = densities
        n_alpha_beta = np.einsum("i,j -> ij", n_alpha, n_alpha)
        sqrt_alpha_beta = np.sqrt(n_alpha_beta)

        return n_alpha, sqrt_alpha_beta

    def get_F_matrix(self, theta_points: int = 100) -> tuple:
        """
        Construct the vector scattering amplitude tensor F for each size, angle, and azimuth.

        This builds a 2 component vector amplitude in the (e_theta, e_phi) basis using
        the convention in your code:

        F(θ_azimuth, φ_polar) = (i / k) * [ e_theta * S2(φ) * cos(θ) - e_phi * S1(φ) * sin(θ) ]

        Parameters
        ----------
        theta_points
            Number of azimuthal sampling points in [0, 2*pi]. This is a rotation
            around the incident axis.

        Returns
        -------
        F : NDArray, complex, shape (2, N, Pphi, Ptheta)
            Vector amplitude tensor.

            Index meaning:
            * first axis (size 2): polarization component (theta, phi)
            * second axis: size class index (N)
            * third axis: polar scattering angle samples (Pphi), from `self.phi`
            * fourth axis: azimuthal angle samples (Ptheta), generated internally

        theta : NDArray, shape (Ptheta,)
            Azimuthal angle grid in radians spanning [0, 2*pi].

        Notes
        -----
        The name `theta` here refers to an azimuthal angle. The polar scattering angle
        is stored in `self.phi`.
        """
        if theta_points < 2:
            raise ValueError("theta_points must be at least 2.")
        self.process()
        theta = np.linspace(0, 2 * np.pi, theta_points)

        _cos = np.cos(theta)
        _sin = np.sin(theta)

        prefactor = 1j / self.k

        term_0 = np.einsum("i, jk, l -> ijkl", self.e_theta, self.S2, _cos)
        term_1 = np.einsum("i, jk, l -> ijkl", self.e_phi,   self.S1, _sin)

        F = prefactor * (term_0 - term_1)

        return F, theta

    def get_mu_independant(self, densities: NDArray):
        """
        Compute the independent scattering contribution to the attenuation coefficient.

        Parameters
        ----------
        densities : NDArray, shape (N,)
            Species number densities n_alpha. Units should be 1/length**3.

        Returns
        -------
        mu_independant_scattering
            Scalar attenuation coefficient contribution with units 1/length.

        Notes
        -----
        This implements a standard independent scattering approximation term:

        mu_s,ind = sum_alpha n_alpha * Csca_alpha

        The supplied density vector is validated against the number of species.
        """
        self.process()
        try:
            n_alpha = densities.to("1 / meter**3")
        except AttributeError as error:
            raise TypeError("densities must be a unit-bearing quantity.") from error
        if np.asarray(n_alpha.magnitude).shape != (len(self),):
            raise ValueError("densities must contain one value per species.")
        if not np.all(np.isfinite(n_alpha.magnitude)) or np.any(n_alpha.magnitude < 0):
            raise ValueError("densities must be finite and non-negative.")
        mu_independant_scattering = np.einsum("i,i -> ", self.Csca, n_alpha).to("1/meter")

        return mu_independant_scattering

    def get_mu_dependant(self, densities: NDArray, H: NDArray, wavenumber: NDArray, theta_points: int = 150):
        """
        Compute the dependent scattering contribution using inter particle correlations H(wavenumber).

        Parameters
        ----------
        densities : NDArray, shape (N,)
            Species number densities n_alpha. Units should be 1/length**3.
        H : NDArray, shape (N, N, Pp)
            Mixture total correlation function in reciprocal space, evaluated on the
            `wavenumber` grid. This should correspond to the same ordering as size classes.
        wavenumber : NDArray, shape (Pp,)
            Reciprocal space radial grid (1/length). Must cover the range needed to
            evaluate H at wavenumber = 2k sin(φ/2) for φ in [0, pi].
        theta_points
            Number of azimuthal samples for the integral over the rotation angle.

        Returns
        -------
        mu_dependant_scattering
            Scalar attenuation coefficient contribution. Units depend on the exact
            normalization conventions of F and H used in the integrand, but it is
            intended to be in 1/length.

        Notes
        -----
        The core integrand is:

        term[p_index, theta_index] = sum_{a,b} sqrt(n_a n_b) F_a(wavenumber,theta) F*_b(wavenumber,theta) H_ab(wavenumber)

        then integrated over polar angle (self.phi) and azimuth (theta).
        """
        self._validate_mixture_inputs(densities, H, wavenumber)
        n_alpha, sqrt_alpha_beta = self.get_alpha_beta_factor(densities=densities)

        F, theta = self.get_F_matrix(theta_points=theta_points)

        phi, interpolated_H = self.get_interpolated_H(H=H, wavenumber=wavenumber)

        term = np.einsum("ab,iapt,ibpt,abp -> pt", sqrt_alpha_beta, F, np.conj(F), interpolated_H)

        val = np.trapezoid(term.T * np.sin(phi), x=phi, axis=1)

        mu_dependant_scattering = np.trapezoid(val, x=theta)

        return mu_dependant_scattering

    def get_mu(self, densities: NDArray, H: NDArray, wavenumber: NDArray, theta_points: int = 150) -> NDArray:
        """
        Compute total scattering attenuation coefficient mu_s.

        Parameters
        ----------
        densities : NDArray, shape (N,)
            Species number densities.
        H : NDArray, shape (N, N, Pp)
            Correlation function in wavenumber space.
        wavenumber : NDArray, shape (Pp,)
            wavenumber space grid (1/length).
        theta_points
            Azimuthal sampling count.

        Returns
        -------
        mu
            Total scattering attenuation coefficient (intended units 1/length), computed as:

            mu = mu_independant + mu_dependant
        """
        mu_independant = self.get_mu_independant(densities=densities)

        mu_dependant = self.get_mu_dependant(densities=densities, H=H, wavenumber=wavenumber, theta_points=theta_points)

        return mu_independant + mu_dependant

    def interpolate_last_axis_linear(
        self,
        H: np.ndarray,
        wavenumber: np.ndarray,
        evaluation_wavenumber: np.ndarray,
    ) -> np.ndarray:
        """
        Linear interpolation of H(..., wavenumber) onto evaluation_wavenumber, along the last axis.

        Parameters
        ----------
        H : ndarray, shape (..., Pp)
            Values sampled on grid wavenumber along the last axis.
        wavenumber : ndarray, shape (Pp,)
            Strictly increasing sample points.
        evaluation_wavenumber : ndarray, shape (Pq,)
            Query points.

        Returns
        -------
        ndarray, shape (..., Pq)
            Interpolated values.
        """
        H = np.asarray(H)
        wavenumber = np.asarray(wavenumber, dtype=float)
        evaluation_wavenumber = np.asarray(evaluation_wavenumber, dtype=float)

        if wavenumber.ndim != 1 or wavenumber.size < 2:
            raise ValueError("wavenumber must be a 1D grid with at least two points.")
        if not np.all(np.isfinite(wavenumber)) or np.any(np.diff(wavenumber) <= 0):
            raise ValueError("wavenumber must be finite and strictly increasing.")
        if not np.all(np.isfinite(evaluation_wavenumber)):
            raise ValueError("evaluation_wavenumber must contain only finite values.")
        if evaluation_wavenumber.min() < wavenumber[0] or evaluation_wavenumber.max() > wavenumber[-1]:
            raise ValueError("wavenumber does not cover the scattering wavevector range.")
        if H.shape[-1] != wavenumber.size:
            raise ValueError(f"H last axis ({H.shape[-1]}) must match wavenumber size ({wavenumber.size}).")

        H_2d = H.reshape(-1, wavenumber.size)  # (M, Pp)
        out_2d = np.empty((H_2d.shape[0], evaluation_wavenumber.size), dtype=H_2d.dtype)

        for row in range(H_2d.shape[0]):
            out_2d[row] = np.interp(evaluation_wavenumber, wavenumber, H_2d[row])

        return out_2d.reshape(H.shape[:-1] + (evaluation_wavenumber.size,))

    def get_interpolated_H(self, H, wavenumber):
        """
        Interpolate H(wavenumber) onto the scattering wavevector magnitude wavenumber = 2k sin(φ/2).

        Parameters
        ----------
        H : array like, shape (N, N, Pp)
            Correlation tensor sampled on the input wavenumber grid.
        wavenumber : array like, shape (Pp,)
            Radial reciprocal space grid (1/length).

        Returns
        -------
        interpolated_H : NDArray, shape (N, N, Pphi)
            H evaluated at the scattering wavenumber for φ in [0, pi], where:

            wavenumber(φ) = 2 k sin(φ/2)

        Notes
        -----
        This method assumes:
        * `self.phi` spans [0, pi] and is the polar scattering angle
        * `wavenumber` is provided in increasing order
        """
        self._validate_mixture_inputs(
            densities=np.ones(len(self)) / ureg.meter**3,
            H=H,
            wavenumber=wavenumber,
        )
        phi_evaluate = self._validated_phi()

        evaluation_wavenumber = 2 * self.k * np.sin(phi_evaluate / 2)

        wavenumber_grid = wavenumber.to("1/meter").magnitude
        evaluation_wavenumber_magnitude = evaluation_wavenumber.to("1/meter").magnitude

        interpolated_H = self.interpolate_last_axis_linear(
            H,
            wavenumber_grid,
            evaluation_wavenumber_magnitude,
        )
        return phi_evaluate, interpolated_H

    def get_phase_function(
        self, densities: NDArray, H: NDArray, wavenumber: NDArray, theta_points: int = 150
    ) -> tuple[NDArray, NDArray, NDArray]:
        """
        Compute the mixture phase function including dependent scattering corrections.

        Parameters
        ----------
        densities : NDArray, shape (N,)
            Species number densities n_alpha.
        H : NDArray, shape (N, N, Pp)
            Correlation tensor in reciprocal space on the wavenumber grid.
        wavenumber : NDArray, shape (Pp,)
            wavenumber space grid (1/length).
        theta_points
            Number of azimuthal samples.

        Returns
        -------
        phi : ndarray, shape (Pphi,)
            Polar angle grid in radians.
        theta : ndarray, shape (Ptheta,)
            Azimuthal angle grid in radians.
        phase_function : NDArray, shape (Pphi, Ptheta)
            Phase function indexed first by polar angle, then by azimuth.

        Notes
        -----
        The phase function is assembled as:

        * Independent term: sum of ``n |F|^2`` over species.
        * Dependent term: weighted cross terms involving ``F``, its conjugate,
          and the correlation tensor ``H``.
        """
        self._validate_mixture_inputs(densities, H, wavenumber)
        n_alpha, sqrt_alpha_beta = self.get_alpha_beta_factor(densities=densities)

        phi, interpolated_H = self.get_interpolated_H(H=H, wavenumber=wavenumber)

        F, theta = self.get_F_matrix(theta_points=theta_points)

        phase_function: NDArray = (
            np.einsum("a, jatp, jatp -> tp", n_alpha, F, np.conj(F))
            + np.einsum("ab, japt, jbpt, abp -> pt", sqrt_alpha_beta, F, np.conj(F), interpolated_H)
        )

        return phi, theta, phase_function
