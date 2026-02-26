from typing import Literal, Tuple, Dict
from functools import cached_property
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from MPSPlots import helper
from TypedUnit.units import ureg

def _minimum_image_displacement(delta: np.ndarray, box_length: float) -> np.ndarray:
    return delta - box_length * np.round(delta / box_length)


class Result():
    """
    Output container for an RSA simulation run.

    Holds arrays plus domain metadata, computed statistics, and plotting helpers.
    """
    def __init__(self, binding):
        """
        Initialize the Result object.

        Parameters
        ----------
        binding : Binding
            The binding object containing simulation data.
        """
        self.binding = binding

    @cached_property
    def positions(self) -> np.ndarray:
        """
        Get the sphere center positions as a NumPy array of shape (N, 3).

        Returns
        -------
        np.ndarray
            Array of sphere center positions.
        """
        return self.sphere_configuration.positions

    @cached_property
    def radii(self) -> np.ndarray:
        """
        Get the sphere radii as a NumPy array of shape (N,).

        Returns
        -------
        np.ndarray
            Array of sphere radii.
        """
        return self.sphere_configuration.radii

    @cached_property
    def statistics(self):
        """
        Get the simulation statistics.

        Returns
        -------
        Statistics
            The simulation statistics.
        """
        return self.binding.statistics

    @cached_property
    def sphere_configuration(self):
        """
        Get the sphere configuration object.

        Returns
        -------
        SphereConfiguration
            The configuration of spheres in the simulation.
        """
        return self.binding.sphere_configuration

    @cached_property
    def domain(self):
        """
        Get the simulation domain.

        Returns
        -------
        Domain
            The domain (box) dimensions and boundary conditions.
        """
        return self.binding.domain

    def compute_partial_pair_correlation_function(self, **kwargs) -> tuple:
        """
        Compute the partial pair correlation function g_ij(r).

        Parameters
        ----------
        **kwargs : dict
            Keyword arguments passed to the binding method.

        Returns
        -------
        tuple
            A tuple containing:
            - centers (Quantity): The radial distances with units.
            - g_ij (np.ndarray): The partial pair correlation values.
        """
        centers, g_ij = self.binding.compute_partial_pair_correlation_function(**kwargs)
        return centers * ureg.meter, g_ij

    @cached_property
    def partial_volume_fractions(self) -> np.ndarray:
        """
        Get the partial volume fractions of each component.

        Returns
        -------
        np.ndarray
            Array of partial volume fractions.
        """
        return np.asarray(self.binding.partial_volume_fractions)

    @cached_property
    def partial_volumes(self) -> np.ndarray:
        """
        Get the partial volumes of each component.

        Returns
        -------
        np.ndarray
            Array of partial volumes.
        """
        return np.asarray(self.binding.partial_volumes)

    def compute_pair_correlation_function(self, **kwargs) -> None:
        """
        Compute the total pair correlation function g(r).

        Parameters
        ----------
        **kwargs : dict
            Keyword arguments passed to the binding method.
        """
        return self.binding.compute_pair_correlation_function(**kwargs)

    @cached_property
    def pair_correlation_centers(self) -> np.ndarray:
        """
        Get the radial centers for the pair correlation function.

        Returns
        -------
        np.ndarray
            Array of radial centers.
        """
        return np.asarray(self.binding.pair_correlation_centers)

    @cached_property
    def pair_correlation_values(self) -> np.ndarray:
        """
        Get the values of the pair correlation function g(r).

        Returns
        -------
        np.ndarray
            Array of pair correlation values.
        """
        return np.asarray(self.binding.pair_correlation_values)


    @helper.post_mpl_plot
    def plot_centers_3d(self, maximum_points_3d: int = 10_000) -> plt.Figure:
        """
        Plot the sphere centers in a 3D scatter plot.

        Parameters
        ----------
        maximum_points_3d : int
            Maximum number of points to plot (subsampling if necessary).

        Returns
        -------
        plt.Figure
            The matplotlib Figure object containing the 3D scatter plot.
        """
        n = self.positions.shape[0]
        random_generator = np.random.default_rng()

        if n > maximum_points_3d:
            selected = random_generator.choice(n, size=maximum_points_3d, replace=False)
        else:
            selected = np.arange(n)

        figure, axes = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(8, 6))

        axes.scatter(
            self.positions[selected, 0],
            self.positions[selected, 1],
            self.positions[selected, 2],
            s=6,
            alpha=0.6,
        )
        axes.set_xlim(0, self.domain.length_x)
        axes.set_ylim(0, self.domain.length_y)
        axes.set_zlim(0, self.domain.length_z)
        axes.set_xlabel("x")
        axes.set_ylabel("y")
        axes.set_zlabel("z")
        axes.set_title("RSA centers (subsampled)")


        return figure

    @helper.post_mpl_plot
    def plot_radius_distribution(self, bins: int = 40, density: bool = True, alpha: float = 0.85) -> plt.Figure:
        """
        Plot the distribution of sphere radii.

        Parameters
        ----------
        bins : int
            Number of histogram bins.
        density : bool
            Whether to normalize the histogram to form a probability density.
        alpha : float
            Transparency level for the histogram bars.
        """
        radii = self.radii

        figure, axes = plt.subplots()

        axes.hist(radii, bins=bins, density=density, alpha=alpha)
        axes.set_xlabel("radius")
        axes.set_ylabel("density" if density else "count")
        axes.set_title("Radius distribution")

        return figure

    def _compute_slice_mask(self, coord: np.ndarray, slice_center: float, slice_thickness: float, box_length: float) -> np.ndarray:
        """
        Compute the mask for particles within the slice thickness.
        """
        if self.domain.use_periodic_boundaries:
            delta = _minimum_image_displacement(coord - slice_center, box_length)
            return np.abs(delta) <= 0.5 * slice_thickness
        else:
            return np.abs(coord - slice_center) <= 0.5 * slice_thickness

    @helper.post_mpl_plot
    def plot_slice_2d(self, slice_axis: Literal["x", "y", "z"] = "z", slice_center_fraction: float = 0.5, slice_thickness_fraction: float = 0.08, maximum_circles_in_slice: int = 2500) -> plt.Figure:
        """
        Plot a 2D slice of the sphere configuration.

        Parameters
        ----------
        slice_axis : Literal["x", "y", "z"]
            Axis along which to take the slice.
        slice_center_fraction : float
            Fractional position along the slice axis where the slice is centered (0.0 to 1.0).
        slice_thickness_fraction : float
            Fractional thickness of the slice relative to the box length along the slice axis (0.0 to 1.0).
        maximum_circles_in_slice : int
            Maximum number of circles to plot in the slice (subsampling if necessary).
        """
        if not (0.0 <= slice_center_fraction <= 1.0):
            raise ValueError("slice_center_fraction must be between 0.0 and 1.0")
        if not (0.0 <= slice_thickness_fraction <= 1.0):
            raise ValueError("slice_thickness_fraction must be between 0.0 and 1.0")

        axis_map = {"x": 0, "y": 1, "z": 2}
        slice_idx = axis_map[slice_axis]

        # Determine plot axes indices (a, b) such that slice_idx is excluded
        # For 'z' (2) -> x(0), y(1)
        # For 'y' (1) -> x(0), z(2)
        # For 'x' (0) -> y(1), z(2)
        remaining_axes = [i for i in range(3) if i != slice_idx]
        a_idx, b_idx = remaining_axes[0], remaining_axes[1]

        labels = ["x", "y", "z"]
        a_label, b_label = labels[a_idx], labels[b_idx]

        # box_lengths = np.array([self.domain.length_x, self.domain.length_y, self.domain.length_z])
        box_lengths = [self.domain.length_x, self.domain.length_y, self.domain.length_z]

        # Dimensions for plot
        a_max = box_lengths[a_idx]
        b_max = box_lengths[b_idx]

        # Slice parameters
        L_slice = box_lengths[slice_idx]
        slice_center = slice_center_fraction * L_slice
        slice_thickness = slice_thickness_fraction * L_slice

        # Compute mask
        coord = self.positions[:, slice_idx]
        slice_mask = self._compute_slice_mask(coord, slice_center, slice_thickness, L_slice)

        slice_positions = self.positions[slice_mask]
        slice_radii = self.radii[slice_mask]

        if getattr(a_max, 'to', None):
             a_max_val = a_max.to('meter').magnitude
             b_max_val = b_max.to('meter').magnitude
             slice_pos_val = slice_positions.to('meter').magnitude
             slice_radii_val = slice_radii.to('meter').magnitude
        else:
             a_max_val = a_max
             b_max_val = b_max
             slice_pos_val = slice_positions
             slice_radii_val = slice_radii

        # Subsampling
        n_spheres = slice_positions.shape[0]
        if n_spheres > maximum_circles_in_slice:
            rng = np.random.default_rng()
            indices = rng.choice(n_spheres, size=maximum_circles_in_slice, replace=False)
            slice_pos_val = slice_pos_val[indices]
            slice_radii_val = slice_radii_val[indices]

        figure, axes = plt.subplots()
        axes.set_aspect("equal", adjustable="box")
        axes.set_xlim(0, a_max_val)
        axes.set_ylim(0, b_max_val)
        axes.set_xlabel(a_label)
        axes.set_ylabel(b_label)
        axes.set_title(
            f"2D slice at {slice_axis}≈{slice_center:.2f}, thickness {slice_thickness:.2f}"
            + f" | showing {int(np.sum(slice_mask))} spheres"
        )

        for center, radius in zip(slice_pos_val, slice_radii_val):
            axes.add_patch(
                Circle(
                    (center[a_idx], center[b_idx]),
                    radius=radius,
                    fill=False,
                    linewidth=0.8,
                    alpha=0.7,
                )
            )

        axes.plot([0, a_max_val, a_max_val, 0, 0], [0, 0, b_max_val, b_max_val, 0], linewidth=1.2, color='black')

        return figure

    @helper.post_mpl_plot
    def plot_pair_correlation(self, n_bins: int = 80, maximum_pairs: int = 300_000) -> plt.Figure:
        """
        Plot the partial pair correlation functions g_ij(r) obtained from the RSA
        configuration, overlaying all (i, j) curves on a single axis.

        Parameters
        ----------
        n_bins : int
            Number of radial distance bins.
        maximum_pairs : int
            Number of Monte Carlo sampled pairs used for estimation.

        Returns
        -------
        matplotlib.figure.Figure
            Figure containing the overlaid plot of all g_ij(r) curves.
        """
        # Call C++ to compute (centers, g_ij)
        centers, g_matrix = self.binding.compute_partial_pair_correlation_function(
            n_bins=n_bins,
            maximum_pairs=maximum_pairs,
        )

        centers = np.asarray(centers)
        g_matrix = np.asarray(g_matrix)  # shape (K, K, bins)
        K = g_matrix.shape[0]

        figure, ax = plt.subplots(1, 1)

        for i in range(K):
            for j in range(K):
                gij = g_matrix[i, j]
                ax.plot(
                    centers,
                    gij,
                    linewidth=1.3,
                    alpha=0.9,
                    label=f"g[{i},{j}](r)",
                )

        ax.set_xlabel("r")
        ax.set_ylabel("g_ij(r)")
        ax.grid(alpha=0.2)

        ax.set_title(
            f"Partial pair correlation functions g_ij(r), "
            f"K={K}, periodic={self.domain.use_periodic_boundaries}",
            fontsize=12,
            weight="bold",
        )

        # Legend can get large when K grows, so keep it compact
        ax.legend(
            fontsize=8,
            ncol=min(K * K, 4),
            frameon=False,
        )

        return figure
