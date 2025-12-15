from typing import Literal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from MPSPlots import helper

from PackLab.monte_carlo.utils import _minimum_image_displacement
from PackLab.binary.interface_domain import Domain
from PackLab.binary.interface_result import Result



class Result():
    """
    Output container for an RSA simulation run.

    Holds arrays plus domain metadata, computed statistics, and plotting helpers.
    """
    def __init__(self, binding):
        self.binding = binding

    @property
    def positions(self) -> np.ndarray:
        """
        Get the sphere center positions as a NumPy array of shape (N, 3).

        Returns
        -------
        np.ndarray
            Array of sphere center positions.
        """
        return self.sphere_configuration.positions_numpy()

    @property
    def radii(self) -> np.ndarray:
        """
        Get the sphere radii as a NumPy array of shape (N,).

        Returns
        -------
        np.ndarray
            Array of sphere radii.
        """
        return self.sphere_configuration.radii_numpy()

    @property
    def domain(self) -> Domain:
        """
        Get the simulation domain.

        Returns
        -------
        Domain
            The simulation domain.
        """
        return self.binding.domain

    @property
    def statistics(self):
        """
        Get the simulation statistics.

        Returns
        -------
        Statistics
            The simulation statistics.
        """
        return self.binding.statistics

    def compute_pair_correlation_function(self, **kwargs) -> None:
        return self.binding.compute_pair_correlation_function(**kwargs)

    @property
    def pair_correlation_centers(self) -> np.ndarray:
        return np.asarray(self.binding.pair_correlation_centers)

    @property
    def pair_correlation_values(self) -> np.ndarray:
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
        box_lengths = [self.domain.length_x, self.domain.length_y, self.domain.length_z]

        axis_to_index = {"x": 0, "y": 1, "z": 2}
        slice_axis_index = axis_to_index[slice_axis]

        if not (0.0 <= slice_center_fraction <= 1.0):
            raise ValueError("slice_center_fraction must be between 0.0 and 1.0")
        if not (0.0 <= slice_thickness_fraction <= 1.0):
            raise ValueError("slice_thickness_fraction must be between 0.0 and 1.0")

        slice_center = slice_center_fraction * box_lengths[slice_axis_index]
        slice_thickness = slice_thickness_fraction * box_lengths[slice_axis_index]

        coord = self.positions[:, slice_axis_index]
        if self.domain.use_periodic_boundaries:
            delta = _minimum_image_displacement(coord - slice_center, box_lengths[slice_axis_index])
            slice_mask = np.abs(delta) <= 0.5 * slice_thickness
        else:
            slice_mask = np.abs(coord - slice_center) <= 0.5 * slice_thickness

        slice_positions = self.positions[slice_mask]
        slice_radii = self.radii[slice_mask]

        if slice_axis == "z":
            a_index, b_index = 0, 1
            a_label, b_label = "x", "y"
            a_max, b_max = self.domain.length_x, self.domain.length_y
        elif slice_axis == "y":
            a_index, b_index = 0, 2
            a_label, b_label = "x", "z"
            a_max, b_max = self.domain.length_x, self.domain.length_z
        else:
            a_index, b_index = 1, 2
            a_label, b_label = "y", "z"
            a_max, b_max = self.domain.length_y, self.domain.length_z

        random_generator = np.random.default_rng()
        if slice_positions.shape[0] > maximum_circles_in_slice:
            chosen = random_generator.choice(slice_positions.shape[0], size=maximum_circles_in_slice, replace=False)
            slice_positions = slice_positions[chosen]
            slice_radii = slice_radii[chosen]

        figure, axes = plt.subplots()

        axes.set_aspect("equal", adjustable="box")
        axes.set_xlim(0, a_max)
        axes.set_ylim(0, b_max)
        axes.set_xlabel(a_label)
        axes.set_ylabel(b_label)
        axes.set_title(
            f"2D slice at {slice_axis}≈{slice_center:.2f}, thickness {slice_thickness:.2f}"
            + f" | showing {int(np.sum(slice_mask))} spheres"
        )

        for center, radius in zip(slice_positions, slice_radii):
            axes.add_patch(
                Circle(
                    (center[a_index], center[b_index]),
                    radius=radius,
                    fill=False,
                    linewidth=0.8,
                    alpha=0.7,
                )
            )

        axes.plot([0, a_max, a_max, 0, 0], [0, 0, b_max, b_max, 0], linewidth=1.2)

        return figure

    @helper.post_mpl_plot
    def plot_pair_correlation(
        self,
        bins: int = 80,
        maximum_pairs: int = 300_000,
    ) -> plt.Figure:
        """
        Plot the partial pair correlation functions g_ij(r) obtained from the RSA
        configuration. Produces a K by K panel grid, where K is the number of size
        classes defined by the radius sampler.

        Parameters
        ----------
        number_of_distance_bins : int
            Number of radial distance bins.
        maximum_pairs : int
            Number of Monte Carlo sampled pairs used for estimation.

        Returns
        -------
        matplotlib.figure.Figure
            Figure containing the grid plot of g_ij(r).
        """

        # Call C++ to compute (centers, g_ij)
        centers, g_matrix = self.binding.compute_partial_pair_correlation_function(
            number_of_distance_bins=bins,
            maximum_pairs=maximum_pairs,
        )

        centers = np.asarray(centers)
        g_matrix = np.asarray(g_matrix)  # shape (K, K, bins)
        K = g_matrix.shape[0]

        # Create the panel grid
        figure, axes = plt.subplots(
            nrows=K,
            ncols=K,
            sharex=True,
            sharey=True,
        )

        # If K == 1, axes is not a 2D array
        if K == 1:
            axes = np.array([[axes]])

        for i in range(K):
            for j in range(K):
                ax = axes[i, j]
                gij = g_matrix[i, j]

                ax.plot(centers, gij, color='black', linewidth=1.3)

                # Axis labels on left column and bottom row
                if j == 0:
                    ax.set_ylabel(f"g[{i},{j}](r)")
                if i == K - 1:
                    ax.set_xlabel("r")

                # Light grid
                ax.grid(alpha=0.2)

        figure.suptitle(
            f"Partial pair correlation functions g_ij(r)\n"
            f"K={K} size classes, periodic={self.domain.use_periodic_boundaries}",
            fontsize=14,
            weight="bold",
        )
        figure.tight_layout(rect=[0, 0, 1, 0.96])

        return figure
