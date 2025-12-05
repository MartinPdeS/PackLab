from typing import Literal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from MPSPlots import helper

from PackLab.utils import _minimum_image_displacement
from PackLab.binary.interface_domain import Domain
from PackLab.binary.interface_result import Result



class Result(Result):
    """
    Output container for an RSA simulation run.

    Holds arrays plus domain metadata, computed statistics, and plotting helpers.
    """
    def __init__(self, positions: np.ndarray, radii: np.ndarray, domain: Domain, statistics):
        """
        Initialize the Result object with positions, radii, domain, and statistics.

        Parameters
        ----------
        positions : np.ndarray
            Array of shape (N, 3) containing the sphere center positions.
        radii : np.ndarray
            Array of shape (N,) containing the sphere radii.
        domain : Domain
            The simulation domain.
        statistics : Statistics
            The statistics object containing simulation metrics.
        """
        super().__init__(
            positions=positions,
            radii=radii,
            domain=domain
        )

        if self.positions.ndim != 2 or self.positions.shape[1] != 3:
            raise ValueError(f"positions must have shape (N, 3), got {self.positions.shape}")
        if self.radii.ndim != 1:
            raise ValueError(f"radii must have shape (N,), got {self.radii.shape}")
        if self.positions.shape[0] != self.radii.shape[0]:
            raise ValueError("positions and radii must have matching N")


        self._random_generator = np.random.default_rng

        self.statistics = statistics

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

        random_generator = self._random_generator()
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
        pair_correlation_bins: int = 90,
        maximum_number_of_pairs: int = 2_000_000,
        maximum_distance: float = 0.0,
    ) -> plt.Figure:
        """
        Plot the radial pair correlation function g(r) using particles in the
        simulation domain.

        Parameters
        ----------
        pair_correlation_bins : int
            Number of histogram bins for the radial distance.
        pair_correlation_maximum_pairs : int
            Maximum number of random particle pairs used to approximate g(r).
        maximum_distance : float
            Maximum distance to consider for g(r). If zero or negative, it defaults to half the smallest box length.

        Returns
        -------
        matplotlib.figure.Figure
            The Matplotlib figure containing the plot.
        """
        self.compute_pair_correlation_function(
            bins=pair_correlation_bins,
            maximum_number_of_pairs=maximum_number_of_pairs,
            random_seed=0,
            maximum_distance=maximum_distance
        )

        figure, axes = plt.subplots()

        axes.plot(self.pair_correlation_centers, self.pair_correlation_values, linewidth=1.6)

        axes.set_xlabel("radial distance r")
        axes.set_ylabel("g(r)")
        axes.set_title(
            "Pair correlation function"
            + (" (minimum image)" if self.domain.use_periodic_boundaries else "")
        )

        return figure

    @helper.post_mpl_plot
    def plot_pair_correlation(
        self,
        pair_correlation_bins: int = 90,
        maximum_number_of_pairs: int = 2_000_000,
        maximum_distance: float = 0.0,
        plot_mean_and_std: bool = False,
        repeats: int = 8,
        random_seed: int = 0,
    ) -> plt.Figure:
        """
        Plot the radial pair correlation function g(r).

        If plot_mean_and_std is True, the pair correlation function is computed
        multiple times and the mean and standard deviation are plotted. The
        standard deviation is shown as a shaded band.

        Parameters
        ----------
        pair_correlation_bins : int
            Number of histogram bins for the radial distance.
        maximum_number_of_pairs : int
            Maximum number of random particle pairs used per repetition.
        maximum_distance : float
            Maximum distance to consider for g(r). If zero or negative, it defaults
            to half the smallest box length.
        plot_mean_and_std : bool
            If True, compute g(r) repeats times and plot mean and standard deviation.
            If False, compute and plot a single g(r).
        repeats : int
            Number of repetitions used when plot_mean_and_std is True.
        random_seed : int
            Base random seed used for reproducible Monte Carlo sampling.

        Returns
        -------
        matplotlib.figure.Figure
            The Matplotlib figure containing the plot.
        """
        if plot_mean_and_std:
            self.compute_pair_correlation_function_mean_and_std(
                bins=pair_correlation_bins,
                maximum_number_of_pairs=maximum_number_of_pairs,
                repeats=repeats,
                random_seed=random_seed,
                maximum_distance=maximum_distance,
            )
            centers = self.pair_correlation_centers
            mean_values = self.pair_correlation_mean_values
            std_values = self.pair_correlation_std_values
        else:
            self.compute_pair_correlation_function(
                bins=pair_correlation_bins,
                maximum_number_of_pairs=maximum_number_of_pairs,
                random_seed=random_seed,
                maximum_distance=maximum_distance,
            )
            centers = self.pair_correlation_centers
            mean_values = self.pair_correlation_values
            std_values = None

        figure, axes = plt.subplots()

        axes.plot(centers, mean_values, linewidth=1.6)

        if std_values is not None:
            lower = mean_values - std_values
            upper = mean_values + std_values
            axes.fill_between(centers, lower, upper, alpha=0.25)

        axes.set_xlabel("radial distance r")
        axes.set_ylabel("g(r)")
        axes.set_title(
            "Pair correlation function"
            + (" (minimum image)" if self.domain.use_periodic_boundaries else "")
        )

        return figure
