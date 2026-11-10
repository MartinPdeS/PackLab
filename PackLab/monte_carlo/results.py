from functools import cached_property
from functools import wraps
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from TypedUnit.units import ureg


def post_mpl_plot(plot_function):
    """Return a plotting wrapper with the former ``show`` convenience option."""

    @wraps(plot_function)
    def wrapper(*args, show: bool = True, **kwargs):
        figure = plot_function(*args, **kwargs)
        if show:
            plt.show()
        return figure

    return wrapper


class PackingResult:
    """
    Output container for an RSA or Metropolis simulation run.

    Holds arrays plus domain metadata, computed statistics, and plotting helpers.
    """

    def __init__(self, binding, source: str | None = None, run_metadata: dict | None = None):
        """
        Initialize the packing result.

        Parameters
        ----------
        binding : Binding
            The binding object containing simulation data.
        source : {"rsa", "metropolis"} or None, optional
            Workflow that produced this result when known.
        run_metadata : dict or None, optional
            Reproducibility-relevant options retained by the producing workflow.
        """
        self.binding = binding
        self.source = source or "unknown"
        self.run_metadata = dict(run_metadata or {})

    def save(self, path, **kwargs) -> None:
        """
        Save this configuration in PackLab's portable ``.npz`` format.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination archive ending in ``.npz``.
        **kwargs
            Additional metadata accepted by :func:`save_packing`.
        """
        from PackLab.monte_carlo.persistence import save_packing

        save_packing(self, path, **kwargs)

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
            ``n_bins``, ``maximum_pairs`` (zero examines all pairs), and
            ``random_seed`` (default zero, deterministic) passed to the binding.
            Capped pairs are sampled uniformly without replacement.

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

    @post_mpl_plot
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

        positions = self.positions
        if hasattr(positions, "to"):
            positions = positions.to("meter").magnitude

        lengths = (self.domain.length_x, self.domain.length_y, self.domain.length_z)
        lengths = [
            length.to("meter").magnitude if hasattr(length, "to") else length for length in lengths
        ]

        axes.scatter(
            positions[selected, 0],
            positions[selected, 1],
            positions[selected, 2],
            s=6,
            alpha=0.6,
        )
        axes.set_xlim(0, lengths[0])
        axes.set_ylim(0, lengths[1])
        axes.set_zlim(0, lengths[2])
        axes.set_xlabel("x (m)")
        axes.set_ylabel("y (m)")
        axes.set_zlabel("z (m)")
        axes.set_title("RSA centers (subsampled)")

        return figure

    @post_mpl_plot
    def plot_radius_distribution(
        self, bins: int = 40, density: bool = True, alpha: float = 0.85
    ) -> plt.Figure:
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

    def plot_packing_3d(self, *, show: bool = True) -> go.Figure:
        """Show all spheres and the simulation box in an interactive 3D view.

        Parameters
        ----------
        show : bool, default=True
            Open the Plotly figure after constructing it.

        Returns
        -------
        plotly.graph_objects.Figure
            A figure in micrometres, with sphere radii represented geometrically.
            Spheres are shown whole at their stored centres; periodic images are
            not drawn, so particles at a periodic boundary can extend outside the box.

        Notes
        -----
        Every particle is rendered as a triangulated sphere. Large packings can
        require substantial browser memory and take time to render.
        """
        positions = self.positions
        radii = self.radii
        positions = np.asarray(
            positions.to("micrometer").magnitude if hasattr(positions, "to") else positions * 1e6
        )
        radii = np.asarray(
            radii.to("micrometer").magnitude if hasattr(radii, "to") else radii * 1e6
        )
        lengths = np.asarray(
            [
                length.to("micrometer").magnitude if hasattr(length, "to") else length * 1e6
                for length in (self.domain.length_x, self.domain.length_y, self.domain.length_z)
            ]
        )
        longitude = np.linspace(0, 2 * np.pi, 13)[:-1]
        latitude = np.linspace(0, np.pi, 9)
        theta, phi = np.meshgrid(longitude, latitude)
        unit_vertices = np.column_stack(
            (
                (np.sin(phi) * np.cos(theta)).ravel(),
                (np.sin(phi) * np.sin(theta)).ravel(),
                np.cos(phi).ravel(),
            )
        )
        vertices_per_sphere = len(unit_vertices)
        vertices = (positions[:, None, :] + radii[:, None, None] * unit_vertices).reshape(-1, 3)
        triangles = []
        for row in range(len(latitude) - 1):
            for column in range(len(longitude)):
                first = row * len(longitude) + column
                second = row * len(longitude) + (column + 1) % len(longitude)
                triangles.extend(
                    (
                        (first, second, first + len(longitude)),
                        (second, second + len(longitude), first + len(longitude)),
                    )
                )
        triangles = np.asarray(triangles)
        triangles = (
            triangles[None, :, :] + np.arange(len(radii))[:, None, None] * vertices_per_sphere
        ).reshape(-1, 3)
        figure = go.Figure()
        if len(radii):
            figure.add_trace(
                go.Mesh3d(
                    x=vertices[:, 0],
                    y=vertices[:, 1],
                    z=vertices[:, 2],
                    i=triangles[:, 0],
                    j=triangles[:, 1],
                    k=triangles[:, 2],
                    color="#2583b8",
                    opacity=1.0,
                    flatshading=False,
                    lighting=dict(ambient=0.55, diffuse=0.8, specular=0.25),
                    hoverinfo="skip",
                    name="Particles",
                )
            )
        corners = np.array(
            [[x, y, z] for x in (0, lengths[0]) for y in (0, lengths[1]) for z in (0, lengths[2])]
        )
        edges = [(i, j) for i in range(8) for j in range(i + 1, 8) if bin(i ^ j).count("1") == 1]
        coordinates = [
            [value for i, j in edges for value in (corners[i, axis], corners[j, axis], None)]
            for axis in range(3)
        ]
        figure.add_trace(
            go.Scatter3d(
                x=coordinates[0],
                y=coordinates[1],
                z=coordinates[2],
                mode="lines",
                line=dict(color="#263746", width=4),
                hoverinfo="skip",
                showlegend=False,
            )
        )
        figure.update_layout(
            scene=dict(
                xaxis_title="x (µm)", yaxis_title="y (µm)", zaxis_title="z (µm)", aspectmode="data"
            ),
            margin=dict(l=0, r=0, t=35, b=0),
            title=f"Packing: {len(radii)} spheres",
        )
        if show:
            figure.show()
        return figure

    @post_mpl_plot
    def plot_pair_correlation(
        self, n_bins: int = 80, maximum_pairs: int = 300_000, random_seed: int = 0
    ) -> plt.Figure:
        """
        Plot the partial pair correlation functions g_ij(r) obtained from the RSA
        configuration, overlaying all (i, j) curves on a single axis.

        Parameters
        ----------
        n_bins : int
            Number of radial distance bins.
        maximum_pairs : int
            Maximum number of unordered pairs sampled uniformly without replacement.
            Zero examines all pairs.
        random_seed : int, default=0
            Seed for pair sampling; zero is deterministic.

        Returns
        -------
        matplotlib.figure.Figure
            Figure containing the overlaid plot of all g_ij(r) curves.
        """
        # Call C++ to compute (centers, g_ij)
        centers, g_matrix = self.binding.compute_partial_pair_correlation_function(
            n_bins=n_bins,
            maximum_pairs=maximum_pairs,
            random_seed=random_seed,
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
