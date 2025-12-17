from dataclasses import dataclass
from typing import Optional, Sequence, Literal
import numpy as np
from numpy import pi
from TypedUnit import ureg

from MPSPlots import helper
from tabulate import tabulate


@dataclass
class Domain:
    """
    Cubic simulation domain containing a polydisperse population of spherical particles.

    You must specify the domain ``size`` and the particle radii array ``particle_radii``.
    You can then specify the population using exactly one of the following mechanisms:

    1) ``number_fractions``
       Fractions per radius bin that sum to 1 (dimensionless). Total occupied volume is set
       by ``volume_fraction``.

    2) ``volume_fractions``
       Volume fraction per radius bin. If provided, it must sum to ``volume_fraction``.

    3) ``particle_densities``
       Number density per radius bin (for example in 1 / liter). In this case, ``volume_fraction``
       can be provided for consistency checks but is not required.

    Notes:
        1. This class mirrors the monodisperse ``Domain`` interface, but returns vector valued
           quantities for per bin results and provides totals where applicable.
        2. Integer particle counts are inherently a discretization. This class provides deterministic
           rounding and a small correction to preserve the total count when ``number_fractions`` is used.

    Parameters
    ----------
    size : object
        Domain size (length of one cubic side).
    particle_radii : object
        Array of particle radii defining the polydisperse bins.
    volume_fraction : Optional[object]
        Total occupied volume fraction (between 0 and 1). Required if ``number_fractions`` or ``volume_fractions`` are used.
    number_fractions : Optional[Sequence[float]]
        Number fractions per radius bin that sum to 1 (dimensionless). Mutually exclusive with ``volume_fractions`` and ``particle_densities``.
    rounding_mode : Literal["floor", "round"]
        Rounding mode used when converting non integer particle counts to integers.
        "floor" always rounds down, "round" uses standard rounding rules.
    """
    size: object
    particle_radii: object
    volume_fraction: Optional[object] = None
    number_fractions: Optional[Sequence[float]] = None

    rounding_mode: Literal["floor", "round"] = "floor"

    def __post_init__(self) -> None:
        """Validate inputs and normalize or reshape fields.

        Raises
        ------
        ValueError
            If input validation fails.
        """

        # self.size = 1 * ureg.micrometer
        self.particle_radii = np.asarray(self.particle_radii.magnitude) * self.particle_radii.units


        if self.particle_radii.ndim != 1:
            raise ValueError("particle_radii must be a one dimensional array.")

        self.particle_densities = None
        self.volume_fractions = None

        self.number_fractions = np.asarray(self.number_fractions, dtype=float)

        if self.number_fractions.shape != (self.particle_radii.size,):
            raise ValueError("number_fractions must have the same length as particle_radii.")

        if np.any(self.number_fractions < 0.0):
            raise ValueError("number_fractions must be non negative.")

        total = float(np.sum(self.number_fractions))
        if total <= 0.0:
            raise ValueError("number_fractions must sum to a positive value.")

        self.number_fractions = self.number_fractions / total

        if self.rounding_mode not in ("floor", "round"):
            raise ValueError("rounding_mode must be either 'floor' or 'round'.")

    @property
    def volume(self):
        """Total domain volume.

        Returns
        -------
        float
            Domain volume.
        """
        return self.size ** 3

    @property
    def particle_volumes(self):
        """Per bin particle volumes.

        Returns
        -------
        np.ndarray
            Array of per bin particle volumes.
        """
        return (4.0 / 3.0) * pi * (self.particle_radii ** 3)

    @property
    def total_particle_volume(self):
        """Total occupied particle volume.

        If ``volume_fraction`` is provided, the total occupied volume is ``volume_fraction * volume``.
        If ``particle_densities`` are used, the occupied volume is inferred from the particle counts.

        Returns
        -------
        float
            Total occupied particle volume.
        """
        return self.volume_fraction * self.volume

    @property
    def mean_particle_volume_number_weighted(self):
        """Number weighted mean particle volume.

        Returns
        -------
        float
            Mean particle volume weighted by number fractions.

        Raises
        ------
        AttributeError
            If number_fractions is not defined.
        """
        if self.number_fractions is None:
            raise AttributeError(
                "mean_particle_volume_number_weighted is defined only when number_fractions is used."
            )
        return np.sum(self.number_fractions * self.particle_volumes)

    @property
    def number_of_particles_total(self) -> int:
        """Total number of particles in the domain, inferred from the chosen specification.

        Returns
        -------
        int
            Total integer particle count.

        Raises
        ------
        RuntimeError
            If the internal specification state is invalid.
        """
        if self.particle_densities is not None:
            counts = self.particle_densities * self.volume
            return int(np.sum(self._apply_rounding(counts)))

        total_occupied_volume = self.volume_fraction * self.volume

        if self.volume_fractions is not None:
            counts = (self.volume_fractions * self.volume) / self.particle_volumes
            return int(np.sum(self._apply_rounding(counts)))

        if self.number_fractions is not None:
            mean_volume = self.mean_particle_volume_number_weighted
            estimated_total = total_occupied_volume / mean_volume
            return int(self._apply_rounding(estimated_total))

        raise RuntimeError("Invalid specification state.")

    @property
    def number_of_particles_per_radius(self):
        """Integer particle count per radius bin.

        For ``number_fractions``, a deterministic correction is applied so the per bin counts
        sum exactly to ``number_of_particles_total``.

        Returns
        -------
        np.ndarray
            Array of integer particle counts per radius bin.

        Raises
        ------
        RuntimeError
            If the internal specification state is invalid.
        """
        if self.particle_densities is not None:
            counts = self.particle_densities * self.volume
            return self._apply_rounding(counts).astype(int)

        if self.volume_fractions is not None:
            counts = (self.volume_fractions * self.volume) / self.particle_volumes
            return self._apply_rounding(counts).astype(int)

        if self.number_fractions is not None:
            total_count = self.number_of_particles_total
            raw_counts = self.number_fractions * total_count

            counts = self._apply_rounding(raw_counts).astype(int)
            difference = int(total_count - np.sum(counts))

            if difference != 0:
                residuals = raw_counts - counts
                if difference > 0:
                    indices = np.argsort(residuals)[::-1]
                    counts[indices[:difference]] += 1
                else:
                    indices = np.argsort(residuals)
                    removable = np.where(counts > 0)[0]
                    indices = [i for i in indices if i in set(removable)]
                    counts[np.array(indices[: abs(difference)], dtype=int)] -= 1

            return counts

        raise RuntimeError("Invalid specification state.")

    @property
    def particle_densities_per_radius(self):
        """Number density per radius bin.


        Returns
        -------
        np.ndarray
            Array of per bin number densities.
        """
        return self.number_of_particles_per_radius / self.volume

    @property
    def particle_density_total(self):
        """Total number density over all bins.

        Returns
        -------
        float
            Total number density.
        """
        return self.number_of_particles_total / self.volume

    @property
    def volume_fraction_per_radius(self):
        """Per bin occupied volume fraction.

        Returns
        -------
        np.ndarray
            Array of per bin volume fractions.
        """
        occupied_volume_per_radius = self.number_of_particles_per_radius * self.particle_volumes
        return occupied_volume_per_radius / self.volume

    def _apply_rounding(self, values):
        """
        Apply configured rounding mode.

        Parameters
        ----------
        values : np.ndarray
            Array of values to round.

        Returns
        -------
        np.ndarray
            Rounded values.
        """
        if self.rounding_mode == "floor":
            return np.floor(values)
        return np.round(values)

    def sample_particle_radii(self, number_of_samples: int, random_generator: Optional[np.random.Generator] = None):
        """Draw radii samples according to number fractions.

        Parameters
        ----------
        number_of_samples : int
            Number of radii samples to draw.
        random_generator : Optional[np.random.Generator]
            Optional random generator to use for sampling. If None, a new default generator is created.

        Returns
        -------
        np.ndarray
            Array of sampled particle radii.
        """
        if self.number_fractions is None:
            raise ValueError("Sampling radii requires number_fractions.")
        if random_generator is None:
            random_generator = np.random.default_rng()

        indices = random_generator.choice(
            a=self.particle_radii.size,
            size=int(number_of_samples),
            replace=True,
            p=self.number_fractions,
        )
        return self.particle_radii[indices]

    @staticmethod
    def _format_value(value: object, precision: int = 6) -> str:
        """Format a value for human readable printing.

        This is a best effort formatter that works for Pint quantities, numpy scalars,
        and plain python scalars.

        Parameters
        ----------
        value : object
            The value to format.
        precision : int
            Number of significant digits to use.

        Returns
        -------
        str
            The formatted string.
        """
        if value is None:
            return "None"

        try:
            value_compact = value.to_compact()
        except Exception:
            value_compact = value

        try:
            magnitude = float(value_compact.magnitude)
            units = str(value_compact.units)
            return f"{magnitude:.{precision}g} {units}".strip()
        except Exception:
            pass

        try:
            return f"{float(value):.{precision}g}"
        except Exception:
            return str(value)

    def print_bins(self, precision: int = 6) -> None:
        """Print a per bin table describing the polydisperse population.

        The table includes radii, per particle volumes, per bin number fractions
        (if defined), per bin volume fractions, inferred integer particle counts,
        and inferred number densities.

        Parameters
        ----------
        precision : int
            Number of significant digits to use when printing values.

        Returns
        -------
        None
            Prints the bin table to standard output.
        """
        n_bins = int(self.particle_radii.size)

        indices = np.arange(n_bins, dtype=int)

        radii = self.particle_radii[indices]
        particle_volumes = self.particle_volumes[indices]
        number_of_particles = self.number_of_particles_per_radius[indices]
        particle_densities = self.particle_densities_per_radius[indices]
        volume_fractions = self.volume_fraction_per_radius[indices]

        if self.number_fractions is not None:
            number_fractions = self.number_fractions[indices]
        else:
            number_fractions = [None for _ in indices]

        bin_rows = []
        for k, idx in enumerate(indices):
            bin_rows.append(
                [
                    int(idx),
                    self._format_value(radii[k], precision=precision),
                    self._format_value(particle_volumes[k], precision=precision),
                    self._format_value(number_fractions[k], precision=precision) if number_fractions[k] is not None else "None",
                    self._format_value(volume_fractions[k], precision=precision),
                    int(number_of_particles[k]),
                    self._format_value(particle_densities[k], precision=precision),
                ]
            )

        print(
            tabulate(
                bin_rows,
                headers=[
                    "bin",
                    "radius",
                    "particle_volume",
                    "number_fraction",
                    "volume_fraction",
                    "particle_count",
                    "number_density",
                ],
            )
        )

    @helper.post_mpl_plot
    def plot_radius_distribution(self, normalize: bool = True, show_volume_weighted: bool = False) -> None:
        """
        This is a diagnostic plot intended to quickly verify that the binning and
        mixture specification match your expectations.

        The plot uses the current domain specification:

        - If ``number_fractions`` is defined, it is plotted directly.
        - Otherwise the plotted distribution is inferred from the current per bin particle counts (``number_of_particles_per_radius``) and normalized.


        Optionally, you can overlay a volume weighted view (each bin weighted by
        particle volume), which is often useful when comparing number weighted and
        volume weighted interpretations of the same mixture.

        Parameters
        ----------
        normalize : bool
            Whether to normalize the plotted distribution to sum to 1.
        show_volume_weighted : bool
            Whether to overlay a volume weighted version of the distribution.

        Returns
        -------
        plt.Figure
            Displays a matplotlib figure.

        """
        import matplotlib.pyplot as plt

        radii = self.particle_radii
        x = radii.magnitude

        if self.number_fractions is not None:
            y = np.asarray(self.number_fractions, dtype=float)
            label = "Number fractions (input)"
        else:
            counts = np.asarray(self.number_of_particles_per_radius, dtype=float)
            y = counts
            label = "Number fractions (inferred from counts)"

        if normalize:
            total = float(np.sum(y))
            if total > 0:
                y = y / total

        figure, ax = plt.subplots(1, 1)

        ax.bar(x, y, width=np.diff(x).min() if x.size > 1 else 1.0, align="center", alpha=0.6, label=label, edgecolor='black')


        if show_volume_weighted:
            weights = np.asarray((self.particle_volumes / self.particle_volumes.max()).to("").magnitude, dtype=float)

            if self.number_fractions is not None:
                y_vol = np.asarray(self.number_fractions, dtype=float) * weights
                vol_label = "Volume weighted (from number_fractions)"
            else:
                counts = np.asarray(self.number_of_particles_per_radius, dtype=float)
                y_vol = counts * weights
                vol_label = "Volume weighted (from counts)"

            if normalize:
                total_vol = float(np.sum(y_vol))
                if total_vol > 0:
                    y_vol = y_vol / total_vol

            ax.plot(x, y_vol, marker="s", linestyle="--", label=vol_label)

        ax.set_xlabel(f"Radius [{radii.units:~P}]")
        ax.set_ylabel("Normalized weight" if normalize else "Weight")
        ax.set_title("Discretized radius distribution")
        ax.legend()

        return figure
