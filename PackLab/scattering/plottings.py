import matplotlib.pyplot as plt
import numpy as np
from typing import Literal

def _coerce_phase_function_axes(
    phase_function: np.ndarray,
    phi: np.ndarray,
    theta: np.ndarray,
) -> np.ndarray:
    """
    Ensure phase_function has shape (phi.size, theta.size).

    Accepts:
    - (phi, theta)
    - (theta, phi) and transposes
    """
    phase_function = np.asarray(phase_function)

    if phase_function.shape == (phi.size, theta.size):
        return phase_function
    if phase_function.shape == (theta.size, phi.size):
        return phase_function.T

    raise ValueError(
        "phase_function has incompatible shape. "
        f"Got {phase_function.shape}, expected ({phi.size}, {theta.size}) or ({theta.size}, {phi.size})."
    )


def plot_phase_function_3d(
    phi: np.ndarray,
    theta: np.ndarray,
    phase_function: np.ndarray,
    *,
    mode: Literal["spherical", "surface"] = "spherical",
    normalize: bool = True,
    use_magnitude: bool = True,
) -> plt.Figure:
    """
    Plot a 3D representation of the phase function P(phi, theta).

    Parameters
    ----------
    phi
        Polar scattering angle grid in radians, expected in [0, pi].
    theta
        Azimuthal angle grid in radians, expected in [0, 2*pi].
    phase_function
        Phase function sampled on (phi, theta) or (theta, phi).
    mode
        "spherical": map intensity to radius on a sphere and render a 3D surface.
        "surface": plot a 3D surface with axes (theta, phi, P).
    normalize
        If True, normalize the plotted values by their maximum for visual clarity.
    use_magnitude
        If True, plot abs(P). If False, plot real(P).

    Returns
    -------
    figure
        Matplotlib figure containing the 3D plot.

    Notes
    -----
    In "spherical" mode, the surface is parameterized as:
        x = r * sin(phi) * cos(theta)
        y = r * sin(phi) * sin(theta)
        z = r * cos(phi)
    where r is the (optionally normalized) phase function.
    """
    phi = np.asarray(phi)
    theta = np.asarray(theta)

    P = _coerce_phase_function_axes(phase_function=phase_function, phi=phi, theta=theta)

    if use_magnitude:
        values = np.abs(P)
    else:
        values = np.real(P)

    values = np.asarray(values, dtype=float)

    if normalize:
        max_value = np.nanmax(values)
        if max_value > 0:
            values = values / max_value

    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="xy")

    figure = plt.figure()
    ax = figure.add_subplot(111, projection="3d")

    if mode == "surface":
        ax.plot_surface(theta_grid, phi_grid, values, linewidth=0, antialiased=True)
        ax.set_xlabel("theta (rad)")
        ax.set_ylabel("phi (rad)")
        ax.set_zlabel("P")
        ax.set_title("Phase function surface P(phi, theta)")
        return figure

    if mode == "spherical":
        r = np.clip(values, a_min=0.0, a_max=None)
        x = r * np.sin(phi_grid) * np.cos(theta_grid)
        y = r * np.sin(phi_grid) * np.sin(theta_grid)
        z = r * np.cos(phi_grid)

        ax.plot_surface(x, y, z, linewidth=0, antialiased=True)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.set_title("Phase function mapped on sphere (radius proportional to P)")
        ax.set_box_aspect([1, 1, 1])
        ax.set(
            xlim=[-1, 1],
            ylim=[-1, 1],
            zlim=[-1, 1],
        )
        return figure

    raise ValueError('mode must be "spherical" or "surface".')


def plot_phase_function_2d_projection(
    phi: np.ndarray,
    theta: np.ndarray,
    phase_function: np.ndarray,
    *,
    projection: Literal["azimuth_average", "heatmap"] = "azimuth_average",
    use_magnitude: bool = True,
    normalize: bool = False,
) -> plt.Figure:
    """
    Plot a 2D representation of the phase function.

    Parameters
    ----------
    phi
        Polar scattering angle grid in radians, expected in [0, pi].
    theta
        Azimuthal angle grid in radians, expected in [0, 2*pi].
    phase_function
        Phase function sampled on (phi, theta) or (theta, phi).
    projection
        "azimuth_average": plot P_avg(phi) = (1/2pi) * integral P(phi,theta) dtheta.
        "heatmap": plot a 2D image of P(phi,theta) with axes theta and phi.
    use_magnitude
        If True, plot abs(P). If False, plot real(P).
    normalize
        If True, normalize plotted values by their maximum.

    Returns
    -------
    figure
        Matplotlib figure containing the 2D plot.

    Notes
    -----
    The azimuth averaged projection is the closest analogue to typical S1 and S2
    plots because it reduces the 2D angular dependence into a 1D curve versus phi.
    """
    phi = np.asarray(phi)
    theta = np.asarray(theta)

    P = _coerce_phase_function_axes(phase_function=phase_function, phi=phi, theta=theta)

    if use_magnitude:
        values = np.abs(P)
    else:
        values = np.real(P)

    values = np.asarray(values, dtype=float)

    if normalize:
        max_value = np.nanmax(values)
        if max_value > 0:
            values = values / max_value

    figure, ax = plt.subplots(1, 1)

    if projection == "heatmap":
        image = ax.imshow(
            values,
            aspect="auto",
            origin="lower",
            extent=[theta.min(), theta.max(), phi.min(), phi.max()],
        )
        ax.set_xlabel("theta (rad)")
        ax.set_ylabel("phi (rad)")
        ax.set_title("Phase function heatmap P(phi, theta)")
        figure.colorbar(image, ax=ax, label="P")
        return figure

    if projection == "azimuth_average":
        P_phi = np.trapezoid(values, x=theta, axis=1) / (theta.max() - theta.min())
        ax.plot(phi, P_phi)
        ax.set_xlabel("phi (rad)")
        ax.set_ylabel("Azimuth averaged P(phi)")
        ax.set_title("Azimuth averaged phase function")
        return figure

    raise ValueError('projection must be "azimuth_average" or "heatmap".')
