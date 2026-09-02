"""Shared physical dimensions and typography for manuscript figures."""

from matplotlib.axes import Axes


# The article's 6.5-inch text block scales these 130%-width source figures
# down on inclusion. Typography is deliberately 150% larger at the source,
# yielding a modest, consistent increase in the printed figures.
MANUSCRIPT_WIDTH = 8.45
AXIS_LABEL_SIZE = 15
TICK_LABEL_SIZE = 12
LEGEND_FONT_SIZE = 10.5
PANEL_LABEL_SIZE = 21


def style_axes(axes: Axes, *, grid: bool = True) -> None:
    """Apply the shared manuscript typography to an axes."""
    if grid:
        axes.grid(alpha=0.18)
    axes.tick_params(labelsize=TICK_LABEL_SIZE)
    axes.xaxis.label.set_size(AXIS_LABEL_SIZE)
    axes.yaxis.label.set_size(AXIS_LABEL_SIZE)


def add_panel_label(axes: Axes, label: str, *, x: float = -0.20) -> None:
    """Place a bold panel label just outside the upper-left axes corner."""
    axes.annotate(
        label,
        xy=(x, 1.00),
        xycoords="axes fraction",
        xytext=(0, 6),
        textcoords="offset points",
        fontsize=PANEL_LABEL_SIZE,
        fontweight="bold",
        va="bottom",
        ha="left",
        color="#102a43",
        annotation_clip=False,
        zorder=10,
    )
