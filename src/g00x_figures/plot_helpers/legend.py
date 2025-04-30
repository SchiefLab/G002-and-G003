from typing import Tuple

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D


def plot_legend(
    pallete: dict[str, str],
    ax: plt.Axes,  # type: ignore
    n: int | None = None,
    loc: str = "lower center",
    bbox_to_anchor: Tuple[float, float] = (0.5, -0.15),
    handlelength: float = 0,
    labelspacing: float = 0.3,
    fontsize: int = 16,
    *args,
    **kwargs,
) -> None:
    """Plot a legend with custom circle markers.
    Designed for: O G001 O G002 O G003

    Parameters
    ----------
    pallete : dict[str, str]
        A dictionary with the labels and colors to be used in the legend.
    ax : plt.Axes
        Matplotlib axis object.
    n : int, optional
        Number of columns in the legend, by default None
    loc : str, optional
        Matplotlib axis location, by default "lower center"
    bbox_to_anchor : Tuple[int, int], optional
        Legend origin to keep it outside of plot, by default (0.5, -0.15)
        The default position is below the plot.
    fontsize : int, optional
        Marker and label font size, by default 14
    args : Any
        Matplotlib legend arguments
    kwargs : Any
        Matplotlib legend keyword arguments

    Returns
    -------
    None
        the legend is plotted via object reference.
    """
    custom_lines = []
    for x in pallete:
        custom_lines.append(
            Line2D(
                [0],
                [0],
                marker="s",  # Change marker to square
                color="w",
                markerfacecolor=pallete[x],
                markeredgecolor="black",
                markersize=fontsize - 2,
                label=x,
            )
        )
    ax.legend(
        custom_lines,
        [i.get_label() for i in custom_lines],
        loc=loc,
        frameon=False,
        handlelength=handlelength,
        ncol=len(pallete.keys()) if n is None else n,
        bbox_to_anchor=bbox_to_anchor,
        labelspacing=labelspacing,
        *args,
        **kwargs,
        fontsize=fontsize,
    )
