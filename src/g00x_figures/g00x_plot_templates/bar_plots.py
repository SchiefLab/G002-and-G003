import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import ticker as mtick


def stacked_bar_plot_pivot_df(
    pivot_df: pd.DataFrame,
    title: str,
    ylabel: str | None = None,
    xlabel: str | None = None,
    yaxis_mtick_major: float | None = None,
    yaxis_mtick_minor: float | None = None,
    remove_yticklabels: bool = False,
    xticklabels: list[str] | None = None,
    xticklabel_rotation: int = 0,
    remove_legend: bool = True,
    yaxis_mtick_normalize: bool = True,
    colors: list[str] | None = None,
    ax: plt.Axes | None = None,  # type: ignore
    *args,  # type: ignore
    **kwargs,  # type: ignore
) -> plt.Axes:  # -> Any | Axes | None:
    """stacked bar plot from a pivot dataframe.
    should be limited indexes and columns to keep size down.

    Parameters
    ----------
    pivot_df : pd.DataFrame
        dataframe with multiindex and columns
    title : str
        title of the plot
    ylabel : str | None, optional
        yaxis label, by default None
    xlabel : str | None, optional
        xaxis label, by default None
    yaxis_mtick_major : float | None, optional
        physical yaxis tick - combines with minor, by default None
    yaxis_mtick_minor : float | None, optional
        physical yaxis tick - combines with major, by default None
    remove_yticklabels : bool, optional
        remove values for each ytick, by default False
    xticklabels : list[str] | None, optional
        xtick values, by default None
    xticklabel_rotation : int, optional
        xtick value rotation - 0 is horizonal; 90 is vertical, by default 0
    remove_legend : bool, optional
        removes legen incase you have a shared legend you need to custom make later, by default True
    yaxis_mtick_normalize : bool, optional
        i.e. make ytick value .5 to 50 if False, by default True
    colors : list[str] | None, optional
        hex or color name list with exact number of elements as xaxis, by default None
    ax : plt.Axes | None, optional
        matplotlib graph index location to edit specific subgraph, by default None

    Returns
    -------
    plt.Axes
        matplotlib graph index location to edit specific subgraph
    >>> stacked_bar_plot_pivot_df(
            pivot_df=g00x_seq_prime_igg_vrc01_class_pivot_trial_df,
            title="VRC01-class",
            ylabel="% BCRs with VRC01-class \n" + r"bnAb $\mathregular{V_{K/L}}$",
            yaxis_mtick_major=0.2,
            yaxis_mtick_minor=0.1,
            yaxis_mtick_normalize=False,
            colors=v_call_top_gene_color_mapping.values(),
            ax=axes[0],)
    >>> stacked_bar_plot_pivot_df(
            pivot_df=deskosky_ready,
            title="Control",
            xticklabels=["Dekosky\nVH1-2"],
            remove_yticklabels=True,
            yaxis_mtick_major=0.2,
            yaxis_mtick_minor=0.1,
            colors=v_call_top_gene_color_mapping.values(),
            ax=axes[1],)
    """
    ax: plt.Axes = pivot_df.plot(
        kind="bar",
        stacked=True,
        color=colors,
        width=0.85,
        linewidth=1,
        edgecolor="black",
        ax=ax,
        *args,
        **kwargs,
    )
    ax.set_title(title)
    ax.set_ylabel(ylabel or "")
    ax.set_xlabel(xlabel or "")
    ax.set_xticklabels(xticklabels or ax.get_xticklabels(), rotation=xticklabel_rotation)
    if remove_legend:
        ax.legend_.remove()  # type: ignore
    if yaxis_mtick_major:
        ax.yaxis.set_major_locator(mtick.MultipleLocator(yaxis_mtick_major))
    if yaxis_mtick_minor:
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(yaxis_mtick_minor))
    if remove_yticklabels:
        ax.set_yticklabels([])
    if not yaxis_mtick_normalize:
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
    return ax
