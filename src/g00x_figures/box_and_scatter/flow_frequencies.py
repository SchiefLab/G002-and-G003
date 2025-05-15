import logging
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import patchworklib as pw
import seaborn as sns
from matplotlib import ticker as mtick
from matplotlib.patches import Patch

from g00x_figures.data import Data
from g00x_figures.plot_helpers.font import apply_global_font_settings
from g00x_figures.plot_helpers.legend import plot_legend as legend

# init logger
logger = logging.getLogger()
apply_global_font_settings()


class MinorSymLogLocator(mtick.Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.
    """

    def __init__(self, linthresh):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically
        """
        self.linthresh = linthresh

    def __call__(self):
        "Return the locations of the ticks"
        majorlocs = self.axis.get_majorticklocs()

        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in range(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i - 1]
            if abs(majorlocs[i - 1] + majorstep / 2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i - 1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minorlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError("Cannot get tick locations for a " "%s type." % type(self))


def format_y_axis(y: float, lt: float):
    """Format Y axis for log scale."""
    if y > lt:
        return "{:g}".format(y)
    else:
        if y == 0:
            return "0"
        if y == lt:
            return r"$\leq$" + r"$10^{" + str(round(math.log(y, 10))) + r"}$"
        else:
            return r"$10^{" + str(round(math.log(y, 10))) + r"}$"


def adjust_boxplot(ax: plt.Axes, median_linewidth: float = 2) -> plt.Axes:
    """For box and whisker plots, adjust the appearance of the boxplot.

    Parameters
    ----------
    ax : plt.axis

    Returns
    -------
    plt.axis
        return back axis
    """
    # for boxplots we can make adjustments
    list_of_groups = zip(*(iter(ax.get_lines()),) * 6)
    for line_set in list_of_groups:
        for line_index, line in enumerate(line_set):
            # the median line of the box plot
            if line_index == 4:
                line.set_color("black")
                line.set_linewidth(median_linewidth)
                line.set_solid_capstyle("butt")

    for mybox in ax.artists:
        # Change the appearance of the surrounding boxes, like the sides
        mybox.set_edgecolor("black")
        mybox.set_linewidth(0.4)
    return ax


def generic_plot_box_and_whisker(
    ax: plt.Axes,
    g002_flow_and_seq: pd.DataFrame,
    g001_flow_and_seq: pd.DataFrame,
    g002_label: str,
    g001_label: str,
    y_label: str,
    lt: float = 1e-4,
    pal={"G001": "#6C65FF", "G002": "#91FCC0"},
    y_min: float = 8e-5,
    y_max: float | int = 10,
    plot_legend: bool = True,
    adjust_data_to_lt: bool = True,
    plot_log: str = "log",
    plot_xticklabels: bool = True,
    ylabelpad: int = 0,
    strip_size: int = 5,
) -> plt.Axes:
    # combine g001 and g002 data axis
    plot_df = (
        pd.concat(
            [
                g002_flow_and_seq,
                g001_flow_and_seq.rename({g001_label: g002_label}, axis=1),
            ]
        )
        .astype({"weeks": int})
        .reset_index(drop=True)
    )

    data_for_strip = plot_df.copy()
    data_for_box = plot_df.copy()

    if adjust_data_to_lt:
        # any data less than
        di = plot_df[plot_df[g002_label] < lt].index

        # any data less than, just make it that lt number
        data_for_strip.loc[di, g002_label] = lt

        # but for box, just remove it
        data_for_box.loc[di, g002_label] = np.nan

    sns.stripplot(
        x="weeks",
        y=g002_label,
        hue="trial",
        edgecolor="black",
        linewidth=1,
        s=strip_size,
        jitter=0.1,  # type: ignore
        palette=pal,
        dodge=True,
        hue_order=["G001", "G002"],
        data=data_for_strip,
        ax=ax,
    )

    sns.boxplot(
        x="weeks",
        y=g002_label,
        hue="trial",
        dodge=True,
        palette=pal,
        fliersize=0,
        hue_order=["G001", "G002"],
        data=data_for_box,
        whis=[10, 90],
        ax=ax,
    )
    ax.legend_.set_visible(False)  # type: ignore
    ax.set_ylim(y_min, y_max)
    if plot_log == "log":
        ax.set_yscale("log")

    elif plot_log == "symlog":
        ax.set_yscale("symlog", linthresh=10)

    else:
        ax.set_yscale("linear")

    adjust_boxplot(ax)
    ax.set_ylabel(y_label, labelpad=ylabelpad)
    sns.despine()

    if plot_log == "log":
        ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: format_y_axis(y, lt)))
        ax.yaxis.set_major_locator(mtick.LogLocator(base=10, numticks=10))
        # ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
    elif plot_log == "symlog":
        ax.yaxis.set_minor_locator(MinorSymLogLocator(10))  # custom class to handle sym log ticks.
        ax.yaxis.set_major_formatter(mtick.ScalarFormatter())  # 10 and 100 are major ticks on a sym log
    else:
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())

    if plot_legend:
        custom_lines = []
        for x in pal:
            custom_lines.append(Patch(facecolor=pal[x], edgecolor="black", linewidth=1, label=x))
        ax.legend(
            custom_lines,
            [i._label for i in custom_lines],
            loc="upper center",
            frameon=False,
            handlelength=0.8,
            ncol=2,
            bbox_to_anchor=(0.5, -0.25),
            labelspacing=0.1,
        )
    if plot_xticklabels:
        ax.set_xticklabels(["-4/-5", "4", "8", "16", "24"])
        ax.set_xlabel("Weeks post vaccination")
    else:
        ax.set_xticklabels([])
        ax.set_xlabel("")
    return ax


def plot_flow_frequencies(data: Data, outpath: Path) -> None:
    """Plot flow frequencies."""
    logging.info("Plot flow frequeniencies")

    # G002 Flow and Seq
    g002_flow_and_seq = data.get_g002_flow_and_seq_prime()
    g001_flow_and_seq = data.get_g001_flow_and_seq_prime()

    figure, axes = plt.subplots(4, 2, figsize=(8.518, 8.5))

    # Plot 1 antigen specificity
    g002_label = "Percent antigen-specific among IgG^{+}"
    g001_label = "Percent of IgG+ B cells that are GT8++ (without regard to KO binding status)"
    y_label = "% GT8$^{++}$ among\nIgG$^+$ B cells"
    generic_plot_box_and_whisker(
        axes[0][0],
        g002_flow_and_seq,
        g001_flow_and_seq,
        g002_label,
        g001_label,
        y_label,
        plot_legend=False,
        plot_xticklabels=False,
    )

    # Plot 2 - Epitope Specificity
    g002_label = "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}"
    g001_label = "Percent of IgG+ B cells that are epitope-specific (KO-GT8++)"
    y_label = "% CD4bs-specific among\nIgG$^{+}$ B cells"
    generic_plot_box_and_whisker(
        axes[0][1],
        g002_flow_and_seq,
        g001_flow_and_seq,
        g002_label,
        g001_label,
        y_label,
        plot_legend=False,
        plot_xticklabels=False,
    )

    # Plot 3 specificity among KO
    g002_label = "Percent IgG^{+}KO^- among Ag^{++}"
    g001_label = "Percent of GT8++IgG+ B cells that are KO-"
    y_label = "% of GT8$^{++}$IgG$^{+}$\nB cells that are KO$^{-}$"
    generic_plot_box_and_whisker(
        axes[1][0],
        g002_flow_and_seq,
        g001_flow_and_seq,
        g002_label,
        g001_label,
        y_label,
        plot_legend=False,
        plot_xticklabels=False,
        y_min=0,
        y_max=105,
        plot_log="linear",
        ylabelpad=20,
    )

    # Plot 4 - Count VRC01 class
    g002_label = "Number of IGHG sequences that are VRC01-class"
    g001_label = "Number of epitope-specific (KO-GT8++) sequenced IgG BCRs that are VRC01-class"
    # y_label = "% of IgG$^{+}$ memory B cells\ndetected as VRC01-class"
    y_label = "Number of VRC01-class\nIgG$^{+}$ B cells"
    generic_plot_box_and_whisker(
        axes[1][1],
        g002_flow_and_seq,
        g001_flow_and_seq,
        g002_label,
        g001_label,
        y_label,
        plot_log="symlog",
        y_min=0,
        y_max=1300,
        plot_legend=False,
        plot_xticklabels=False,
        adjust_data_to_lt=False,
        ylabelpad=13,
    )

    # Plot 5 - VRC01 class
    g002_label = "Percent of VRC01-class sequences among IgG"
    g001_label = "Percent of IgG+ B cells detected as VRC01-class (missing seq to 0)"
    y_label = "% of IgG$^{+}$ memory B cells\ndetected as VRC01-class"
    generic_plot_box_and_whisker(
        axes[2][0],
        g002_flow_and_seq,
        g001_flow_and_seq,
        g002_label,
        g001_label,
        y_label,
        y_max=1,
        plot_legend=False,
        plot_xticklabels=False,
        lt=1e-5,
        y_min=8e-6,
    )

    # Plot 6 - Response
    combined = pd.concat(
        [
            g002_flow_and_seq,
            g001_flow_and_seq.rename({"Response (missing seq to 0)": "Response x"}, axis=1),
        ]
    )
    sns.pointplot(
        data=combined,
        hue="trial",
        y="Response x",
        x="weeks",
        join=False,
        ci="wilson",
        ax=axes[3][1],
        dodge=0.25,
        hue_order=["G001", "G002"],
        palette={"G001": "#6C65FF", "G002": "#91FCC0"},
    )
    axes[3][1].set_xlabel("Weeks post vaccination")
    axes[3][1].set_xticklabels(["-4/-5", "4", "8", "16", "24"])
    axes[3][1].legend_.remove()
    # axes[3][1].set_xticklabels([])
    axes[3][1].set_ylabel("% VRC01-class\nIgG$^{+}$ responders", labelpad=20)
    axes[3][1].yaxis.set_major_locator(mtick.MultipleLocator(base=0.2))
    axes[3][1].yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: str(int(float(y) * 100))))
    custom_lines = []
    pal = {"G001": "#6C65FF", "G002": "#91FCC0"}
    for x in pal:
        custom_lines.append(Patch(facecolor=pal[x], edgecolor="black", linewidth=1, label=x))
    axes[3][1].legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=2,
        bbox_to_anchor=(0.5, -0.25),
        labelspacing=0.1,
    )

    # Plot 7 - Efficiency
    a = "Percent of IGHG sequences that are VRC01-class"
    b = "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}"
    c = "Percent antigen-specific among IgG^{+}"
    d = "Percent VRC01-class among eOD-specific IgG+ memory BCR sequences"

    # this is how we derive D.
    g002_flow_and_seq[d] = g002_flow_and_seq[a] * (g002_flow_and_seq[b] / g002_flow_and_seq[c])

    g002_label = a
    g001_label = "Percent of epitope-specific (KO-GT8++) sequenced IgG BCRs that are VRC01-class"
    generic_plot_box_and_whisker(
        axes[2][1],
        g002_flow_and_seq,
        g001_flow_and_seq,
        g002_label,
        g001_label,
        "% VRC01-class among\nCD4bs-specific\nIgG$^{+}$ memory sequences",
        plot_xticklabels=False,
        plot_log="linear",
        y_max=104,
        y_min=0,
        plot_legend=False,
        ylabelpad=10,
    )
    # axes[1][1].legend(
    #     custom_lines,
    #     [i._label for i in custom_lines],
    #     loc="upper center",
    #     frameon=False,
    #     handlelength=0.8,
    #     ncol=2,
    #     bbox_to_anchor=(0.5, -0.25),
    #     labelspacing=0.1,
    # )
    # axes[0][0].legend(
    #     custom_lines,
    #     [i._label for i in custom_lines],
    #     loc="upper center",
    #     frameon=False,
    #     handlelength=0.8,
    #     ncol=2,
    #     bbox_to_anchor=(0.5, -0.25),
    #     labelspacing=0.1,
    # )
    # axes[1][0].legend(
    #     custom_lines,
    #     [i._label for i in custom_lines],
    #     loc="upper center",
    #     frameon=False,
    #     handlelength=0.8,
    #     ncol=2,
    #     bbox_to_anchor=(0.5, -0.25),
    #     labelspacing=0.1,
    # )
    # plot 8, antigen efficiency
    g002_label = d
    g001_label = "Percent of GT8++ IgG+ B cells detected as VRC01-class (missing seq to 0)"
    y_label = "% VRC01-class among\nGT8-specific\nIgG$^{+}$ memory sequences"
    generic_plot_box_and_whisker(
        axes[3][0],
        g002_flow_and_seq,
        g001_flow_and_seq,
        g002_label,
        g001_label,
        y_label,
        plot_xticklabels=True,
        plot_log="linear",
        y_max=104,
        y_min=0,
        plot_legend=True,
    )
    plt.tight_layout()

    figure.savefig(str(outpath) + ".png", dpi=300)


def plot_response_summary(data: Data, outpath: str):
    """plot the response summary between g001 and g002"""
    df = data.get_response_count_dataframe()
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(5.3, 2.6),
        sharey=True,
        gridspec_kw={"wspace": 0.1, "left": 0.15, "right": 0.99, "top": 0.83},
    )

    # first shot
    ax = axes[0]
    sns.barplot(
        x="trial",
        y="percentage_responders",
        data=df.query("shot=='first'"),
        linewidth=1,
        edgecolor="black",
        ax=ax,
        hue="trial",
        dodge=False,
        palette={"G001": "#6C65FF", "G002": "#91FCC0"},
    )
    ax.yaxis.set_major_locator(mtick.MultipleLocator(base=0.1))
    ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: str(int(float(y) * 100))))
    ax.set_ylim(0, 1)
    sns.despine()
    ax.legend_.remove()

    ax = axes[1]
    sns.barplot(
        x="trial",
        y="percentage_responders",
        data=df.query("shot=='all'"),
        linewidth=1,
        edgecolor="black",
        ax=ax,
        hue="trial",
        dodge=False,
        palette={"G001": "#6C65FF", "G002": "#91FCC0"},
    )
    sns.despine()
    ax.set_ylabel("")
    ax.legend_.remove()

    for ax_index, ax in enumerate(axes):
        for x in ax.get_xticklabels():
            (x1, _) = x.get_position()
            lookup = x.get_text()
            if ax_index == 0:
                value = "first"
            else:
                value = "all"
            sub_df = df.query("shot==@value").query("trial==@lookup")
            responders = sub_df["responders"].iloc[0]
            total = sub_df["total"].iloc[0]
            percentage = sub_df["percentage_responders"].iloc[0]
            ax.text(x1, percentage + 0.01, f"{responders}/{total}", ha="center")

    axes[0].set_title("After 1 vaccination\n(wk 4,8)", pad=10)
    axes[1].set_title("After 1 or 2 vaccinations\n(wk 4,8,16,24)", pad=10)
    axes[0].set_ylabel("% VRC01-class\nIgG$^{+}$ responders", labelpad=3)
    axes[0].set_xlabel("")
    axes[1].set_xlabel("")
    fig.savefig(str(outpath) + ".png", dpi=300)


def plot_boost_frequences(data: Data, outpath: Path, fig_num: str = "fig4") -> None:
    apply_global_font_settings(fontsize=14)
    boosting_data = data.get_g002_flow_and_seq_boost()
    week_colors = data.get_week_palette()
    figure, axes = plt.subplots(4, 2, figsize=(8.5, 11))
    axes = [
        [pw.Brick(figsize=(4, 2)), pw.Brick(figsize=(4, 2))],
        [pw.Brick(figsize=(4, 2)), pw.Brick(figsize=(4, 2))],
        [pw.Brick(figsize=(4, 2)), pw.Brick(figsize=(4, 2))],
        [pw.Brick(figsize=(4, 2)), pw.Brick(figsize=(4, 2))],
    ]
    a = "Percent of IGHG sequences that are VRC01-class"
    b = "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}"
    c = "Percent antigen-specific among IgG^{+}"
    d = "Percent VRC01-class among Core-specific IgG+ memory BCR sequences"

    outpath_img = outpath / "Main"
    outpath_metric = outpath / f"Main-Metrics/{fig_num}"

    # outpath_img.mkdir(parents=True, exist_ok=True)
    outpath_metric.mkdir(parents=True, exist_ok=True)

    # this is how we derive D.
    boosting_data[d] = boosting_data[a] * (boosting_data[b] / boosting_data[c])
    boosting_data = boosting_data[~boosting_data["pseudogroup"].isna()]

    def _plot_boost_frequency_figure(
        g002_label: str,
        y_label: str,
        ax: plt.Axes,
        letter: str,
        low: float = 0.005,
        high: float = 1,
        plot_log: str = "log",
        plot_lt: bool = False,
        lt: float = 1e-5,
        ylabelpad=3,
        use_threshold: bool = False,
    ) -> pd.DataFrame:
        data_for_strip = boosting_data.copy(deep=True)
        data_for_box = boosting_data.copy(deep=True)

        if plot_lt:
            # any data less than
            di = boosting_data[boosting_data[g002_label] < lt].index

            # any data less than, just make it that lt number
            data_for_strip.loc[di, g002_label] = lt

            # but for box, just remove it
            data_for_box.loc[di, g002_label] = np.nan

        # print(data_for_box.groupby("pseudogroup")[g002_label].median())
        median_df = data_for_box.groupby("pseudogroup")[g002_label].median()

        sns.stripplot(
            x="pseudogroup",
            y=g002_label,
            data=data_for_strip,
            hue="weeks",
            edgecolor="black",
            linewidth=1,
            size=7,
            dodge=False,
            palette=week_colors,
            ax=ax,
        )

        g = sns.boxplot(
            x="pseudogroup",
            y=g002_label,
            data=data_for_box,
            hue="weeks",
            # showmeans=True,
            linewidth=1,
            whis=[10, 90],
            fliersize=0,
            dodge=False,
            palette=week_colors,
            ax=ax,
        )
        data_for_box = boosting_data.copy(deep=True)
        yname = y_label.replace("\n", " ").strip().replace("%", "percent")
        data_for_box[["pubID", "trial", "pseudogroup", "weeks", g002_label]].to_csv(
            outpath_metric / f"fig{letter}_{yname}.csv", index=False
        )
        # from IPython import embed

        # # embed()
        # print(y_label)
        median_lines = [line for line in g.get_lines() if line.get_linestyle() == "-"]
        # for line in median_lines:
        #     # The x and y data for median lines
        #     xdata, ydata = line.get_data()
        #     # Median value is constant, so we can use the first ydata point
        #     median_value = ydata[0]

        #     # Calculate the x position. This assumes that median lines are vertical.
        #     x_position = np.mean(xdata)
        #     print(median_value, x_position)
        #     plt.text(
        #         x_position,
        #         median_value,
        #         f"{median_value:.2f}",
        #         # ha="center",
        #         # va="bottom",
        #         fontweight="bold",
        #         color="black",
        #         fontsize=6,
        #         transform=ax.transAxes,
        #     )
        ax.legend_.remove()
        adjust_boxplot(ax)
        if plot_log == "log":
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: format_y_axis(y, lt)))
            ax.yaxis.set_major_locator(mtick.LogLocator(base=10, numticks=10))
        elif plot_log == "symlog":
            ax.set_yscale("symlog", linthresh=10)
            ax.yaxis.set_minor_locator(MinorSymLogLocator(10))  # custom class to handle sym log ticks.
            ax.yaxis.set_major_formatter(mtick.ScalarFormatter())  # 10 and 100 are major ticks on a sym log
        else:
            ax.yaxis.set_major_locator(mtick.MultipleLocator(25))
            ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        ax.set_ylim(low, high)
        ax.xaxis.set_ticklabels([])
        ax.set_xlabel("")
        ax.set_ylabel(y_label, labelpad=ylabelpad)
        return median_df

    median_dfs = []

    # antigen specific
    # TODO: create csv
    g002_label = "Percent antigen-specific among IgG^{+}"
    y_label = "% core$^{++}$ among\nIgG$^+$ B cells"
    ax = axes[0][0]
    median_df = _plot_boost_frequency_figure(g002_label, y_label, ax, ylabelpad=20, letter="A")
    # add 'A_' to series column
    median_df.name = "figA"
    median_dfs.append(median_df)

    # TODO: create csv
    # Plot 2 - Epitope Specificity
    g002_label = "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}"
    y_label = "% CD4bs-specific among\nIgG$^{+}$ B cells"
    ax = axes[0][1]
    median_df = _plot_boost_frequency_figure(g002_label, y_label, ax, letter="B")
    median_df.name = "figB"
    median_dfs.append(median_df)

    # TODO: create csv
    # Plot 3 specificity among KO
    g002_label = "Percent IgG^{+}KO^- among Ag^{++}"
    y_label = "% of core$^{++}$IgG$^{+}$\nB cells that are KO$^{-}$"
    ax = axes[1][0]
    median_df = _plot_boost_frequency_figure(
        g002_label,
        y_label,
        ax,
        letter="F",
        low=0,
        high=101,
        plot_log="linear",
        ylabelpad=22,
    )
    median_df.name = "figF"
    median_dfs.append(median_df)

    # TODO: create csv
    # Plot 4, count VRC01 class among sequences
    g002_label = "Number of IGHG sequences that are VRC01-class"
    y_label = "Number of core$^{++}$KO$^{-}$\nVRC01-class IgG$^{+}$ B cells"
    ax = axes[1][1]
    median_df = _plot_boost_frequency_figure(g002_label, y_label, ax, low=0, high=1000, plot_log="symlog", letter="C")
    median_df.name = "figC"
    median_dfs.append(median_df)

    # TODO: create csv
    # plot 5, VRC01 among cells
    g002_label = "Percent of VRC01-class sequences among IgG"
    # y_label = "% of IgG$^{+}$ memory B cells\ndetected as VRC01-class"
    y_label = "% core$^{++}$KO$^{-}$\nVRC01-class B cells\namong IgG$^{+}$ memory B cells"
    ax = axes[2][0]
    median_df = _plot_boost_frequency_figure(
        g002_label,
        y_label,
        ax,
        letter="D",
        plot_lt=True,
        low=8e-6,
        lt=1e-5,
        high=1,
        plot_log="log",
    )
    median_df.name = "figD"
    median_dfs.append(median_df)

    # TODO: create csv
    # plot 6, response
    y_label = "% of participants\nwith core$^{++}$KO$^{-}$\nVRC01-class IgG$^{+}$ B cells"
    name = y_label.replace("\n", " ").strip().replace("%", "percent")
    boosting_data = boosting_data[~boosting_data[d].isna()]
    boosting_data[name] = boosting_data["Response y"]
    boosting_data[["pubID", "trial", "pseudogroup", "weeks", name]].to_csv(
        outpath_metric / f"figE_{name}.csv", index=False
    )
    sns.pointplot(
        data=boosting_data,
        hue="weeks",
        y="Response y",
        x="pseudogroup",
        join=False,
        ci="wilson",
        ax=axes[2][1],
        # dodge=0.25,
        palette=week_colors,
    )
    median_df = boosting_data.groupby("pseudogroup")[name].median()
    median_df.name = "figE"

    # removed medians for pointplot
    # median_dfs.append(median_df)
    axes[2][1].set_xlabel("")
    axes[2][1].set_ylim(bottom=0)
    axes[2][1].set_xticklabels([])
    axes[2][1].legend_.remove()
    axes[2][1].set_ylabel(
        "% of participants\nwith core$^{++}$KO$^{-}$\nVRC01-class IgG$^{+}$ B cells",
        labelpad=22,
    )
    axes[2][1].yaxis.set_major_locator(mtick.MultipleLocator(base=0.1))
    axes[2][1].yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: str(int(float(y) * 100))))

    # TODO: create csv
    g002_label = "Percent of IGHG sequences that are VRC01-class"
    y_label = "% VRC01-class among\nCD4bs-specific IgG$^{+}$\nmemory sequences"
    median_df = _plot_boost_frequency_figure(
        g002_label,
        y_label,
        axes[3][0],
        letter="G",
        low=-3,
        high=103,
        plot_log="linear",
        ylabelpad=10,
    )
    median_df.name = "figG"
    median_dfs.append(median_df)

    # TODO: create csv
    g002_label = "Percent VRC01-class among Core-specific IgG+ memory BCR sequences"
    y_label = "% VRC01-class among\ncore-specific IgG$^{+}$\nmemory sequences"
    median_df = _plot_boost_frequency_figure(
        g002_label,
        y_label,
        axes[3][1],
        low=-3,
        high=103,
        plot_log="linear",
        letter="H",
    )
    median_df.name = "figH"
    median_dfs.append(median_df)

    # set the bottom
    for ax in axes[3]:
        ax.set_xticklabels(
            [
                "eOD",
                "core",
                r"eOD$\rightarrow$eOD",
                r"eOD$\rightarrow$core",
                r"eOD$\rightarrow$eOD",
                r"eOD$\rightarrow$core",
                r"eOD$\rightarrow$eOD$\rightarrow$core",
            ],
            rotation=45,
            ha="right",
            rotation_mode="anchor",
        )
    custom_lines = []

    pal = {
        # 8: "grey",
        # 8: "#91FCC0",
        # 8: "#FFB3B2",
        8: "gold",
        16: "#E377C2",
        24: "#2078B4",
    }
    legend(
        ax=axes[3][0],
        pallete={"wk 8": "gold", "wk 16": "#E377C2", "wk 24": "#2078B4"},
        # fontsize=12,
        bbox_to_anchor=(0.5, -1),
    )
    legend(
        ax=axes[3][1],
        pallete={"wk 8": "gold", "wk 16": "#E377C2", "wk 24": "#2078B4"},
        # fontsize=12,
        bbox_to_anchor=(0.5, -1),
    )
    # flatten list of lists
    faxes = [item for sublist in axes for item in sublist]
    letters = "ABCDEFGH"
    letters = "ABFCDEGH"
    for ax, letter in zip(faxes, letters):
        ax.set_title(letter, fontsize=18, fontweight="bold", x=-0.3, y=1.1)

    for ax in faxes:
        sns.despine(ax=ax)
    # sns.despine()
    plt.tight_layout()
    median_df = pd.concat(median_dfs, axis=1)
    # median_df = median_df.round(3)
    median_df = median_df.reset_index()
    median_df = data.populate_psname(median_df)
    median_df.to_csv(outpath_metric / f"all_medians.csv", index=False)
    g1, g2, g3, g4, g5, g6, g7, g8 = faxes

    # g1 g6
    # g2 g4
    [
        [
            g1,
            g4,
            g2,
            g5,
            g3,
            g6,
            g7,
            g8,
        ]
    ]

    c1 = pw.stack([g1, g4, g6, g7], operator="/", margin=0.2)
    c2 = pw.stack([g2, g5, g3, g8], operator="/", margin=0.2)
    g = pw.hstack(c1, c2, margin=0.1)
    # g = (g1 | g2) / (g3 | g4) / (g5 | g6) / (g7 | g8)
    print(outpath_img)
    g.savefig(outpath_img / "fig4.png", dpi=700)


def plot_efficiency(data: Data, outpath: Path) -> None:
    """Plot efficiency."""
    boosting_data = data.get_g002_flow_and_seq_boost()
    week_colors = data.get_week_palette()
    figure, axes = plt.subplots(2, 1, figsize=(8.518 / 2, 7.5))

    a = "Percent of IGHG sequences that are VRC01-class"
    b = "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}"
    c = "Percent antigen-specific among IgG^{+}"
    d = "Percent VRC01-class among Core-specific IgG+ memory BCR sequences"

    # this is how we derive D.
    boosting_data[d] = boosting_data[a] * (boosting_data[b] / boosting_data[c])

    sns.stripplot(
        x="pseudogroup",
        y=a,
        data=boosting_data,
        hue="weeks",
        edgecolor="black",
        linewidth=1,
        size=7,
        dodge=False,
        palette=week_colors,
        ax=axes[0],
    )
    sns.boxplot(
        x="pseudogroup",
        y=a,
        data=boosting_data,
        hue="weeks",
        linewidth=1,
        whis=[10, 90],
        fliersize=0,
        dodge=False,
        palette=week_colors,
        ax=axes[0],
    )
    axes[0].legend_.remove()
    axes[0].set_ylabel("% VRC01-class among CD4bs-specific\nIgG+ memory BCR sequences")
    axes[0].set_xlabel("")
    axes[0].set_xticklabels([])
    adjust_boxplot(axes[0])

    sns.stripplot(
        x="pseudogroup",
        y=d,
        data=boosting_data,
        hue="weeks",
        edgecolor="black",
        linewidth=1,
        size=7,
        dodge=False,
        palette=week_colors,
        ax=axes[1],
    )
    sns.boxplot(
        x="pseudogroup",
        y=d,
        data=boosting_data,
        hue="weeks",
        linewidth=1,
        whis=[10, 90],
        fliersize=0,
        dodge=False,
        palette=week_colors,
        ax=axes[1],
    )
    adjust_boxplot(axes[1])
    axes[1].legend_.remove()
    axes[1].set_ylabel("% VRC01-class among core-specific\nIgG+ memory BCR sequences")
    axes[1].set_xticklabels(
        [
            "core",
            r"eOD$\rightarrow$eOD",
            r"eOD$\rightarrow$core",
            r"eOD$\rightarrow$eOD",
            r"eOD$\rightarrow$core",
            r"eOD$\rightarrow$eOD$\rightarrow$core",
        ],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
    )
    axes[1].set_xlabel("")
    custom_lines = []

    pal = {
        8: "#91FCC0",
        16: "#E377C2",
        24: "#2078B4",
    }
    for x in pal:
        custom_lines.append(Patch(facecolor=pal[x], edgecolor="black", linewidth=1, label=x))

    axes[1].legend(
        custom_lines,
        ["wk " + i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=3,
        bbox_to_anchor=(0.5, -0.4),
        labelspacing=0.1,
    )
    sns.despine()
    plt.tight_layout()

    figure.savefig(str(outpath) + ".png", dpi=300)


def plot_pre_post_frequency(data: Data):
    dfs = []

    apply_global_font_settings(fontsize=10)
    pre_post = data.get_pre_post_core_df()
    fig, axes = plt.subplots(3, 2, figsize=(8.5, 11 * 0.75))
    pal = data.get_week_palette()

    def clean_label(label: str) -> str:
        for i in ["\n", "^", "{", "}", "^"]:
            label = label.replace(i, "")
        label = label.replace("%", "percent")
        label = label.replace(" ", "_")
        return label

    def annotate_df(df: pd.DataFrame, ax: plt.Axes, key: list[Any], label: str) -> pd.DataFrame:
        df.name = ax.letter + "_" + clean_label(label)
        df.key = key
        return df

    def _plot_boost_frequency_figure(
        g002_label: str,
        y_label: str,
        ax: plt.Axes,
        low: float = 0.005,
        high: float = 1,
        plot_log: str = "log",
        plot_lt: bool = False,
        lt: float = 1e-4,
    ) -> pd.DataFrame:
        data_for_strip = pre_post.copy(deep=True)
        data_for_box = pre_post.copy(deep=True)

        if plot_lt:
            # any data less than
            di = pre_post[pre_post[g002_label] < lt].index

            # any data less than, just make it that lt number
            data_for_strip.loc[di, g002_label] = lt

            # but for box, just remove it
            data_for_box.loc[di, g002_label] = np.nan

        sns.stripplot(
            x="pseudogroup",
            y=g002_label,
            data=data_for_strip,
            edgecolor="black",
            hue="weeks",
            linewidth=1,
            dodge=False,
            size=7,
            ax=ax,
            palette=pal,
        )

        sns.boxplot(
            x="pseudogroup",
            y=g002_label,
            data=data_for_box,
            linewidth=1,
            hue="weeks",
            dodge=False,
            whis=[10, 90],
            fliersize=0,
            palette=pal,
            ax=ax,
        )
        # return actual files
        data_for_box = pre_post.copy(deep=True)
        ax.legend_.remove()
        adjust_boxplot(ax)
        if plot_log == "log":
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: format_y_axis(y, lt)))
            ax.yaxis.set_major_locator(mtick.LogLocator(base=10, numticks=10))
        elif plot_log == "symlog":
            ax.set_yscale("symlog", linthresh=10)
            ax.yaxis.set_minor_locator(MinorSymLogLocator(10))  # custom class to handle sym log ticks.
            ax.yaxis.set_major_formatter(mtick.ScalarFormatter())  # 10 and 100 are major ticks on a sym log
        ax.set_ylim(low, high)
        ax.xaxis.set_ticklabels([])
        ax.set_xlabel("")
        ax.set_ylabel(y_label)

        return data_for_box

    axes[0][0].letter = "A"
    axes[0][1].letter = "B"
    axes[1][0].letter = "C"
    axes[1][1].letter = "D"
    axes[2][0].letter = "E"
    axes[2][1].letter = "F"

    ax_A = axes[0][0]
    ax_B = axes[0][1]
    ax_C = axes[1][0]
    ax_D = axes[1][1]
    ax_E = axes[2][0]
    ax_F = axes[2][1]

    # antigen specific
    g002_label = "Percent antigen-specific among IgG^{+}"
    y_label = "% core$^{++}$ among\nIgG$^+$ B cells"
    ax = ax_A
    df = _plot_boost_frequency_figure(g002_label, y_label, ax)
    dfs.append(annotate_df(df=df, ax=ax, key=[g002_label], label=y_label))

    # Plot 2 - Epitope Specificity
    g002_label = "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}"
    y_label = "% CD4bs-specific among\nIgG$^{+}$ B cells"
    ax = ax_B
    df = _plot_boost_frequency_figure(g002_label, y_label, ax)
    dfs.append(annotate_df(df=df, ax=ax, key=[g002_label], label=y_label))

    # Plot 3 specificity among KO
    g002_label = "Percent IgG^{+}KO^- among Ag^{++}"
    y_label = "% of core$^{++}$IgG$^{+}$\nB cells that are KO$^{-}$"
    ax = ax_F
    ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
    df = _plot_boost_frequency_figure(g002_label, y_label, ax, low=0, high=101, plot_log="linear")
    dfs.append(annotate_df(df=df, ax=ax, key=[g002_label], label=y_label))

    # Plot 4, count VRC01 class among sequences
    g002_label = "Number of IGHG sequences that are VRC01-class"
    y_label = "Number of core$^{++}$KO$^{-}$\nVRC01-class IgG$^{+}$ B cells"
    ax = ax_C
    df = _plot_boost_frequency_figure(g002_label, y_label, ax, low=0, high=1000, plot_log="symlog")
    dfs.append(annotate_df(df=df, ax=ax, key=[g002_label], label=y_label))

    # plot 5, VRC01 among cells
    g002_label = "Percent of VRC01-class sequences among IgG"
    y_label = "% core$^{++}$KO$^{-}$\nVRC01-class B cells\namong IgG$^{+}$ memory B cells"
    ax = ax_D
    df = _plot_boost_frequency_figure(
        g002_label,
        y_label,
        ax,
        plot_lt=True,
        low=8e-6,
        lt=1e-5,
        high=1,
        plot_log="log",
    )
    dfs.append(annotate_df(df=df, ax=ax, key=[g002_label], label=y_label))

    ax = ax_E
    y_label = "% of participants\nwith core$^{++}$KO$^{-}$\nVRC01-class IgG$^{+}$ B cells"
    # plot 6, response
    sns.pointplot(
        data=pre_post,
        y="Response y",
        x="pseudogroup",
        hue="weeks",
        join=False,
        ci="wilson",
        ax=ax,
        dodge=0.25,
        palette=pal,
    )

    ax.set_xlabel("")
    ax.set_ylabel(
        y_label,
        labelpad=10,
    )
    ax.yaxis.set_major_locator(mtick.MultipleLocator(base=0.1))
    ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: str(int(float(y) * 100))))

    dfs.append(annotate_df(df=pre_post, ax=ax, key=["Response y"], label=y_label))

    # # set the bottom
    for ax in axes[2]:
        ax.set_xticklabels(
            [
                "Gp 4\nwk -5",
                "Gp 4\nwk 4",
                "Gp 4\nwk 8",
                "Gp 2\nwk 8",
            ],
            # rotation=45,
            # verticalalignment="top",
            # horizontalalignment="right",
        )
    custom_lines = []
    pal = {
        -5: "#9567BD",
        4: "#17BFD0",
        8: "gold",
    }

    for x in pal:
        custom_lines.append(Patch(facecolor=pal[x], edgecolor="black", linewidth=1, label=x))

    for ax in axes[2]:
        ax.legend(
            custom_lines,
            ["wk " + i._label for i in custom_lines],
            loc="upper center",
            frameon=False,
            handlelength=0.8,
            ncol=3,
            bbox_to_anchor=(0.5, -0.25),
            labelspacing=0.1,
        )

    for ax in axes.flatten():
        ax.text(-0.4, 1.1, ax.letter, transform=ax.transAxes, size=14, weight="bold")
    sns.despine()
    plt.tight_layout()
    return fig, dfs
    # figure.savefig(str(outpath) + ".png", dpi=300)


def plot_cp_frequency(data: Data, img_outpath: Path, metric_outdir: Path) -> None:
    """Sup 25"""
    letter_map = {
        0: {0: "A", 1: "B"},
        0: {0: "C", 1: "D"},
    }

    def plot(df, y, ax, name: str):
        ax = sns.stripplot(
            data=df,
            x="pseudogroup",
            hue="weeks",
            y=y,
            edgecolor="black",
            dodge=False,
            linewidth=1,
            palette=data.get_week_palette(),
            ax=ax,
            order=range(1, 8),
        )
        ax = sns.boxplot(
            data=df,
            x="pseudogroup",
            hue="weeks",
            y=y,
            whis=[10, 90],
            dodge=False,
            fliersize=0,
            palette=data.get_week_palette(),
            ax=ax,
            order=range(1, 8),
        )
        ax.legend_ = None
        ax.set_xlabel("")
        sns.despine()
        ax.set_xticklabels(
            [],
            rotation=45,
            ha="right",
            rotation_mode="anchor",
        )
        adjust_boxplot(ax)
        return (
            ax,
            df.groupby("pseudogroup")[y].median().to_frame(f"median {name}"),
            data.populate_psname(df[["pubID", y, "pseudogroup", "weeks"]]).sort_values(["pseudogroup", y]),
        )

    def dplot(df, y, ax, adjust_data_to_lt: bool = True, name: str = None):
        data_for_strip = df.copy(deep=True)
        data_for_box = df.copy(deep=True)
        lt = 1e-5
        if adjust_data_to_lt:
            print("hit", y)
            # any data less than
            di = df[df[y] < lt].index

            # any data less than, just make it that lt number
            data_for_strip.loc[di, y] = lt

            # but for box, just remove it
            data_for_box.loc[di, y] = np.nan

        ax = sns.stripplot(
            data=data_for_strip,
            x="pseudogroup",
            hue="weeks",
            y=y,
            edgecolor="black",
            dodge=False,
            linewidth=1,
            palette=data.get_week_palette(),
            ax=ax,
            order=range(1, 8),
        )
        ax = sns.boxplot(
            data=data_for_box,
            x="pseudogroup",
            hue="weeks",
            y=y,
            whis=[10, 90],
            dodge=False,
            fliersize=0,
            palette=data.get_week_palette(),
            ax=ax,
            order=range(1, 8),
        )
        # return actual values for metric files
        data_for_box = df.copy(deep=True)
        ax.legend_ = None
        ax.set_xlabel("")
        sns.despine()
        ax.set_xticklabels(
            [],
            rotation=45,
            ha="right",
            rotation_mode="anchor",
        )
        adjust_boxplot(ax)
        return (
            ax,
            data_for_box.groupby("pseudogroup")[y].median().to_frame(f"median {name}"),
            data.populate_psname(data_for_box[["pubID", y, "pseudogroup", "weeks"]]).sort_values(["pseudogroup", y]),
        )

    def pplot(df, y, ax, name):
        ax = sns.pointplot(
            data=df,
            y="Response y",
            x="pseudogroup",
            hue="weeks",
            ci="wilson",
            join=False,
            palette=data.get_week_palette(),
            ax=ax,
        )
        ax.legend_ = None
        ax.set_xlabel("")
        ax.set_xticklabels([])
        sns.despine()

        medians = (
            df.groupby("pseudogroup")["Response y"]
            .value_counts(normalize=True, dropna=False)
            .to_frame(f"median {name}")
            .reset_index()
            .set_index("pseudogroup")
            .query("`Response y`==1.0")
        )
        medians.drop(columns="Response y", inplace=True)
        return (
            ax,
            medians,
            data.populate_psname(df[["pubID", "Response y", "pseudogroup", "weeks"]]).sort_values(
                ["pseudogroup", "Response y"]
            ),
        )

    def calculate_resonse_boost(
        df: pd.DataFrame,
        option: str,
    ) -> pd.DataFrame:
        df["Response y "] = 0
        print(df.shape)
        print(
            df[df["Percent of VRC01-class sequences among IgG"].isna()][
                ["pubID", "pseudogroup", "Percent of VRC01-class sequences among IgG"]
            ]
        )
        df = df[~df["Percent of VRC01-class sequences among IgG"].isna()]
        print(df.shape)
        """calculate response rate given our critera from paper 1"""
        for index, data in df.iterrows():
            if (data["Percent of VRC01-class sequences among IgG"] > 0) & (data[option] > 0):
                df.loc[index, "Response y"] = 1  # type: ignore
            elif data["Percent of VRC01-class sequences among IgG"] == 0:
                df.loc[index, "Response y"] = 0  # type: ignore
            else:
                df.loc[index, "Response y"] = 0  # type: ignore
        return df

    def calc_clonality(df):
        df = df.query(f"is_vrc01_class==True").query("top_c_call=='IGHG'")
        clonality_df = (
            df.groupby(["pubID", "ptid", "pseudogroup", "weeks"])
            .apply(
                lambda x: pd.Series(
                    {
                        "clonality": len(x["cluster"].unique()) / len(x),
                        "num_seqs": len(x),
                        "num_clusters": len(x["cluster"].unique()),
                    }
                )
            )
            .reset_index()
        )
        return data.populate_psname(clonality_df)

    flow_seq = data.get_g002_flow_and_seq_boost_plus()
    flow_seq["weeks"] = flow_seq["weeks"].astype(int)
    flow_seq = flow_seq.query("weeks!=4")
    flow_seq = flow_seq[~flow_seq["pseudogroup"].isna()]
    seq = data.get_g002_sequences_boost_plus().query("is_vrc01_class==True").query("top_c_call=='IGHG'")
    seq["weeks"] = seq["weeks"].astype(int)
    seq = seq.query("weeks!=4")
    # breakpoint()

    fig, axes = plt.subplots(4, 2, figsize=(10, 14))

    ax_A = axes[0, 0]
    ax_B = axes[1, 0]
    ax_C = axes[2, 0]
    ax_D = axes[3, 0]
    ax_E = axes[0, 1]
    ax_F = axes[1, 1]
    ax_G = axes[2, 1]
    ax_H = axes[3, 1]

    letter = "C"
    stem = (
        f"fig{letter}"
        + " percent of IgG$^{+}$ memory B cells detected as VRC01-class with key HC residues > 4".replace(" ", "_")
    )
    ax, medians, df = dplot(
        flow_seq,
        "Percent VRC01 gt n among IgG",
        ax=ax_C,
        adjust_data_to_lt=True,
        name="% of IgG$^{+}$ memory B cells\ndetected as VRC01-class\nwith key HC residues > 4",
    )
    print(medians)
    df.to_csv(metric_outdir / f"{stem}.csv", index=False)
    ax.set_title(letter, fontsize=18, fontweight="bold", x=-0.35, y=1.05)
    ax.set_ylim(8e-6, 0.15)
    ax.set_yscale("log")
    ax.set_ylabel(
        "% of IgG$^{+}$ memory B cells\ndetected as core$^{++}$KO$^{-}$\nVRC01-class\nwith key HC residues > 4"
    )
    ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: format_y_axis(y, 1e-5)))
    ax.yaxis.set_major_locator(mtick.LogLocator(base=10, numticks=10))

    letter = "G"
    stem = (
        f"fig{letter}"
        + " percent of IgG$^{+}$ memory B cells detected as VRC01-class with key HCDR2 residues > 1".replace(" ", "_")
    )
    ax, medians2, df = dplot(
        flow_seq,
        "Percent VRC01 gt n among IgG hcdr2",
        ax=ax_G,
        adjust_data_to_lt=True,
        name="% of IgG$^{+}$ memory B cells\ndetected as VRC01-class\nwith key HCDR2 residues > 1",
    )
    df.to_csv(metric_outdir / f"{stem}.csv", index=False)
    ax.set_title(letter, fontsize=18, fontweight="bold", x=-0.35, y=1.05)
    ax.set_ylim(8e-6, 0.15)
    ax.set_yscale("log")
    ax.set_ylabel(
        "% of IgG$^{+}$ memory B cells\ndetected as core$^{++}$KO$^{-}$\nVRC01-class\nwith key HCDR2 residues > 1"
    )
    ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: format_y_axis(y, 1e-5)))
    ax.yaxis.set_major_locator(mtick.LogLocator(base=10, numticks=10))
    # _ = ax.set_yticklabels([])

    letter = "B"
    stem = f"fig{letter}" + " percent of participants with core$^{++}$KO$^{-}$ VRC01-class IgG$^{+}$ B cells with key HC residues > 4".replace(
        " ", "_"
    ).replace(
        "-", "_"
    ).replace(
        "%", "percent"
    )
    ax, medians5, df = pplot(
        calculate_resonse_boost(flow_seq.copy(deep=True), "fraction_of_cottrell_gt_n"),
        "Response y",
        ax=ax_B,
        name="% VRC01-class responders\nwith key HC residues > 4",
    )
    ax.set_title(letter, fontsize=18, fontweight="bold", x=-0.35, y=1.05)
    df.to_csv(metric_outdir / f"{stem}.csv", index=False)
    ax.set_ylabel(
        "% of participants\nwith core$^{++}$KO$^{-}$\nVRC01-class IgG$^{+}$ B cells\nwith key HC residues > 4"
    )
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))
    ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
    ax.set_ylim(0, 1.05)

    letter = "F"
    stem = f"fig{letter}" + " percent VRC01-class responders with key HCDR2 residues > 1".replace(" ", "_")
    ax, medians6, df = pplot(
        calculate_resonse_boost(flow_seq.copy(deep=True), "fraction_of_cottrell_gt_hcdr2_n"),
        "Response y",
        ax=ax_F,
        name="% VRC01-class responders\nwith key HCDR2 residues > 1",
    )
    df.to_csv(metric_outdir / f"{stem}.csv", index=False)
    ax.set_ylabel(
        "% of participants\nwith core$^{++}$KO$^{-}$\nVRC01-class IgG$^{+}$ B cells\nwith key HCDR2 residues > 1"
    )
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))
    ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
    ax.set_ylim(0, 1.05)
    ax.set_title(letter, fontsize=18, fontweight="bold", x=-0.35, y=1.05)

    letter = "D"
    clonality_metrics = calc_clonality(seq[seq["fraction_of_cottrell_gt_n"] > 0])
    stem = f"fig{letter}" + " fraction unique clones for VRC01-class IgG$^{+}$ BCRs with key HC residues > 4".replace(
        " ", "_"
    )
    ax, medians3, df = plot(
        clonality_metrics,
        "clonality",
        ax=ax_D,
        name="Fraction Unique Clones\nfor VRC01-class IgG$^{+}$ BCRs\nwith key HC residues > 4",
    )
    df.to_csv(metric_outdir / f"{stem}.csv", index=False)
    ax.set_title(letter, fontsize=18, fontweight="bold", x=-0.35, y=1.05)
    ax.set_ylabel(
        "Fraction Unique Clones\nfor core$^{++}$KO$^{-}$\nVRC01-class IgG$^{+}$ BCRs\nwith key HC residues > 4"
    )
    _ = ax.set_xticklabels(
        [
            "eOD",
            "core",
            r"eOD$\rightarrow$eOD",
            r"eOD$\rightarrow$core",
            r"eOD$\rightarrow$eOD",
            r"eOD$\rightarrow$core",
            r"eOD$\rightarrow$eOD$\rightarrow$core",
        ],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
    )

    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))
    ax.set_ylim(0, 1.05)

    letter = "H"
    stem = (
        f"fig{letter}"
        + " Fraction Unique Clones for core$^{++}$KO$^{-}$ VRC01-class IgG$^{+}$ BCRs with key HCDR2 residues > 1".replace(
            " ", "_"
        )
    )

    ax, medians4, df = plot(
        calc_clonality(seq[seq["fraction_of_cottrell_gt_hcdr2_n"] > 0]),
        "clonality",
        name="Fraction Unique Clones\nfor VRC01-class IgG$^{+}$ BCRs\nwith key HCDR2 residues > 1",
        ax=ax_H,
    )
    df.to_csv(metric_outdir / f"{stem}.csv", index=False)
    ax.set_title(letter, fontsize=18, fontweight="bold", x=-0.35, y=1.05)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))
    ax.set_ylabel(
        "Fraction Unique Clones\nfor core$^{++}$KO$^{-}$\nVRC01-class IgG$^{+}$ BCRs\nwith key HCDR2 residues > 1"
    )
    ax.set_ylim(0, 1.05)
    _ = ax.set_xticklabels(
        [
            "eOD",
            "core",
            r"eOD$\rightarrow$eOD",
            r"eOD$\rightarrow$core",
            r"eOD$\rightarrow$eOD",
            r"eOD$\rightarrow$core",
            r"eOD$\rightarrow$eOD$\rightarrow$core",
        ],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
    )

    sub_seq = (
        seq.query("is_vrc01_class")
        .query("top_c_call=='IGHG'")
        .sort_values(["run_date", "pubID"], ascending=False)
        .drop_duplicates(
            [
                # "run_date",
                "pubID",
                "ptid",
                "group",
                "weeks",
                "visit_id",
                "probe_set",
                "sample_type",
            ]
        )
    )
    letter = "A"
    stem = f"fig{letter}" + " percent of VRC01-class IgG$^{+}$ sequences with key HC residues > 4".replace(" ", "_")
    ax, _, df = plot(
        sub_seq,
        "fraction_of_cottrell_gt_n",
        ax=ax_A,
        name="percent of VRC01-class IgG$^{+}$ sequences\nwith key HC residues > 4",
    )
    # breakpoint()
    df.to_csv(metric_outdir / f"{stem}.csv", index=False)
    ax.set_title(letter, fontsize=18, fontweight="bold", x=-0.35, y=1.05)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))
    ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
    ax.set_ylabel("% of VRC01-class IgG$^{+}$ sequences\nwith key HC residues > 4")
    ax.set_ylim(0, 0.3)
    # breakpoint()

    letter = "E"
    stem = f"fig{letter}" + " percent of VRC01-class IgG$^{+}$ sequences with key HCDR2 residues > 1".replace(" ", "_")
    ax, _, df = plot(
        sub_seq,
        "fraction_of_cottrell_gt_hcdr2_n",
        ax=ax_E,
        name="% of VRC01-class IgG$^{+}$ sequences\nwith key HCDR2 residues > 1",
    )
    df.to_csv(metric_outdir / f"{stem}.csv", index=False)
    ax.set_title(letter, fontsize=18, fontweight="bold", x=-0.35, y=1.05)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))
    ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
    # ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
    ax.set_ylabel("% of VRC01-class IgG$^{+}$ sequences\nwith key HCDR2 residues > 1")
    ax.set_ylim(0, 0.3)

    # custom_lines = []
    # pal = {
    #     8: "gold",
    #     16: "#E377C2",
    #     24: "#2078B4",
    # }
    # for x in pal:
    #     custom_lines.append(
    #         Patch(facecolor=pal[x], edgecolor="black", linewidth=1, label=x)
    #     )
    # ax.legend(
    #     custom_lines,
    #     ["wk " + i._label for i in custom_lines],
    #     loc="upper center",
    #     frameon=False,
    #     handlelength=0.8,
    #     ncol=3,
    #     bbox_to_anchor=(0.5, -0.3),
    #     labelspacing=0.1,
    # )
    # DEBUG: sanity check medians
    all_medians = pd.concat([medians5, medians6, medians, medians2, medians3, medians4], axis=1).sort_index()
    legend(
        ax=axes[3, 0],
        pallete={
            "wk 8": "gold",
            "wk 16": "#E377C2",
            "wk 24": "#2078B4",
        },
        fontsize=12,
        bbox_to_anchor=(0.5, -0.65),
    )
    legend(
        ax=axes[3, 1],
        pallete={
            "wk 8": "gold",
            "wk 16": "#E377C2",
            "wk 24": "#2078B4",
        },
        fontsize=12,
        bbox_to_anchor=(0.5, -0.65),
    )
    plt.tight_layout()
    # breakpoint()
    fig.savefig(img_outpath, dpi=700)


def plot_group_4_responders(data: Data, outpath: Path) -> None:
    """Supp 20

    Parameters
    ----------
    data : Data
        Data object of source data paths and dataframes
    outpath : Path
        Output path
    """
    plt.rcParams.update({"font.size": 40})
    df = data.get_pre_post_core_df()
    count = 0
    bad_count = 0
    for _, gf in df.query("group==4").query("weeks!=-5").groupby(["ptid"]):
        if (gf["Response y"] == 1).any():
            count += 1
        else:
            bad_count += 1

    fig, axes = plt.subplots(1, 1, figsize=(5, 4))
    ax = pd.Series({"VRC01-class Responders": count, "Non-Responders": bad_count}).plot(
        kind="bar",
        stacked=True,
        linewidth=1,
        edgecolor="black",
        color=["#91FCC0", "grey"],
        ax=axes,
    )
    ax.set_ylabel("Number of Participants")
    ax.set_xticklabels(
        ["VRC01-class\nResponders", "VRC01-class\nNon-Responders"],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
    )
    sns.despine()
    plt.tight_layout()
    fig.savefig(str(outpath) + ".png", bbox_inches="tight", dpi=300)
