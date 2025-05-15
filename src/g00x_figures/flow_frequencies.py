import contextlib
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import ticker as mtick
from matplotlib.patches import Patch

from g00x_figures.data import Data, calculate_resonse
from g00x_figures.plot_helpers.boxplot import adjust_boxplot, format_y_axis
from g00x_figures.plot_helpers.font import apply_global_font_settings
from g00x_figures.plot_helpers.legend import plot_legend as legend
from g00x_figures.plot_helpers.symlog import MinorSymLogLocator
from g00x_figures.wilson import wilsonify

pal = {
    "G001": "#6C65FF",
    "G002": "#91FCC0",
    "G003": "#E377C2",
}
pallete = pal
week_ticks = ["-5", "4", "8", "10", "16", "24/21"]
week_order = [-5, 4, 8, 10, 16, 21]
yaxis_fontsize = 16
xaxis_fontsize = 16
apply_global_font_settings()
set_letters = True
plot_xticklabels = False
plot_legend = False


def normalize_threshold_to_value(
    df: pd.DataFrame,
    col: str,
    eq: str,
    threshold: float | int,
    value: float | int,
):
    norm_df = df.copy(deep=True)
    di = norm_df.query(f"`{col}` {eq} {threshold}").index
    norm_df.loc[di, col] = value
    return norm_df


def unsync_weeks(_df: pd.DataFrame) -> pd.DataFrame:
    df = _df.copy(deep=True)
    df.loc[(df["trial"] == "G002") & (df["weeks"] == 21), "weeks"] = 24
    return df


def generic_plot_box_and_whisker(
    df: pd.DataFrame,
    x: str,
    y: str,
    hue: str,
    y_label: str,
    letter: str,
    metric_outdir: Path | str,
    x_label: str = "",
    lt: float = 1e-4,
    pal: dict[str, str] = {
        "G001": "#6C65FF",
        "G002": "#91FCC0",
        "G003": "#E377C2",
    },
    y_min: float = 8e-5,
    y_max: float | int = 10,
    plot_legend: bool = False,
    adjust_data_to_lt: bool = True,
    plot_log: str = "log",
    plot_xticklabels: bool = False,
    ylabelpad: int = 0,
    strip_size: int = 5,
    ax: plt.Axes = None,
    use_raw_df: bool = False,
    use_threshold: bool = False,
) -> tuple[plt.Axes, pd.DataFrame, pd.DataFrame]:
    # combine g001 and g002 data axis
    week_ticks = ["-5", "4", "8", "10", "16", "24/21"]
    week_order = [-5, 4, 8, 10, 16, 21]

    # sync trials weeks
    df.loc[(df["trial"] == "G002") & (df["weeks"] == 24), "weeks"] = 21

    # make all data under threshold the same value
    if use_threshold:
        data_for_strip = normalize_threshold_to_value(
            df=df,
            col=y,
            eq="<",
            threshold=lt,
            value=lt,
        )

        # We do not want to include data under threshold in the box plot to misspresent the median
        data_for_box = normalize_threshold_to_value(
            df=df,
            col=y,
            eq="<",
            threshold=lt,
            value=np.nan,
        )
    else:
        data_for_strip = df.copy(deep=True)
        data_for_box = df.copy(deep=True)
    # if "Response x" is not null except else nan
    # if use_resonders:
    #     df["is_responder"] = df["Response x"].apply(
    #         lambda x: True if x > 0 else False
    #     )
    #     data_for_box = df[df["is_responder"] == True]

    yname = y_label.replace(" ", "_").replace("\n", "").replace("%", "percent").strip().lower()
    if use_raw_df:
        median_metric_df = unsync_weeks(df.copy(deep=True))[["trial", "pubID", x, y]].groupby(["trial", x])[y].median()
    else:
        median_metric_df = (
            unsync_weeks(data_for_box.copy(deep=True))[["trial", "pubID", x, y]].groupby(["trial", x])[y].median()
        )
    # breakpoint()
    median_metric_df.name = letter + " " + median_metric_df.name

    # print(median_metric_df)
    # median_metric_df.to_csv(
    #     Path(metric_outdir) / f"fig{letter}_median_{yname}.csv",
    #     index=True,
    # )

    pal = {k: v for k, v in pal.items() if k in data_for_strip[hue].unique()}

    ax = sns.stripplot(
        x=x,
        y=y,
        hue=hue,
        edgecolor="black",
        linewidth=1,
        s=strip_size,
        jitter=0.1,  # type: ignore
        palette=pal,
        dodge=True,
        hue_order=pal.keys(),
        data=data_for_strip,
        order=week_order,
        ax=ax if ax else None,
    )

    ax = sns.boxplot(
        x=x,
        y=y,
        hue=hue,
        dodge=True,
        palette=pal,
        fliersize=0,
        hue_order=pal.keys(),
        order=week_order,
        data=data_for_box,
        whis=[10, 90],  # type: ignore
        ax=ax,
    )
    data_for_box = df.copy(deep=True)
    metric_df = unsync_weeks(data_for_box.copy(deep=True))[["trial", "pubID", x, y]]
    metric_df.to_csv(Path(metric_outdir) / f"fig{letter}_{yname}.csv", index=False)
    ax.legend_.set_visible(False)  # type: ignore
    ax.set_ylim(y_min, y_max)
    if plot_log == "log":
        ax.set_yscale("log")

    elif plot_log == "symlog":
        ax.set_yscale("symlog", linthresh=10)

    else:
        ax.set_yscale("linear")

    adjust_boxplot(ax)
    ax.set_ylabel(y_label, labelpad=ylabelpad, fontsize=yaxis_fontsize)
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
            ncol=5,
            bbox_to_anchor=(0.5, -0.2),
            labelspacing=0.1,
            fontsize=xaxis_fontsize,
        )
    if x_label:
        ax.set_xlabel(x_label, fontsize=yaxis_fontsize)
    else:
        ax.set_xlabel("")
    if plot_xticklabels:
        ax.set_xticklabels(week_ticks, fontsize=xaxis_fontsize)
        ax.set_xlabel("Weeks post vaccination", fontsize=yaxis_fontsize)
    else:
        ax.set_xticklabels([])
        ax.set_xlabel("")

    ax.spines[["right", "top"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=xaxis_fontsize)
    ax.tick_params(axis="y", labelsize=xaxis_fontsize)
    adjust_boxplot(ax)

    # print(f"fig{letter}_{yname}.csv")

    # median_metric_df = median_metric_df.set_index(["trial", x])

    return ax, metric_df, median_metric_df


def plot_flow_frequencies(data: Data, img_outdir: Path, metric_outdir: Path, use_alt: bool = False) -> None:
    g00x_flow_and_seq = data.get_g00x_flow_and_seq()

    @contextlib.contextmanager
    def figure_context(*args, **kwargs):
        fig, axes = plt.subplots(*args, **kwargs)
        yield fig, axes
        plt.close("all")

    with figure_context(4, 2, figsize=(14, 18)) as (figure, axes):
        plt.subplots_adjust(wspace=0.3, hspace=0.3)

        plottable_df = g00x_flow_and_seq.copy(deep=True)

        median_metric_dfs = []

        y = "Percent antigen-specific among IgG^{+}"
        ylabel = "% GT8$^{++}$ among\nIgG$^+$ B cells"
        ax, metric_df, median_metric_df = generic_plot_box_and_whisker(
            df=plottable_df,
            x="weeks",
            y=y,
            hue="pseudogroup",
            letter="A",
            y_label=ylabel,
            ax=axes[0][0],
            plot_xticklabels=plot_xticklabels,
            plot_legend=plot_legend,
            metric_outdir=metric_outdir,
        )
        median_metric_dfs.append(median_metric_df)
        # TODO: replace y with ylabel and pull all y into single CSV except for last fig

        y = "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}"
        ylabel = "% CD4bs-specific among\nIgG$^{+}$ B cells"
        ax, metric_df, median_metric_df = generic_plot_box_and_whisker(
            df=plottable_df,
            x="weeks",
            y=y,
            hue="pseudogroup",
            y_label=ylabel,
            letter="B",
            ax=axes[0][1],
            plot_xticklabels=plot_xticklabels,
            plot_legend=plot_legend,
            metric_outdir=metric_outdir,
        )
        median_metric_dfs.append(median_metric_df)

        y = "Percent IgG^{+}KO^- among Ag^{++}"
        ylabel = "% of GT8$^{++}$IgG$^{+}$\nB cells that are KO$^{-}$"
        ax, metric_df, median_metric_df = generic_plot_box_and_whisker(
            df=plottable_df,
            x="weeks",
            y=y,
            hue="pseudogroup",
            letter="C",
            y_label=ylabel,
            ax=axes[1][0],
            y_min=0,
            y_max=105,
            plot_log="linear",
            ylabelpad=20,
            plot_xticklabels=plot_xticklabels,
            plot_legend=plot_legend,
            metric_outdir=metric_outdir,
        )
        median_metric_dfs.append(median_metric_df)

        y = "Number of IGHG sequences that are VRC01-class"
        ylabel = "Number of VRC01-class\nIgG$^{+}$ B cells"
        ax, metric_df, median_metric_df = generic_plot_box_and_whisker(
            df=plottable_df,
            x="weeks",
            y=y,
            hue="pseudogroup",
            letter="D",
            y_label=ylabel,
            ax=axes[1][1],
            plot_log="symlog",
            y_min=0,
            y_max=1800,
            adjust_data_to_lt=False,
            ylabelpad=13,
            plot_xticklabels=plot_xticklabels,
            plot_legend=plot_legend,
            metric_outdir=metric_outdir,
        )
        median_metric_dfs.append(median_metric_df)

        # TODO: A version of Fig. 2E that shows
        # “% of all B cells detected as VRC01-class” (rather than % of IgG memory B cells).
        # Original value in 2E divided by the fraction-of-all-B-cells-that-are-IgG
        if use_alt:
            y = "Percent of VRC01-class sequences among IgD-"
            ylabel = "% of IgD^- B cells\ndetected as VRC01-class"
        else:
            y = "Percent of VRC01-class sequences among IgG"
            ylabel = "% of IgG$^{+}$ memory B cells\ndetected as VRC01-class"

        ax, metric_df, median_metric_df = generic_plot_box_and_whisker(
            ax=axes[2][0],
            df=plottable_df,
            x="weeks",
            y=y,
            hue="pseudogroup",
            letter="E",
            y_label=ylabel,
            y_max=1,
            plot_xticklabels=plot_xticklabels,
            plot_legend=plot_legend,
            lt=1e-5,
            use_threshold=True,
            y_min=8e-6,  # TODO: finish this so normal wont be changed
            metric_outdir=metric_outdir,
            use_raw_df=False,
        )
        median_metric_dfs.append(median_metric_df)

        y_label = "%VRC01-class\nrepsonders"
        yname = "Percent of VRC01-class responders".replace(" ", "_").lower()
        _df = plottable_df = plottable_df[~plottable_df["Percent of VRC01-class sequences among IgG"].isna()]
        _df[yname] = _df["Response x"]
        # TODO:remove ci for -5
        # _df.loc[_df[_df["weeks"] == -5].index, yname] = 10
        ax = sns.pointplot(
            data=_df,
            hue="pseudogroup",
            y=yname,
            x="weeks",
            join=False,
            ci="wilson",
            hue_order=pal.keys(),
            # order=["-5", "4", "8", "10", "16", "21"],
            ax=axes[3][1],
            plot_legend=plot_legend,
            palette=pal,
            dodge=0.25,
        )
        # wilsonify(
        #     g=ax,
        #     data=_df.query("weeks != -5"),
        #     x="weeks",
        #     y=yname,
        #     hue="pseudogroup",
        # )
        # median to csv
        median_metric_df = (
            unsync_weeks(_df)
            .groupby(["trial", "weeks"])[yname]
            .median()
            # .reset_index()
        )
        # breakpoint()
        # median_metric_df = median_metric_df.rename({"Response x": yname})
        # removing pointplot medians
        # median_metric_dfs.append(median_metric_df)
        # median_metric_df.to_csv(
        #     metric_outdir / f"figH_median_{yname}.csv", index=False
        # )
        metric_df = _df[["pubID", "trial", yname, "weeks"]]
        metric_df.to_csv(metric_outdir / f"figH_{yname}.csv", index=False)
        axes[3][1].set_xlabel("Weeks post vaccination", fontsize=xaxis_fontsize)
        axes[3][1].set_xticklabels(
            labels=["-5", "4", "8", "10", "16", "24/21"],
            fontsize=xaxis_fontsize,
        )
        axes[3][1].legend_.remove()
        # axes[3][1].set_xticklabels([])
        axes[3][1].set_ylabel(y_label, labelpad=20, fontsize=yaxis_fontsize)
        axes[3][1].yaxis.set_major_locator(mtick.MultipleLocator(base=0.2))
        axes[3][1].yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: str(int(float(y) * 100))))
        legend(
            ax=axes[3][1],
            loc="upper center",
            pallete=pal,
            bbox_to_anchor=(0.5, -0.2),
            handlelength=0.8,
            labelspacing=0.1,
        )
        axes[3][1].tick_params(axis="x", labelsize=xaxis_fontsize)
        axes[3][1].tick_params(axis="y", labelsize=xaxis_fontsize)

        if use_alt:
            a = "Percent of IGHD sequences that are VRC01-class"
            b = "Percent of epitope-specific (KO^-Ag^{++}) among IgD^-"
            c = "Percent antigen-specific among IgD^-"
            d = "Percent VRC01-class among eOD-specific IgD- memory BCR sequences"
            y_label = "% VRC01-class among\nCD4bs-specific\nIgD$^{-}$ memory sequences"
        else:
            a = "Percent of IGHG sequences that are VRC01-class"
            b = "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}"
            c = "Percent antigen-specific among IgG^{+}"
            d = "Percent VRC01-class among eOD-specific IgG+ memory BCR sequences"
            y_label = "% VRC01-class among\nGT8-specific\nIgG$^{+}$ memory sequences"
        # this is how we derive D.
        plottable_df[d] = plottable_df[a] * (plottable_df[b] / plottable_df[c])

        # g002_label = a
        # g001_label = "Percent of epitope-specific (KO-GT8++) sequenced IgG BCRs that are VRC01-class"
        ax, metric_df, median_metric_df = generic_plot_box_and_whisker(
            df=plottable_df,
            x="weeks",
            y=a,
            hue="pseudogroup",
            letter="F",
            y_label="% VRC01-class among\nCD4bs-specific\nIgG$^{+}$ memory sequences",
            plot_log="linear",
            y_max=104,
            y_min=0,
            plot_xticklabels=plot_xticklabels,
            plot_legend=plot_legend,
            ylabelpad=10,
            ax=axes[2][1],
            metric_outdir=metric_outdir,
        )
        median_metric_dfs.append(median_metric_df)

        y_label = "% VRC01-class among\nGT8-specific\nIgG$^{+}$ memory sequences"
        ax, metric_df, median_metric_df = generic_plot_box_and_whisker(
            df=plottable_df,
            x="weeks",
            y=d,
            hue="pseudogroup",
            letter="G",
            y_label=y_label,
            x_label="Weeks post vaccination",
            plot_log="linear",
            y_max=104,
            y_min=0,
            plot_xticklabels=True,
            plot_legend=False,
            ax=axes[3][0],
            metric_outdir=metric_outdir,
        )
        legend(
            ax=axes[3][0],
            loc="upper center",
            pallete=pal,
            bbox_to_anchor=(0.5, -0.2),
            handlelength=0.8,
            labelspacing=0.1,
        )
        median_metric_dfs.append(median_metric_df)

        if set_letters:
            for ax, l in zip(axes.flat, ["A", "B", "C", "D", "E", "F", "G", "H"]):
                ax.set_title(l, fontsize=20, fontweight="bold", x=-0.24, y=1.1)

        median_metric_df = pd.concat(median_metric_dfs, axis=1)
        median_metric_df.to_csv(metric_outdir / f"fig2_median_all.csv", index=True)
        figure.savefig(
            img_outdir / f"percent_x_among_igg_bcells.png",
            bbox_inches="tight",
            dpi=500,
        )
        # plt.show()
