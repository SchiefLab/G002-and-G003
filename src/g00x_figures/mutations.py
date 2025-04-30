## Plotting somatic mutation frequencies

import contextlib
from pathlib import Path
from typing import Any, List, Literal, Tuple, Union

import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
import patchworklib as pw
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

from g00x_figures.box_and_scatter.flow_frequencies import adjust_boxplot
from g00x_figures.data import Data
from g00x_figures.plot_helpers.font import apply_global_font_settings
from g00x_figures.plot_helpers.legend import plot_legend

apply_global_font_settings()
# TODO: verify if these can be properly added to the global font settings


def unsync_weeks(_df: pd.DataFrame) -> pd.DataFrame:
    df = _df.copy(deep=True)
    df.loc[(df["trial"] == "G002") & (df["weeks"] == 21), "weeks"] = 24
    return df


# TODO: check if live
def plot_somatic_mutation_nt_frequencies(data: Data, outpath: Path) -> None:
    # g002 sequences
    g002_sequences = data.get_g002_sequences_prime()
    g001_sequences = data.get_g001_sequences_prime()
    trial_palette = data.get_trial_palette()
    figure, axes = plt.subplots(2, 2, figsize=(8.518, 6.5))
    for ax, metric in zip(axes[0], ["v_mutation_heavy", "v_mutation_light"]):
        groupby_g001 = (
            g001_sequences.query("weeks!=-5")
            .query("is_vrc01_class==True")
            .groupby(["pubID", "weeks"])[metric]
            .median()
            .reset_index()
            .assign(trial="G001")
        )
        groupby_g002 = (
            g002_sequences.query("weeks!=-5")
            .query("top_c_call=='IGHG'")
            .query("is_vrc01_class==True")
            .groupby(["pubID", "weeks"])[metric]
            .median()
            .reset_index()
            .assign(trial="G002")
        )
        trial_concat = pd.concat([groupby_g001, groupby_g002], axis=0).reset_index(drop=True)
        trial_concat["weeks"] = trial_concat["weeks"].astype(int)

        # Plot 1 - Heavy chain V gene
        sns.stripplot(
            x="weeks",
            y=metric,
            hue="trial",
            edgecolor="black",
            linewidth=1,
            s=7,
            jitter=0.15,
            dodge=True,
            palette=trial_palette,
            hue_order=["G001", "G002"],
            order=[4, 8, 16, 24],
            data=trial_concat,
            ax=ax,
        )
        sns.boxplot(
            x="weeks",
            y=metric,
            hue="trial",
            dodge=True,
            fliersize=0,
            palette=trial_palette,
            order=[4, 8, 16, 24],
            hue_order=["G001", "G002"],
            data=trial_concat,
            ax=ax,
        )
        adjust_boxplot(ax)
        ax.legend_.set_visible(False)
        if "heavy" in metric:
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))
        else:
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=1))
        if metric == "v_mutation_heavy":
            label = r"$\mathregular{V_H}$" + f" gene\n%mutation (nt)"
        else:
            label = r"$\mathregular{V_{K/L}}$" + f" gene\n%mutation (nt)"
        ax.set_ylabel(label)
        ax.set_xticklabels([])
        ax.set_xlabel("")
    for ax, metric in zip(axes[1], ["v_mutation_heavy", "v_mutation_light"]):
        s1 = g001_sequences.query("weeks!=-5").query("is_vrc01_class==True")
        s2 = (
            g002_sequences.query("weeks!=-5")
            .query("is_vrc01_class==True")
            .assign(trial="G002")
            .query("top_c_call=='IGHG'")
        )
        concat_for_violin = pd.concat([s1, s2], axis=0).reset_index(drop=True)
        concat_for_violin["weeks"] = concat_for_violin["weeks"].astype(int)
        sns.violinplot(
            ax=ax,
            data=concat_for_violin,
            bw=0.3,
            x="weeks",
            y=metric,
            hue="trial",
            palette=trial_palette,
            inner="quartile",
            saturation=1,
            scale="area",
            order=[4, 8, 16, 24],
            cut=0,
        )
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))
        ax.legend_.set_visible(False)
        if metric == "v_mutation_heavy":
            label = r"$\mathregular{V_H}$" + f" gene\n%mutation (nt)"
        else:
            label = r"$\mathregular{V_{K/L}}$" + f" gene\n%mutation (nt)"
        ax.set_ylabel(label)
        ax.set_xlabel("Weeks post vaccination")
    sns.despine()
    custom_lines = []
    for x in trial_palette:
        custom_lines.append(
            Patch(
                facecolor=trial_palette[x],
                edgecolor="black",
                linewidth=1,
                label=x,
            )
        )
    for ax in axes[1]:
        ax.legend(
            custom_lines,
            [i._label for i in custom_lines],
            loc="upper center",
            frameon=False,
            handlelength=0.8,
            ncol=2,
            bbox_to_anchor=(0.5, -0.15),
            labelspacing=0.1,
        )

    figure.tight_layout()
    figure.savefig(str(outpath) + ".png", dpi=300)


# TODO: check if live
def plot_somatic_mutation_aa_frequencies(data: Data, outpath: Path) -> None:
    # g002 sequences
    g002_sequences = data.get_g002_sequences_prime()
    g001_sequences = data.get_g001_sequences_prime()
    trial_palette = data.get_trial_palette()
    figure, axes = plt.subplots(2, 2, figsize=(8.518, 6.5))
    for ax, metric in zip(axes[0], ["v_mutation_aa_heavy", "v_mutation_aa_light"]):
        groupby_g001 = (
            g001_sequences.query("weeks!=-5")
            .query("is_vrc01_class==True")
            .groupby(["pubID", "weeks"])[metric]
            .median()
            .reset_index()
            .assign(trial="G001")
        )
        groupby_g002 = (
            g002_sequences.query("weeks!=-5")
            .query("top_c_call=='IGHG'")
            .query("is_vrc01_class==True")
            .groupby(["pubID", "weeks"])[metric]
            .median()
            .reset_index()
            .assign(trial="G002")
        )
        trial_concat = pd.concat([groupby_g001, groupby_g002], axis=0).reset_index(drop=True)
        trial_concat["weeks"] = trial_concat["weeks"].astype(int)

        # Plot 1 - Heavy chain V gene
        sns.stripplot(
            x="weeks",
            y=metric,
            hue="trial",
            edgecolor="black",
            linewidth=1,
            s=7,
            jitter=0.15,
            dodge=True,
            palette=trial_palette,
            hue_order=["G001", "G002"],
            order=[4, 8, 16, 24],
            data=trial_concat,
            ax=ax,
        )
        sns.boxplot(
            x="weeks",
            y=metric,
            hue="trial",
            dodge=True,
            fliersize=0,
            palette=trial_palette,
            order=[4, 8, 16, 24],
            hue_order=["G001", "G002"],
            data=trial_concat,
            ax=ax,
        )
        adjust_boxplot(ax)
        ax.legend_.set_visible(False)
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))
        if metric == "v_mutation_aa_heavy":
            label = r"$\mathregular{V_H}$" + f" gene\n%mutation (aa)"
        else:
            label = r"$\mathregular{V_{K/L}}$" + f" gene\n%mutation (aa)"
        ax.set_ylabel(label)
        ax.set_xticklabels([])
        ax.set_xlabel("")
        # ax.set_xlabel("Weeks post vaccination")
    for ax, metric in zip(axes[1], ["v_mutation_aa_heavy", "v_mutation_aa_light"]):
        s1 = g001_sequences.query("weeks!=-5").query("is_vrc01_class==True")
        s2 = (
            g002_sequences.query("weeks!=-5")
            .query("is_vrc01_class==True")
            .assign(trial="G002")
            .query("top_c_call=='IGHG'")
        )
        concat_for_violin = pd.concat([s1, s2], axis=0).reset_index(drop=True)
        concat_for_violin["weeks"] = concat_for_violin["weeks"].astype(int)
        sns.violinplot(
            ax=ax,
            data=concat_for_violin,
            bw=0.3,
            x="weeks",
            y=metric,
            hue="trial",
            palette=trial_palette,
            inner="quartile",
            saturation=1,
            scale="area",
            order=[4, 8, 16, 24],
            cut=0,
        )
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))
        ax.legend_.set_visible(False)
        if metric == "v_mutation_aa_heavy":
            label = r"$\mathregular{V_H}$" + f" gene\n%mutation (aa)"
        else:
            label = r"$\mathregular{V_{K/L}}$" + f" gene\n%mutation (aa)"
        ax.set_ylabel(label)
        ax.set_xlabel("Weeks post vaccination")
    sns.despine()
    custom_lines = []
    for x in trial_palette:
        custom_lines.append(
            Patch(
                facecolor=trial_palette[x],
                edgecolor="black",
                linewidth=1,
                label=x,
            )
        )
    # for ax in axes[0]:
    #     ax.legend(
    #         custom_lines,
    #         [i._label for i in custom_lines],
    #         loc="upper center",
    #         frameon=False,
    #         handlelength=0.8,
    #         ncol=2,
    #         bbox_to_anchor=(0.5, -0.15),
    #         labelspacing=0.1,
    #     )
    for ax in axes[1]:
        ax.legend(
            custom_lines,
            [i._label for i in custom_lines],
            loc="upper center",
            frameon=False,
            handlelength=0.8,
            ncol=2,
            bbox_to_anchor=(0.5, -0.15),
            labelspacing=0.1,
        )

    figure.tight_layout()
    figure.savefig(str(outpath) + ".png", dpi=300)


# TODO: check if live
def plot_somatic_mutation_aa_frequencies_boost(data: Data, outpath: Path) -> None:
    g002_sequences = data.get_g002_sequences_boost()
    weeks_palette = data.get_week_palette()
    figure, axes = plt.subplots(2, 2, figsize=(8.518, 6.5))
    for ax, metric in zip(axes[0], ["v_mutation_aa_heavy", "v_mutation_aa_light"]):
        groupby_g002 = (
            g002_sequences.query("weeks!=-5")
            .query("top_c_call=='IGHG'")
            .query("is_vrc01_class==True")
            .query("pseudogroup!=2")
            .groupby(["pubID", "ptid", "pseudogroup", "weeks"])[metric]
            .median()
            .reset_index()
        )

        # Plot 1 - Heavy chain V gene
        sns.stripplot(
            x="pseudogroup",
            y=metric,
            hue="weeks",
            edgecolor="black",
            linewidth=1,
            s=7,
            jitter=0.15,
            dodge=False,
            palette=weeks_palette,
            data=groupby_g002,
            ax=ax,
        )
        sns.boxplot(
            x="pseudogroup",
            y=metric,
            hue="weeks",
            dodge=False,
            fliersize=0,
            palette=weeks_palette,
            data=groupby_g002,
            ax=ax,
        )

        # Sample for debuggin
        groupby_g002.rename(columns={metric: "median_" + metric}).query("weeks in [16, 24]").query(
            "pseudogroup in [3, 4, 5, 6, 7]"
        ).to_csv(str(outpath) + f"-{metric}.csv")

        adjust_boxplot(ax)
        ax.legend_.set_visible(False)
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))
        if metric == "v_mutation_aa_heavy":
            label = r"$\mathregular{V_H}$" + f" gene\n%mutation (aa)"
        else:
            label = r"$\mathregular{V_{K/L}}$" + f" gene\n%mutation (aa)"
        ax.set_ylabel(label)
        ax.set_xticklabels([])
        ax.set_xlabel("")

    for ax, metric in zip(axes[1], ["v_mutation_aa_heavy", "v_mutation_aa_light"]):
        s2 = (
            g002_sequences.query("weeks!=-5")
            .query("pseudogroup!=2")
            .query("is_vrc01_class==True")
            .query("top_c_call=='IGHG'")
        )
        sns.violinplot(
            ax=ax,
            data=s2,
            bw=0.3,
            x="pseudogroup",
            y=metric,
            hue="weeks",
            palette=weeks_palette,
            inner="quartile",
            dodge=False,
            saturation=1,
            scale="area",
            cut=0,
        )
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))
        ax.legend_.set_visible(False)
        if metric == "v_mutation_aa_heavy":
            label = r"$\mathregular{V_H}$" + f" gene\n%mutation (aa)"
        else:
            label = r"$\mathregular{V_{K/L}}$" + f" gene\n%mutation (aa)"
        ax.set_ylabel(label)
        ax.set_xlabel("")
    for ax in axes[1]:
        ax.set_xticklabels(
            [
                r"eOD",
                r"eOD$\rightarrow$eOD",
                r"eOD$\rightarrow$core",
                r"eOD$\rightarrow$eOD",
                r"eOD$\rightarrow$core",
                r"eOD$\rightarrow$eOD$\rightarrow$core",
            ],
            rotation=45,
            verticalalignment="top",
            horizontalalignment="right",
        )
    pal = {
        8: "gold",
        16: "#E377C2",
        24: "#2078B4",
    }
    custom_lines = []
    for x in pal:
        custom_lines.append(Patch(facecolor=pal[x], edgecolor="black", linewidth=1, label=x))

    for ax in axes[1]:
        ax.legend(
            custom_lines,
            ["wk " + i._label for i in custom_lines],
            loc="upper center",
            frameon=False,
            handlelength=0.8,
            ncol=3,
            bbox_to_anchor=(0.5, -0.45),
            labelspacing=0.1,
        )
    sns.despine()
    figure.tight_layout()
    figure.savefig(str(outpath) + ".png", dpi=300)


def get_group_median(seq: pd.DataFrame, metric: str) -> pd.DataFrame:
    seq_group_median = seq.groupby(["pubID", "trial", "pseudogroup", "weeks"])[metric].median().reset_index()
    return seq_group_median


def plot_stripbox(
    ax: pw.Brick,
    df: pd.DataFrame,
    x: Literal["weeks", "pseudogroup"],
    y: Literal[
        "v_mutation_aa_heavy",
        "v_mutation_aa_light",
        "v_mutation_heavy",
        "v_mutation_light",
    ],
    hue: Literal["trial", "weeks"],
    pallete: dict,
) -> None:
    """Plot stripbox plot.

    >>> plot_stripbox(ax=ax df=g00x_seq, x=weeks, y=v_mutation_heavy, hue=trial)
    >>> plot_stripbox(ax=ax df=g002_seq, x=pseudogroup, y=v_mutation_heavy, hue=weeks)

    Parameters
    ----------
    ax : pw.Brick
        Patchworklib Brick
    df : pd.DataFrame
        Dataframe
    x : Literal["weeks", "pseudogroup"]
        x-axis variable name
    y : Literal["v_mutation_aa_heavy", "v_mutation_aa_light", "v_mutation_heavy", "v_mutation_light"]
        y-axis variable name
    hue : Literal["trial", "weeks"]
        hue variable name
    pallete : dict
        color pallete
    """

    # If there are multiple hue values for each x value, set dodge to True
    dodge = df.groupby([x])[hue].value_counts().unstack().notna().sum(axis=1).gt(1).any()
    sns.stripplot(
        x=x,
        y=y,
        hue=hue,
        edgecolor="black",
        linewidth=1,
        s=5,
        alpha=0.8,
        jitter=0.1,
        dodge=dodge,
        palette=pallete,
        hue_order=pallete.keys(),
        data=df,
        ax=ax,  # type: ignore
    )
    sns.boxplot(
        x=x,
        y=y,
        hue=hue,
        dodge=dodge,
        fliersize=0,
        palette=pallete,
        hue_order=pallete.keys(),
        data=df,
        ax=ax,  # type: ignore
    )
    adjust_boxplot(ax, median_linewidth=2.5)  # type: ignore
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0, symbol="%"))


def plot_violin(
    ax: pw.Brick,
    df: pd.DataFrame,
    x: Literal["weeks", "pseudogroup"],
    y: Literal[
        "v_mutation_aa_heavy",
        "v_mutation_aa_light",
        "v_mutation_heavy",
        "v_mutation_light",
    ],
    hue: Literal["trial", "weeks"],
    pallete: dict,
) -> None:
    """Plot violin plot.

    >>> plot_violin(ax=ax df=g00x_seq, x=weeks, y=v_mutation_heavy, hue=trial)
    >>> plot_violin(ax=ax df=g002_seq, x=pseudogroup, y=v_mutation_heavy, hue=weeks)

    Parameters
    ----------
    ax : pw.Brick
        Patchworklib Brick
    df : pd.DataFrame
        Dataframe
    x : Literal["weeks", "pseudogroup"]
        x-axis variable name
    y : Literal["v_mutation_aa_heavy", "v_mutation_aa_light", "v_mutation_heavy", "v_mutation_light"]
        y-axis variable name
    hue : Literal["trial", "weeks"]
        hue variable name
    pallete : dict
        color pallete
    """
    dodge = df.groupby([x])[hue].value_counts().unstack().notna().sum(axis=1).gt(1).any()
    sns.violinplot(
        ax=ax,
        data=df,
        bw=0.3,
        x=x,
        y=y,
        hue=hue,
        palette=pallete,
        dodge=dodge,
        inner="quartile",
        saturation=1,
        scale="area",
        hue_order=pallete.keys(),
        cut=0,
    )
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))


def plot_v_mutation(
    ax: pw.Brick,
    data: Any,
    x: Literal["weeks", "pseudogroup"],
    y: Literal[
        "v_mutation_aa_heavy",
        "v_mutation_aa_light",
        "v_mutation_heavy",
        "v_mutation_light",
    ],
    hue: Literal["trial", "weeks"],
    xlabel: str | None = None,
    ylabel: str | None = None,
    data_src: Literal["prime", "boost"] = "prime",
    plot_type: Literal["stripbox", "violin"] = "stripbox",
    rm_week_neg5: bool = True,
    rm_psuedogroup_2: bool = True,
    is_vrc01_class: bool = True,
    is_ighg: bool = True,
    use_median: bool = True,
    use_xticks: bool = True,
    week_xticks: Tuple[str, ...] = ("4", "8", "10", "16", "24/21"),
    psuedogroup_xticks: Tuple[str, ...] = (
        r"eOD",
        r"eOD$\rightarrow$eOD",
        r"eOD$\rightarrow$core",
        r"eOD$\rightarrow$eOD",
        r"eOD$\rightarrow$core",
        r"eOD$\rightarrow$eOD$\rightarrow$core",
    ),
    use_legend: bool = False,
    legend_pallete: dict | None = None,
    letter: str | None = None,
) -> pd.DataFrame:
    """Plot stripbox or violin plot for mutational group."""
    if data_src == "prime":
        df = data.get_g00x_sequences_prime().copy(deep=True)
        # sync trials weeks
        df.loc[(df["trial"] == "G002") & (df["weeks"] == 24), "weeks"] = 21
        pallete = data.get_trial_g001_g002_g003_palette()
        xticks = week_xticks
        xtick_rotation = 0
    elif data_src == "boost":
        df = data.get_g002_sequences_boost().copy(deep=True)
        pallete = data.get_week_palette()
        xticks = psuedogroup_xticks
        xtick_rotation = 45

    if rm_week_neg5:
        df = df.query("weeks != -5")
    if rm_psuedogroup_2:
        df = df.query("pseudogroup != 2")
    if is_vrc01_class:
        df = df.query("is_vrc01_class==True")
    if is_ighg:
        df = df.query("top_c_call=='IGHG'")

    # Use the median value each group instead of all points
    if use_median:
        df = get_group_median(df, y)

    if plot_type == "stripbox":
        plot_stripbox(ax=ax, df=df, x=x, y=y, hue=hue, pallete=pallete)
    elif plot_type == "violin":
        plot_violin(ax=ax, df=df, x=x, y=y, hue=hue, pallete=pallete)

    ax.legend_.set_visible(False)  # type: ignore
    ax.set_xlabel(xlabel)  # type: ignore
    ax.set_ylabel(ylabel)  # type: ignore

    ax.tick_params(axis="x")
    ax.tick_params(axis="y")

    ax.set_xticklabels([])

    sns.despine(ax=ax)

    # if use_legend:
    #     if legend_pallete:
    #         pallete = legend_pallete
    #     plot_legend(pallete=pallete, ax=ax, loc="upper center")

    if use_xticks:
        ax.set_xticklabels(
            xticks,
            ha="right",
            # verticalalignment="top",
            # horizontalalignment="center",
            rotation=xtick_rotation,
        )

    if use_legend:
        if legend_pallete and data_src == "boost":
            plot_legend(
                pallete=legend_pallete,
                ax=ax,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.45),
                fontsize=16,
            )
        elif legend_pallete:
            plot_legend(
                pallete=legend_pallete,
                ax=ax,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.2),
                fontsize=15,
            )
        else:
            plot_legend(pallete=pallete, ax=ax, loc="upper center")

    if letter:
        ax.set_title(letter, fontsize=20, fontweight="bold", x=-0.2, y=1)

    # unsync trials weeks
    df.loc[(df["trial"] == "G002") & (df["weeks"] == 21), "weeks"] = 24

    return df[["pubID", "trial", "pseudogroup", "weeks"] + [y]]


def plot_v_mutations(
    data,
    src_type: Literal["prime", "boost"],
    seq_type: Literal["aa", "nt"],
    mectric_outdir: Path,
):
    apply_global_font_settings()
    figsize = (4, 3)
    xlabel = "Weeks post vaccination"

    if src_type == "boost":
        x = "pseudogroup"
        hue = "weeks"
        legend_pallete = {
            "wk 8": "gold",
            "wk 16": "#E377C2",
            "wk 24": "#2078B4",
        }
        xlabel = None
    elif src_type == "prime":
        x = "weeks"
        hue = "trial"
        legend_pallete = data.get_trial_g001_g002_g003_palette()
        xlabel = "Weeks post vaccination"

    if seq_type == "aa":
        heavy = "aa_heavy"
        light = "aa_light"
    elif seq_type == "nt":
        heavy = "heavy"
        light = "light"

    stripbox_v_mut_heavy_g = pw.Brick(figsize=figsize)
    stripbox_v_mut_heavy_df = plot_v_mutation(
        ax=stripbox_v_mut_heavy_g,
        data=data,
        x=x,
        y=f"v_mutation_{heavy}",
        hue=hue,
        ylabel=r"$\mathregular{V_H}$" + " gene\n% mutation " + f"({seq_type})",
        data_src=src_type,
        plot_type="stripbox",
        use_legend=False,
        use_xticks=False,
        legend_pallete=legend_pallete,
        letter="A",
    )

    violin_v_mut_heavy_g = pw.Brick(figsize=figsize)
    violin_v_mut_heavy_df = plot_v_mutation(
        ax=violin_v_mut_heavy_g,
        data=data,
        x=x,
        y=f"v_mutation_{heavy}",
        hue=hue,
        ylabel=r"$\mathregular{V_H}$" + " gene\n% mutation " + f"({seq_type})",
        xlabel=xlabel,
        data_src=src_type,
        plot_type="violin",
        use_median=False,
        use_legend=True,
        use_xticks=True,
        legend_pallete=legend_pallete,
        letter="C",
    )

    stripbox_v_mut_light_g = pw.Brick(figsize=figsize)
    stripbox_v_mut_light_df = plot_v_mutation(
        ax=stripbox_v_mut_light_g,
        data=data,
        x=x,
        y=f"v_mutation_{light}",
        hue=hue,
        ylabel=r"$\mathregular{V_{K/L}}$" + " gene\n% mutation " + f"({seq_type})",
        data_src=src_type,
        plot_type="stripbox",
        use_legend=False,
        use_xticks=False,
        legend_pallete=legend_pallete,
        letter="B",
    )

    violin_v_mut_light_g = pw.Brick(figsize=figsize)
    violin_v_mut_light_df = plot_v_mutation(
        ax=violin_v_mut_light_g,
        data=data,
        x=x,
        y=f"v_mutation_{light}",
        hue=hue,
        ylabel=r"$\mathregular{V_{K/L}}$" + " gene\n% mutation " + f"({seq_type})",
        xlabel=xlabel,
        data_src=src_type,
        plot_type="violin",
        use_median=False,
        use_legend=True,
        use_xticks=True,
        legend_pallete=legend_pallete,
        letter="D",
    )

    # custom y limits
    if src_type == "prime":
        if seq_type == "aa":
            stripbox_v_mut_heavy_g.set_ylim(0, 0.10)
            stripbox_v_mut_light_g.set_ylim(0, 0.06)

        if seq_type == "nt":
            stripbox_v_mut_heavy_g.set_ylim(0, 0.05)
            stripbox_v_mut_light_g.set_ylim(0, 0.03)
            violin_v_mut_heavy_g.set_ylim(0, 0.12)
            violin_v_mut_light_g.set_ylim(0, 0.13)
            violin_v_mut_heavy_g.set_yticks(np.arange(0, 0.13, 0.01))
            violin_v_mut_light_g.set_yticks(np.arange(0, 0.14, 0.01))

    stripbox_v_mut_heavy_df.to_csv(
        mectric_outdir / f"figA_{src_type}_v_heavy_{seq_type}_mutation.csv",
        index=False,
    )
    violin_v_mut_heavy_df.to_csv(
        mectric_outdir / f"figC_{src_type}_v_heavy_{seq_type}_mutation.csv",
        index=False,
    )
    stripbox_v_mut_light_df.to_csv(
        mectric_outdir / f"figB_{src_type}_v_light_{seq_type}_mutation.csv",
        index=False,
    )
    violin_v_mut_light_df.to_csv(
        mectric_outdir / f"figD_{src_type}_v_light_{seq_type}_mutation.csv",
        index=False,
    )

    # TODO: create median of medians for each group as mutation_medians.csv
    if src_type == "prime":
        cols = ["trial", "pseudogroup", "weeks"]
    elif src_type == "boost":
        cols = ["trial", "pseudogroup"]
    heavy = stripbox_v_mut_heavy_df.groupby(cols)[f"v_mutation_{heavy}"].median()
    light = stripbox_v_mut_light_df.groupby(cols)[f"v_mutation_{light}"].median()
    median_df = pd.concat([heavy, light], axis=1)
    median_df = median_df.round(3)
    if src_type == "boost":
        median_df = median_df.reset_index()
        median_df = data.populate_psname(median_df)
    median_df.to_csv(
        mectric_outdir / f"figAB_medians.csv",
        index=True,
    )

    return (
        stripbox_v_mut_heavy_g,
        stripbox_v_mut_light_g,
        violin_v_mut_heavy_g,
        violin_v_mut_light_g,
    )


def plot_key_mutations_boost(
    data: Data,
    metric_outdir: Path,
    method: Literal[
        "linear",  # Linear interpolation between closest ranks
        "lower",  # Use the lower of the two values
        "higher",  # Use the higher of the two values
        "midpoint",  # Use the average of the two values
        "nearest",  # Use the closest value
        "inverted_cdf",  # Inverse CDF interpolation
        "hazen",  # (rank - 0.5)/nobs interpolation
        "weibull",  # (rank)/(nobs + 1) interpolation
        "median_unbiased",  # Type 8 estimate - median unbiased regardless of distribution
        "normal_unbiased",  # Type 7 estimate - unbiased if data is normally distributed
    ] = "median_unbiased",
) -> list[pw.Brick]:
    """Plot Key Mutations"""
    apply_global_font_settings()
    minimum_set = [
        "12A12",
        "12A21",
        "N6",
        "VRC27",
        "N17",
        "N60P1.1",
        "N60P25.1",
        "N60P23",
        "PCIN63_71I",
        "PCIN63_66B",
        "PCIN63_71G",
        "NIH45-46",
        "VRC07b",
        "VRC23",
        "VRC01",
        "VRC02",
        "VRC18",
        "VRC08",
        "VRC-PG19",
    ]
    g002_seqs = data.get_g002_sequences_boost()
    vrc01_seqs = data.get_vrc01_class_bnabs()
    palette = data.get_week_palette()
    vrc01_seqs = vrc01_seqs[vrc01_seqs["sequence_id"].isin(minimum_set)].copy()

    def _get_top_n_percet_cottrell(df, percent=0.9):
        """Groupby funciton to get qunatile of key residues"""
        return pd.Series(
            {"residues": df["cottrell_focused_v_common_score"].quantile(percent, interpolation="midpoint")}
        )

    def _get_top_n_percet_hcdr2(df, percent=0.9):
        """Groupby funciton to get qunatile of key residues"""
        return pd.Series({"residues": df["num_hcdr2_mutations"].quantile(percent, interpolation="midpoint")})

    # plottable = (
    #     g002_seqs.query("pseudogroup!=2")
    #     .groupby(["pubID", "weeks", "pseudogroup"])
    #     .apply(_get_top_n_percet_cottrell)
    #     .reset_index()
    # )
    plottable = (
        g002_seqs.query("pseudogroup!=2")
        .groupby(["pubID", "weeks", "pseudogroup"])["cottrell_focused_v_common_score"]
        .apply(
            lambda x: np.quantile(
                x,
                0.9,
                method=method,
            )
        )
        .reset_index()
    )
    plottable["residues"] = plottable["cottrell_focused_v_common_score"]

    ylabel = "90th percentile\nnumber of key VRC01-class\nHC residues"
    yname = ylabel.replace("\n", "_").replace(" ", "_")
    data.populate_psname(plottable).to_csv(metric_outdir / f"figE_{yname}.csv", index=False)
    print(metric_outdir / f"figE_{yname}.csv")
    plottable_hcdr2 = (
        g002_seqs.query("pseudogroup!=2")
        .groupby(["pubID", "weeks", "pseudogroup"])["num_hcdr2_mutations"]
        .apply(
            lambda x: np.quantile(
                x,
                0.9,
                method=method,
            )
        )
        .reset_index()
    )
    plottable_hcdr2["residues"] = plottable_hcdr2["num_hcdr2_mutations"]

    ylabel = "90th percentile\nnumber of key VRC01-class\nHCDR2 residues"
    yname = ylabel.replace("\n", "_").replace(" ", "_")
    data.populate_psname(plottable_hcdr2).to_csv(metric_outdir / f"figF_{yname}.csv", index=False)

    axes = [
        pw.Brick(figsize=(4, 4)),
        pw.Brick(figsize=(1, 4)),
        pw.Brick(figsize=(4, 4)),
        pw.Brick(figsize=(1, 4)),
    ]

    sns.boxplot(
        data=plottable,
        x="pseudogroup",
        y="residues",
        hue="weeks",
        dodge=False,
        ax=axes[0],
        linewidth=1,
        whis=[10, 90],
        fliersize=0,
        palette=palette,
    )
    sns.stripplot(
        data=plottable,
        x="pseudogroup",
        y="residues",
        hue="weeks",
        dodge=False,
        ax=axes[0],
        linewidth=1,
        palette=palette,
        size=7,
        alpha=0.8,
        edgecolor="black",
    )
    adjust_boxplot(axes[0])
    axes[0].legend_.remove()  # type: ignore
    axes[0].set_ylabel("90th percentile\nnumber of key VRC01-class\nHC residues", labelpad=14)
    axes[0].set_ylim(0, 8)
    axes[0].yaxis.set_major_locator(mtick.MultipleLocator(2))
    axes[0].yaxis.set_minor_locator(mtick.MultipleLocator(1))
    axes[0].set_xticklabels([])
    axes[0].set_xlabel("")
    axes[0].set_xticklabels(
        [
            r"eOD",
            r"eOD$\rightarrow$eOD",
            r"eOD$\rightarrow$core",
            r"eOD$\rightarrow$eOD",
            r"eOD$\rightarrow$core",
            r"eOD$\rightarrow$eOD$\rightarrow$core",
        ],
        rotation=45,
        verticalalignment="top",
        horizontalalignment="right",
    )

    # here is the control plot
    y = "cottrell_focused_v_common_score"
    sns.boxplot(
        data=vrc01_seqs,
        x="is_vrc01_class",
        y=y,
        color="#F2F0F2",
        linewidth=1,
        ax=axes[1],
        fliersize=0,
        whis=[10, 90],
    )
    sns.stripplot(
        data=vrc01_seqs,
        x="is_vrc01_class",
        y=y,
        color="#F2F0F2",
        linewidth=1,
        alpha=0.8,
        edgecolor="black",
        size=7,
        ax=axes[1],
    )
    adjust_boxplot(axes[1])
    axes[1].set_title("Control")
    axes[1].set_ylabel("")
    axes[1].set_ylim(0, 18)
    axes[1].yaxis.set_major_locator(mtick.MultipleLocator(4))
    axes[1].yaxis.set_minor_locator(mtick.MultipleLocator(2))
    axes[1].set_xlabel("")
    axes[1].set_xticklabels(["VRC01-class\nbnAbs"])

    # # hcdr2
    sns.boxplot(
        data=plottable_hcdr2,
        x="pseudogroup",
        y="residues",
        hue="weeks",
        dodge=False,
        ax=axes[2],
        linewidth=1,
        whis=[10, 90],
        fliersize=0,
        palette=palette,
    )
    sns.stripplot(
        data=plottable_hcdr2,
        x="pseudogroup",
        y="residues",
        hue="weeks",
        dodge=False,
        ax=axes[2],
        linewidth=1,
        palette=palette,
        size=7,
        alpha=0.8,
        edgecolor="black",
    )
    adjust_boxplot(axes[2])
    axes[2].set_ylabel(
        "90th percentile\nnumber of key VRC01-class\nHCDR2 residues",
        labelpad=14,
    )
    axes[2].set_xlabel("")
    axes[2].yaxis.set_major_locator(mtick.MultipleLocator(1))
    axes[2].yaxis.set_minor_locator(mtick.MultipleLocator(0.5))
    axes[2].set_xticklabels(
        [
            r"eOD",
            r"eOD$\rightarrow$eOD",
            r"eOD$\rightarrow$core",
            r"eOD$\rightarrow$eOD",
            r"eOD$\rightarrow$core",
            r"eOD$\rightarrow$eOD$\rightarrow$core",
        ],
        rotation=45,
        verticalalignment="top",
        horizontalalignment="right",
    )
    axes[2].legend_.set_visible(False)  # type: ignore
    plot_legend(
        pallete={"wk 8": "gold", "wk 16": "#E377C2", "wk 24": "#2078B4"},
        ax=axes[0],
        bbox_to_anchor=(0.5, -0.6),
    )
    plot_legend(
        pallete={"wk 8": "gold", "wk 16": "#E377C2", "wk 24": "#2078B4"},
        ax=axes[2],
        bbox_to_anchor=(0.5, -0.6),
    )

    # here is the control plot
    y = "num_hcdr2_mutations"
    sns.boxplot(
        data=vrc01_seqs,
        x="is_vrc01_class",
        y=y,
        color="#F2F0F2",
        linewidth=1,
        ax=axes[3],
        fliersize=0,
        whis=[10, 90],
    )
    sns.stripplot(
        data=vrc01_seqs,
        x="is_vrc01_class",
        y=y,
        color="#F2F0F2",
        linewidth=1,
        alpha=0.8,
        edgecolor="black",
        size=7,
        ax=axes[3],
    )
    adjust_boxplot(axes[3])

    axes[3].yaxis.set_major_locator(mtick.MultipleLocator(1))
    axes[3].yaxis.set_minor_locator(mtick.MultipleLocator(0.5))
    axes[3].set_xlabel("")
    axes[3].set_ylabel("")
    axes[3].set_xticklabels(["VRC01-class\nbnAbs"])
    axes[3].set_xlabel("")
    axes[3].set_title("Control")

    for ax in axes:
        sns.despine(ax=ax)

    axes[0].set_title("E", fontsize=20, fontweight="bold", x=-0.2, y=1)
    axes[2].set_title("F", fontsize=20, fontweight="bold", x=-0.2, y=1)
    plt.tight_layout()
    return axes  # figure.savefig(outpath + ".png", dpi=300, transparent=False)


def run_90_percentile_hc_residues(
    data: Data,
    metric_outdir: Path,
    method: Literal[
        "linear",  # Linear interpolation between closest ranks
        "lower",  # Use the lower of the two values
        "higher",  # Use the higher of the two values
        "midpoint",  # Use the average of the two values
        "nearest",  # Use the closest value
        "inverted_cdf",  # Inverse CDF interpolation
        "hazen",  # (rank - 0.5)/nobs interpolation
        "weibull",  # (rank)/(nobs + 1) interpolation
        "median_unbiased",  # Type 8 estimate - median unbiased regardless of distribution
        "normal_unbiased",  # Type 7 estimate - unbiased if data is normally distributed
    ] = "median_unbiased",
):
    apply_global_font_settings()

    def _get_top_n_percet(df, percent=0.9):
        """Groupby funciton to get qunatile of key residues"""
        # TODO: Why have we always used midpoint?
        return pd.Series(
            {
                "residues": df["cottrell_focused_v_common_score"].quantile(
                    # method="median unbiased",
                    q=percent,
                    interpolation="midpoint",
                )
            }
        )

    minimum_set = [
        "12A12",
        "12A21",
        "N6",
        "VRC27",
        "N17",
        "N60P1.1",
        "N60P25.1",
        "N60P23",
        "PCIN63_71I",
        "PCIN63_66B",
        "PCIN63_71G",
        "NIH45-46",
        "VRC07b",
        "VRC23",
        "VRC01",
        "VRC02",
        "VRC18",
        "VRC08",
        "VRC-PG19",
    ]
    pallete = data.get_trial_g001_g002_g003_palette()
    week_ticks = ["4", "8", "10", "16", "24/21"]

    vrc01_seqs = data.get_vrc01_class_bnabs()
    vrc01_seqs = vrc01_seqs[vrc01_seqs["sequence_id"].isin(minimum_set)].copy()

    set_letters = True
    metric_base_cols = ["pubID", "trial", "pseudogroup", "weeks"]

    axes = [pw.Brick(figsize=(4, 3)), pw.Brick(figsize=(1.6, 3))]

    combined = data.get_g00x_sequences_prime().copy(deep=True)

    combined = combined[~combined.is_vrc01_class.isna()]
    combined = combined.query("is_vrc01_class")

    combined["weeks"] = combined["weeks"].astype(int)
    combined = combined.query("weeks > 0")

    # Top 90th percentile of key residues
    plottable = (
        combined.groupby(metric_base_cols)["cottrell_focused_v_common_score"]
        .apply(
            lambda x: np.quantile(
                x,
                0.9,
                method=method,
            )
        )
        .reset_index()
    )
    plottable["residues"] = plottable["cottrell_focused_v_common_score"]

    ylabel = "90th percentile number of key VRC01-class HC residues".replace(" ", "_")
    plottable.to_csv(metric_outdir / f"figE_{ylabel}.csv", index=False)
    # sync trials weeks
    plottable.loc[(plottable["trial"] == "G002") & (plottable["weeks"] == 24), "weeks"] = 21
    g = sns.boxplot(
        data=plottable,
        x="weeks",
        y="residues",
        hue="trial",
        dodge=True,
        ax=axes[0],
        hue_order=pallete.keys(),
        linewidth=1,
        whis=[10, 90],
        fliersize=0,
        palette=pallete,
    )
    g = sns.stripplot(
        data=plottable,
        x="weeks",
        y="residues",
        hue="trial",
        dodge=True,
        ax=axes[0],
        # hue_order=["G002", "G003", "CFHR", "Aurum", "NAC"],
        hue_order=pallete.keys(),
        linewidth=1,
        jitter=0.1,
        palette=pallete,
        size=5,
        alpha=0.8,
        edgecolor="black",
    )
    adjust_boxplot(axes[0], median_linewidth=2.5)
    axes[0].legend_.remove()  # type: ignore

    # here is the control plot
    y = "cottrell_focused_v_common_score"
    g = sns.boxplot(
        data=vrc01_seqs,
        x="is_vrc01_class",
        y=y,
        color="#F2F0F2",
        linewidth=1,
        ax=axes[1],
        fliersize=0,
        whis=[10, 90],
    )
    g = sns.stripplot(
        data=vrc01_seqs,
        x="is_vrc01_class",
        y=y,
        color="#F2F0F2",
        linewidth=1,
        alpha=0.8,
        edgecolor="black",
        size=7,
        ax=axes[1],
    )

    adjust_boxplot(axes[1])

    axes[0].set_ylabel(
        r"$90^{th}$" + "percentile\nnumber of key VRC01-class\nHC residues",
        labelpad=14,
        # fontsize=yaxis_fontsize,
    )
    axes[0].tick_params(axis="x", labelsize=12)
    axes[0].tick_params(axis="y", labelsize=12)
    axes[0].set_xlabel("Weeks post vaccination")
    # week_ticks = sorted(plottable.weeks.unique())
    week_ticks = ["4", "8", "10", "16", "24/21"]
    axes[0].set_xticklabels(week_ticks)
    axes[0].set_ylim(0, 5)
    axes[0].yaxis.set_major_locator(mtick.MultipleLocator(2))
    axes[0].yaxis.set_minor_locator(mtick.MultipleLocator(1))
    axes[1].set_title("Control")
    # axes[1].set_yticklabels([])
    axes[1].set_ylabel("")
    axes[1].set_ylim(0, 18)
    axes[1].yaxis.set_major_locator(mtick.MultipleLocator(4))
    axes[1].yaxis.set_minor_locator(mtick.MultipleLocator(2))
    axes[1].set_xlabel("")
    axes[1].set_xticklabels(["VRC01-class\nbnAbs"])

    plot_legend(
        pallete=pallete,
        ax=axes[1],
        loc="center",
        bbox_to_anchor=(1.3, 0.5),
        n=1,
        labelspacing=0.5,
        fontsize=16,
    )
    axes[1].tick_params(axis="x")
    axes[1].tick_params(axis="y")
    sns.despine(ax=axes[0])
    sns.despine(ax=axes[1])

    if set_letters:
        axes[0].set_title("E", fontsize=20, fontweight="bold", x=-0.2, y=1)

    plottable = unsync_weeks(plottable)

    return (
        pw.hstack(axes[0], axes[1], margin=0),
        plottable[metric_base_cols + ["residues"]],
    )
