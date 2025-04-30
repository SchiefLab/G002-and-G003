import contextlib
from typing import Any, List, Literal, Tuple, Union

import numpy as np
import pandas as pd
import patchworklib as pw
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import ticker as mtick
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# plt type for type hinting
from matplotlib.pyplot import Figure

from g00x_figures.data import Data, calculate_resonse
from g00x_figures.plot_helpers.boxplot import adjust_boxplot, format_y_axis
from g00x_figures.plot_helpers.font import apply_global_font_settings
from g00x_figures.plot_helpers.legend import plot_legend

pw.overwrite_axisgrid()

metric_base_cols = ["pubID", "trial", "weeks"]
week_ticks = ["-5", "4", "8", "10", "16", "24/21"]
week_order = [-5, 4, 8, 10, 16, 21]
yaxis_fontsize = 16
xaxis_fontsize = 16
letter_fontsize = 20

apply_global_font_settings()
data = Data()
pallete = data.get_trial_g001_g002_g003_palette()
outdir = data.paths.figure_outdir


def get_mutational_group(seq: pd.DataFrame, metric: str) -> pd.DataFrame:
    mutational_seq = (
        seq.groupby(["pubID", "weeks", "trial"])[metric]
        .median()
        .reset_index()
        # .assign(trial=trial)
    )
    return mutational_seq


def patch(func: callable, figsize: tuple = (2, 2), *agrs, **kwargs):
    """Overwrite plot axes with patchworklib Brick

    Parameters
    ----------
    func :
        plot function to be patched
    """
    ax = pw.Brick(figsize=figsize)
    print(ax)
    # replace plot axes with patchworklib Brick
    g = func(*agrs, **kwargs, ax=ax)
    return g


def plot_stripbox_v_mutation_groups(
    data: Any,
    y: Literal[
        "v_mutation_aa_heavy",
        "v_mutation_aa_light",
        "v_mutation_heavy",
        "v_mutation_light",
    ],
    ax: Figure | pw.Brick,
    plot_type: Literal["stripbox", "violin"] = "stripbox",
    legend: bool = False,
):
    """Plot stripbox plot for multiple mutational groups"""

    g00x_seq = data.get_g00x_sequences_prime(sync_weeks=True).copy(deep=True)
    g00x_seq = g00x_seq.query("weeks != -5").query("is_vrc01_class==True").query("top_c_call=='IGHG'")
    seq_mut = get_mutational_group(g00x_seq, y)
    pallete = data.get_trial_g001_g002_g003_palette()
    df = pd.DataFrame()

    if plot_type == "stripbox":
        df = seq_mut
        g = sns.stripplot(
            x="weeks",
            y=y,
            hue="trial",
            edgecolor="black",
            linewidth=1,
            s=7,
            jitter=0.15,
            dodge=True,
            palette=pallete,
            hue_order=pallete.keys(),
            data=seq_mut,
            ax=ax,
        )
        g = sns.boxplot(
            x="weeks",
            y=y,
            hue="trial",
            dodge=True,
            fliersize=0,
            palette=pallete,
            hue_order=pallete.keys(),
            data=seq_mut,
            ax=ax,
        )
        adjust_boxplot(ax)
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=0, symbol="%"))
    elif plot_type == "violin":
        df = g00x_seq
        g = sns.violinplot(
            ax=ax,
            data=g00x_seq,
            bw=0.3,
            x="weeks",
            y=y,
            hue="trial",
            palette=pallete,
            inner="quartile",
            saturation=1,
            scale="area",
            hue_order=pallete.keys(),
            cut=0,
        )
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))

    ax.legend_.set_visible(False)
    ax.set_xlabel("")
    ax.set_xticklabels([])
    ax.tick_params(axis="x", labelsize=xaxis_fontsize)
    ax.tick_params(axis="y", labelsize=xaxis_fontsize)
    sns.despine(ax=ax)

    if legend:
        plot_legend(pallete=pallete, ax=ax, loc="upper center")

    return ax, df[metric_base_cols + [y]]


def _get_top_n_percet(df, percent=0.9):
    """Groupby funciton to get qunatile of key residues"""
    return pd.Series({"residues": df["cottrell_focused_v_common_score"].quantile(percent, interpolation="midpoint")})


def run_90_percentile_hc_residues():
    set_letters = True
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

    axes = [pw.Brick(figsize=(4, 4)), pw.Brick(figsize=(1.6, 4))]

    combined = g00x_seq = data.get_g00x_sequences_prime(sync_weeks=True).copy(deep=True)
    combined = combined[~combined.is_vrc01_class.isna()]
    combined = combined.query("is_vrc01_class")

    combined["weeks"] = combined["weeks"].astype(int)
    combined = combined.query("weeks > 0")

    plottable = combined.groupby(["pubID", "trial", "weeks"]).apply(_get_top_n_percet).reset_index()

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
        palette=pallete,
        size=7,
        alpha=0.8,
        edgecolor="black",
    )
    adjust_boxplot(axes[0])
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
        fontsize=yaxis_fontsize,
    )
    axes[0].tick_params(axis="x", labelsize=12)
    axes[0].tick_params(axis="y", labelsize=12)
    axes[0].set_xlabel("Weeks Post Vaccination", fontsize=14)
    # week_ticks = sorted(plottable.weeks.unique())
    week_ticks = ["4", "8", "10", "16", "24/21"]
    axes[0].set_xticklabels(week_ticks, fontsize=xaxis_fontsize)
    axes[0].set_ylim(0, 5)
    axes[0].yaxis.set_major_locator(mtick.MultipleLocator(2))
    axes[0].yaxis.set_minor_locator(mtick.MultipleLocator(1))
    axes[1].set_title("Control", fontsize=yaxis_fontsize)
    # axes[1].set_yticklabels([])
    axes[1].set_ylabel("")
    axes[1].set_ylim(0, 18)
    axes[1].yaxis.set_major_locator(mtick.MultipleLocator(4))
    axes[1].yaxis.set_minor_locator(mtick.MultipleLocator(2))
    axes[1].set_xlabel("")
    axes[1].set_xticklabels(["VRC01-class\nbnAbs"], fontsize=xaxis_fontsize)

    plot_legend(
        pallete=pallete,
        ax=axes[0],
        loc="upper center",
        bbox_to_anchor=(0.5, 1.1),
    )
    axes[1].tick_params(axis="x", labelsize=xaxis_fontsize)
    axes[1].tick_params(axis="y", labelsize=xaxis_fontsize)
    sns.despine(ax=axes[0])
    sns.despine(ax=axes[1])

    if set_letters:
        axes[0].set_title("E", fontsize=letter_fontsize, fontweight="bold", x=-0.2, y=1)

    return axes[0] | axes[1], plottable[metric_base_cols + ["residues"]]

    return axes[0] | axes[1], plottable[metric_base_cols + ["residues"]]


def main(data, outdir: Path, seq_type: Literal["aa", "nt"]) -> None:
    metric_base_cols = ["pubID", "trial", "weeks"]
    week_ticks = ["4", "8", "10", "16", "24/21"]
    yaxis_fontsize = 16
    xaxis_fontsize = 16
    letter_fontsize = 20

    apply_global_font_settings()
    data = Data()
    pallete = data.get_trial_g001_g002_g003_palette()
    g003_site_palette = data.get_g003_site_palette()

    main_outdir = outdir / "Main/fig3"
    main_metric_outdir = outdir / "Main-Metrics/fig3"
    sup_outdir = outdir / "Sup/fig3"
    sup_metric_outdir = outdir / "Sup-Metrics/fig3"

    if seq_type == "aa":
        v_mut = "v_mutation_aa"
    elif seq_type == "nt":
        v_mut = "v_mutation"

    name = "v_mutation_aa" if seq_type == "aa" else "v_mutation_nt"

    (
        stripbox_v_mut_heavy_g,
        stripbox_v_mut_heavy_df,
    ) = patch(
        plot_stripbox_v_mutation_groups,
        figsize=(4, 4),
        data=data,
        y=f"{v_mut}_heavy",
    )
    (
        stripbox_v_mut_light_g,
        stripbox_v_mut_light_df,
    ) = patch(
        plot_stripbox_v_mutation_groups,
        figsize=(4, 4),
        data=data,
        y=f"{v_mut}_light",
    )
    violin_v_mut_heavy_g, violin_v_mut_heavy_df = patch(
        plot_stripbox_v_mutation_groups,
        figsize=(4, 4),
        data=data,
        plot_type="violin",
        y=f"{v_mut}_heavy",
        legend=True,
    )
    violin_v_mut_light_g, violin_v_mut_light_df = patch(
        plot_stripbox_v_mutation_groups,
        figsize=(4, 4),
        data=data,
        y=f"{v_mut}_light",
        plot_type="violin",
        legend=True,
    )

    # y-axis limits
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

    # y-axis labels
    stripbox_v_mut_heavy_g.set_ylabel(
        r"$\mathregular{V_H}$" + f" gene\n% mutation ({seq_type})",
        fontsize=yaxis_fontsize,
    )
    stripbox_v_mut_light_g.set_ylabel(
        r"$\mathregular{V_{K/L}}$" + f" gene\n% mutation ({seq_type})",
        fontsize=yaxis_fontsize,
    )
    violin_v_mut_heavy_g.set_ylabel(
        r"$\mathregular{V_H}$" + f" gene\n% mutation ({seq_type})",
        fontsize=yaxis_fontsize,
    )
    violin_v_mut_light_g.set_ylabel(
        r"$\mathregular{V_{K/L}}$" + f" gene\n% mutation ({seq_type})",
        fontsize=yaxis_fontsize,
    )

    violin_v_mut_heavy_g.set_xlabel("Weeks post vaccination", fontsize=xaxis_fontsize)
    violin_v_mut_light_g.set_xlabel("Weeks post vaccination", fontsize=xaxis_fontsize)

    stripbox_v_mut_heavy_g.set_title("A", fontsize=letter_fontsize, fontweight="bold", x=-0.2, y=1)
    stripbox_v_mut_light_g.set_title("B", fontsize=letter_fontsize, fontweight="bold", x=-0.2, y=1)
    violin_v_mut_heavy_g.set_title("C", fontsize=letter_fontsize, fontweight="bold", x=-0.2, y=1)
    violin_v_mut_light_g.set_title("D", fontsize=letter_fontsize, fontweight="bold", x=-0.2, y=1)

    violin_v_mut_heavy_g.set_xticklabels(["4", "8", "10", "16", "24/21"], fontsize=xaxis_fontsize)
    violin_v_mut_light_g.set_xticklabels(["4", "8", "10", "16", "24/21"], fontsize=xaxis_fontsize)
    num_key_vrc01_hc_g, num_key_vrc01_hc_df = run_90_percentile_hc_residues()

    stripbox_v_mut_heavy_g.minorticks_on()  # Add minor ticks
    stripbox_v_mut_light_g.minorticks_on()  # Add minor ticks
    violin_v_mut_heavy_g.minorticks_on()  # Add minor ticks
    violin_v_mut_light_g.minorticks_on()  # Add minor ticks
    num_key_vrc01_hc_g.minorticks_on()  # Add minor ticks

    stripbox_v_mut_heavy_df.to_csv(
        main_metric_outdir / f"percent_{name}_heavy_stripbox_G00X.csv",
        index=False,
    )
    stripbox_v_mut_light_df.to_csv(
        main_metric_outdir / f"percent_{name}_light_stripbox_G00X.csv",
        index=False,
    )
    violin_v_mut_heavy_df.to_csv(
        main_metric_outdir / f"percent_{name}_heavy_violin_G00X.csv",
        index=False,
    )
    violin_v_mut_light_df.to_csv(
        main_metric_outdir / f"percent_{name}_light_violin_G00X.csv",
        index=False,
    )

    num_key_vrc01_hc_df.to_csv(
        main_metric_outdir / "90_percentile_num_hc_residues_G00X.csv",
        index=False,
    )

    g1 = stripbox_v_mut_heavy_g | stripbox_v_mut_light_g
    g2 = violin_v_mut_heavy_g | violin_v_mut_light_g
    g3 = num_key_vrc01_hc_g
    g3 = pw.spacer(g3, 0.15) | g3 | pw.spacer(g3, 0.15)
    g = g1 / g2 / g3

    g.savefig(main_outdir / f"percent_{name}_G00x.png", dpi=700)
