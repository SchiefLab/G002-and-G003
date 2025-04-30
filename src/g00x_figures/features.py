from pathlib import Path

import logomaker
import numpy as np
import pandas as pd
import seaborn as sns
from Levenshtein import distance
from matplotlib import pyplot as plt
from matplotlib import ticker as mtick
from matplotlib.patches import Patch

from g00x_figures.box_and_scatter.flow_frequencies import adjust_boxplot
from g00x_figures.data import Data
from g00x_figures.plot_helpers.font import apply_global_font_settings

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

explicit_mapping = {
    "IGKV1-33": "#91FCC0",
    "IGKV1-5": "#E377C2",
    "IGKV3-20": "#2078B4",
    "IGLV2-14": "#17BFD0",
    "IGKV3-15": "#9567BD",
    "IGLV2-11": "#BDBD23",
    "IGLV2-23": "#FF7F0F",
    "Other": "white",
}

sort_order = [
    "IGKV1-33",
    "IGKV3-20",
    "IGKV1-5",
    "IGKV3-15",
    "IGLV2-14",
    "IGLV2-23",
    "IGLV2-11",
    "Other",
]


def plot_sequence_logo(
    ax: plt.Axes,
    show_left_spine: bool = False,
    ytitle: str | None = None,
    title: str | None = None,
    char_set: list[str] = ["ABCDEF"],
    xticks: list[int | None] = [],
    ytitle_font: int = 8,
):
    color_scheme = "NajafabadiEtAl2017"

    # make matrix for character set
    mat_df = logomaker.alignment_to_matrix(char_set)
    mat_df = logomaker.transform_matrix(mat_df, normalize_values=True)

    # make and style logo
    logomaker.Logo(mat_df, ax=ax, color_scheme=color_scheme, show_spines=False)
    if xticks:
        ax.xaxis.set_major_locator(mtick.FixedLocator(range(0, 5)))
    else:
        ax.set_xticks([])
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_xticklabels(xticks, rotation=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if show_left_spine:
        ax.spines["left"].set_visible(True)
    if title:
        ax.set_title(title)
    if ytitle:
        ax.set_ylabel(ytitle, fontsize=ytitle_font, rotation=0, ha="right", va="center")


def process_dekosky(dekosky_df: pd.DataFrame) -> pd.DataFrame:
    # add dekosky data
    dekosky_df.loc[
        dekosky_df[dekosky_df["v_call_top_light"] == "IGKV3D-15"].index,
        "v_call_top_light",
    ] = "IGKV3-15"

    # anything not in explicit mapping is other
    dekosky_df.loc[
        dekosky_df[~dekosky_df["v_call_top_light"].isin(explicit_mapping.keys())].index,
        "v_call_top_light",
    ] = "Other"

    deskosky_ready = (
        dekosky_df.groupby(["Replicate", "donor"])
        .apply(lambda x: x["v_call_top_light"].value_counts(normalize=True))
        .reset_index()
        .rename({"v_call_top_light": "light_plot"}, axis=1)
        .groupby(["level_2"])
        .mean()
        .drop("Replicate", axis=1)
        .reset_index()
    )
    deskosky_ready["x-axis"] = "Dekosky"
    return deskosky_ready


def plot_light_chain_usage(data: Data, outpath: str) -> None:
    """Plot light chain usage for VRC01, non-VRC01 and all and dekosky"""

    figure, axes = plt.subplots(
        1,
        4,
        figsize=(4.83, 2.5),
        gridspec_kw={
            "width_ratios": [1, 1, 0.5, 0.2],
            "hspace": 0.12,
            "wspace": 0.1,
            "bottom": 0.25,
            "top": 0.85,
            "left": 0.18,
            "right": 0.9,
        },
    )
    colors = [explicit_mapping[i] for i in sort_order]

    g001_seq = data.get_g001_sequences_prime()
    g002_seq = data.get_g002_sequences_prime()
    dekosky = data.get_dekosky_vh12()
    g002_seq_ig = g002_seq.query("top_c_call=='IGHG'")
    combined = pd.concat([g001_seq, g002_seq_ig]).reset_index(drop=True)

    # use just the gene and not the allele
    combined["light_plot"] = combined["v_call_top_light"].str.split("*").str.get(0)

    # change IGKV3-15 to IGKV3-15
    combined.loc[
        combined[combined["light_plot"] == "IGKV3D-15"].index,
        "light_plot",
    ] = "IGKV3-15"

    # if light plot not in explicit mapping, set to other
    combined.loc[
        combined[~combined["light_plot"].isin(explicit_mapping.keys())].index,
        "light_plot",
    ] = "Other"

    # do you at least have a vh12?
    combined.loc[
        combined[combined["v_call_top_heavy"].str.split("*").str.get(0) == "IGHV1-2"].index,
        "has_vh1-2",
    ] = True

    vrc01_class = combined.query("is_vrc01_class==True")
    non_vrc01_class = combined.query("is_vrc01_class==False").query("`has_vh1-2`==True")

    pivot_vrc01_class = (
        vrc01_class.groupby("trial")
        .apply(lambda x: x["light_plot"].value_counts(normalize=True))
        .reset_index()
        .set_index("trial")
        .loc[:, sort_order]
    )
    pivot_nonvrc01_class = (
        non_vrc01_class.groupby("trial")
        .apply(lambda x: x["light_plot"].value_counts(normalize=True))
        .reset_index()
        .pivot(index="trial", columns="level_1", values="light_plot")
        .loc[:, sort_order]
    )
    pivot_vrc01_class.plot(
        kind="bar",
        stacked=True,
        color=colors,
        width=0.85,
        linewidth=1,
        edgecolor="black",
        ax=axes[0],
    )
    pivot_nonvrc01_class.plot(
        kind="bar",
        stacked=True,
        color=colors,
        width=0.85,
        linewidth=1,
        edgecolor="black",
        ax=axes[1],
    )
    pivot_dekosky = process_dekosky(dekosky)
    pivot_dekosky = pivot_dekosky.pivot(index="x-axis", columns="level_2", values="light_plot").loc[:, sort_order]
    pivot_dekosky.plot(
        kind="bar",
        stacked=True,
        color=colors,
        width=0.85,
        linewidth=1,
        edgecolor="black",
        ax=axes[2],
    )

    for ax_index, ax in enumerate(axes):
        ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        if ax_index == 0:
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
            label = "% BCRs with VRC01-class \n" + r"bnAb $\mathregular{V_{K/L}}$"
            ax.set_ylabel(label, size=10)
            ax.set_title("VRC01-class")
        ax.set_xlabel("")
        if ax_index > 0:
            ax.set_ylabel("")
            ax.set_yticklabels([])
        if ax_index == 1:
            ax.set_title("Non-VRC01-class\nwith VH1-2")
        if ax_index < 2:
            ax.set_xticklabels(["G001", "G002"], rotation=0)
        if ax_index == 2:
            ax.set_xticklabels(["DeKosky\nVH1-2"], rotation=0)
            ax.set_title("Control")
    axes[0].legend_.remove()  # type: ignore
    axes[1].legend_.remove()  # type: ignore
    axes[2].legend_.remove()  # type: ignore

    custom_lines = []
    for label in sort_order:
        color = explicit_mapping[label]
        if label != "Other":
            crude_label = label[2:].replace("V", "")
        else:
            crude_label = label
        custom_lines.append(
            Patch(
                facecolor=color,
                edgecolor="black",
                linewidth=1,
                label=crude_label,
            )
        )
    axes[-1].legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=1,
        fontsize=8,
        bbox_to_anchor=(0.5, 1.05),
        labelspacing=0.1,
    )

    sns.despine()
    axes[-1].set_axis_off()
    figure.savefig(outpath + ".png", dpi=300, transparent=False)


def plot_light_chain_usage_boost(data: Data, outpath: str) -> pd.DataFrame:
    """Plot light chain usage for VRC01, non-VRC01 and all and dekosky"""

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(4.83, 3.5),
        gridspec_kw={
            "width_ratios": [1, 0.5, 0.2],
            # "hspace": 0.12,
            # "wspace": 0.1,
            # "bottom": 0.25,
            # "top": 0.85,
            # "left": 0.18,
            # "right": 0.9,
        },
    )
    colors = [explicit_mapping[i] for i in sort_order]

    dekosky = data.get_dekosky_vh12()

    combined = data.get_g002_sequences_boost()

    # use just the gene and not the allele
    combined["light_plot"] = combined["v_call_top_light"].str.split("*").str.get(0)

    # change IGKV3-15 to IGKV3-15
    combined.loc[
        combined[combined["light_plot"] == "IGKV3D-15"].index,
        "light_plot",
    ] = "IGKV3-15"

    # if light plot not in explicit mapping, set to other
    combined.loc[
        combined[~combined["light_plot"].isin(explicit_mapping.keys())].index,
        "light_plot",
    ] = "Other"

    # do you at least have a vh12?
    vrc01_class = combined.query("is_vrc01_class==True")

    pivot_vrc01_class = (
        vrc01_class.groupby("pseudogroup").apply(lambda x: x["light_plot"].value_counts(normalize=True)).reset_index()
    )
    pivot_vrc01_class = pivot_vrc01_class.pivot(index="pseudogroup", columns="level_1", values="light_plot").loc[
        :, sort_order
    ]

    pivot_vrc01_class.plot(
        kind="bar",
        stacked=True,
        color=colors,
        width=0.85,
        linewidth=1,
        edgecolor="black",
        ax=axes[0],
    )
    pivot_dekosky = process_dekosky(dekosky)
    pivot_dekosky = pivot_dekosky.pivot(index="x-axis", columns="level_2", values="light_plot").loc[:, sort_order]
    pivot_dekosky.plot(
        kind="bar",
        stacked=True,
        color=colors,
        width=0.85,
        linewidth=1,
        edgecolor="black",
        ax=axes[1],
    )

    for ax_index, ax in enumerate(axes):
        ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        if ax_index == 0:
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
            label = "% BCRs with VRC01-class \n" + r"bnAb $\mathregular{V_{K/L}}$"
            ax.set_ylabel(label)
            ax.set_title("VRC01-class")
            ax.set_xticklabels(
                [
                    r"eOD",
                    r"core",
                    r"eOD$\rightarrow$eOD",
                    r"eOD$\rightarrow$core",
                    r"eOD$\rightarrow$eOD",
                    r"eOD$\rightarrow$core",
                    r"eOD$\rightarrow$eOD$\rightarrow$core",
                ],
                rotation=45,
                fontsize=8,
                ha="right",
                rotation_mode="anchor",
            )
        ax.set_xlabel("")
        if ax_index > 0:
            ax.set_ylabel("")
            ax.set_yticklabels([])
        # if ax_index < 2:
        #     ax.set_xticklabels(["G001", "G002"], rotation=0)
        if ax_index == 1:
            ax.set_xticklabels(["DeKosky\nVH1-2"], rotation=0)
            ax.set_title("Control")
    axes[0].legend_.remove()  # type: ignore
    axes[1].legend_.remove()  # type: ignore
    # axes[2].legend_.remove()  # type: ignore

    custom_lines = []
    for label in sort_order:
        color = explicit_mapping[label]
        if label != "Other":
            crude_label = label[2:].replace("V", "")
        else:
            crude_label = label
        custom_lines.append(
            Patch(
                facecolor=color,
                edgecolor="black",
                linewidth=1,
                label=crude_label,
            )
        )
    axes[-1].legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=1,
        fontsize=8,
        bbox_to_anchor=(0.5, 1.05),
        labelspacing=0.1,
    )

    sns.despine()
    axes[-1].set_axis_off()
    plt.tight_layout()
    figure.savefig(outpath + ".png", dpi=300, transparent=False)
    return vrc01_class, pivot_vrc01_class


def plot_qe_on_light_chain_boost(data: Data, outpath: str) -> pd.DataFrame:
    """is there a Q or E at position 96"""
    oas_5_len = data.get_oas_5_len()
    g002_seq = data.get_g002_sequences_boost()
    g002_seq = g002_seq[~g002_seq["cdr3_aa_light"].isna()]
    g002_seq_ig = g002_seq.query("top_c_call=='IGHG'")

    # do you at least have a vh12?
    g002_seq_ig["has_QE_at_96"] = g002_seq_ig["cdr3_aa_light"].apply(lambda x: True if x[-2] in ["Q", "E"] else False)
    oas_5_len["has_QE_at_96"] = oas_5_len["cdr3_aa"].apply(lambda x: True if x[-2] in ["Q", "E"] else False)
    combined_vrc01 = g002_seq_ig.query("is_vrc01_class==True").reset_index(drop=True)
    pal = data.get_week_palette()
    figure, axes = plt.subplots(
        1,
        3,
        figsize=(4.83, 3.5),
        gridspec_kw={
            "width_ratios": [1, 0.5, 0.2],
            # "hspace": 0.12,
            # "wspace": 0.1,
            # "bottom": 0.25,
            # "top": 0.85,
            # "left": 0.18,
            # "right": 0.9,
        },
    )

    # First is VRC01 class
    def _determine_metric(df: pd.DataFrame) -> pd.Series:
        return pd.Series(
            {
                True: df["has_QE_at_96"].sum() / len(df),
                False: 1 - (df["has_QE_at_96"].sum() / len(df)),
            }
        )

    plottable_vrc01 = (
        combined_vrc01.groupby(["pubID", "pseudogroup", "weeks"])
        .apply(_determine_metric)
        .stack()
        .reset_index()
        .rename({"level_3": "has_QE_at_96", 0: "freq"}, axis=1)
        .query("has_QE_at_96==True")
    )
    sns.boxplot(
        data=plottable_vrc01,
        x="pseudogroup",
        y="freq",
        hue="weeks",
        linewidth=1,
        dodge=False,
        palette=pal,
        ax=axes[0],
        whis=[10, 90],
        fliersize=0,
    )
    sns.stripplot(
        data=plottable_vrc01,
        x="pseudogroup",
        y="freq",
        hue="weeks",
        dodge=False,
        linewidth=1,
        palette=pal,
        ax=axes[0],
        size=4,
        edgecolor="black",
    )
    adjust_boxplot(axes[0])

    control_df = (
        oas_5_len.groupby(["oas_subject", "oas_author"])
        .apply(lambda x: x["has_QE_at_96"].value_counts(normalize=True))
        .to_frame()
        .reset_index()
        .query("level_2")
        .assign(dose_group="control")
    )

    # plot 3 is oas
    sns.boxplot(
        data=control_df,
        x="dose_group",
        y="has_QE_at_96",
        linewidth=1,
        dodge=False,
        color="#F2F0F2",
        ax=axes[1],
        whis=[10, 90],
        fliersize=0,
    )
    sns.stripplot(
        data=control_df,
        x="dose_group",
        y="has_QE_at_96",
        linewidth=1,
        color="#F2F0F2",
        dodge=False,
        jitter=0.15,
        ax=axes[1],
        size=4,
        edgecolor="black",
    )
    adjust_boxplot(axes[1])

    for ax_index, ax in enumerate(axes):
        ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        if ax_index == 0:
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
            label = "% BCRs with E/Q at\nLCDR3 position 96"
            ax.set_ylabel(label, size=10)
            ax.set_title("VRC01-Class")
            ax.set_xticklabels(
                [
                    r"eOD",
                    r"core",
                    r"eOD$\rightarrow$eOD",
                    r"eOD$\rightarrow$core",
                    r"eOD$\rightarrow$eOD",
                    r"eOD$\rightarrow$core",
                    r"eOD$\rightarrow$eOD$\rightarrow$core",
                ],
                rotation=45,
                fontsize=8,
                ha="right",
                rotation_mode="anchor",
            )
        if ax_index == 1:
            ax.set_ylabel("")
            ax.set_yticklabels([])

        if ax_index == 1:
            ax.set_xticklabels(["OAS"], rotation=0)
            ax.set_xticklabels(["OAS 5aa L3"])
            ax.set_title("Control")
        ax.set_ylim(0, 1.03)
        ax.set_xlabel("")

    axes[0].legend_.remove()  # type: ignore
    sns.despine()
    axes[-1].set_axis_off()
    custom_lines = []
    pal = {
        8: "gold",
        16: "#E377C2",
        24: "#2078B4",
    }
    for label in pal:
        color = pal[label]
        crude_label = label
        custom_lines.append(
            Patch(
                facecolor=color,
                edgecolor="black",
                linewidth=1,
                label=crude_label,
            )
        )
    axes[-1].legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=1,
        fontsize=8,
        bbox_to_anchor=(0.5, 1.05),
        labelspacing=0.1,
        title="Weeks",
    )
    plt.tight_layout()
    figure.savefig(outpath + ".png", dpi=300, transparent=False)
    return combined_vrc01, plottable_vrc01


def plot_qe_on_light_chain(data: Data, outpath: str) -> None:
    """is there a Q or E at position 96"""
    oas_5_len = data.get_oas_5_len()
    g001_seq = data.get_g001_sequences_prime()
    g002_seq = data.get_g002_sequences_prime()
    g002_seq_ig = g002_seq_ig = g002_seq.query("top_c_call=='IGHG'")
    combined = pd.concat([g001_seq, g002_seq_ig]).reset_index(drop=True)

    # do you at least have a vh12?
    combined["has_QE_at_96"] = combined["cdr3_aa_light"].apply(lambda x: True if x[-2] in ["Q", "E"] else False)
    oas_5_len["has_QE_at_96"] = oas_5_len["cdr3_aa"].apply(lambda x: True if x[-2] in ["Q", "E"] else False)
    combined_vrc01 = combined.query("is_vrc01_class==True").reset_index(drop=True)
    pal = data.get_trial_palette()
    combined_nonvrc01 = combined.query("is_vrc01_class==False").query("has_5_len==True")
    figure, axes = plt.subplots(
        1,
        4,
        figsize=(4.83, 2.5),
        gridspec_kw={
            "width_ratios": [1, 1, 0.5, 0.2],
            "hspace": 0.12,
            "wspace": 0.1,
            "bottom": 0.25,
            "top": 0.85,
            "left": 0.18,
            "right": 0.9,
        },
    )

    # First is VRC01 class
    def _determine_metric(df: pd.DataFrame) -> pd.Series:
        return pd.Series(
            {
                True: df["has_QE_at_96"].sum() / len(df),
                False: 1 - (df["has_QE_at_96"].sum() / len(df)),
            }
        )

    plottable_vrc01 = (
        combined_vrc01.groupby(["pubID", "trial"])
        .apply(_determine_metric)
        .stack()
        .reset_index()
        .rename({"level_2": "has_QE_at_96", 0: "freq"}, axis=1)
        .query("has_QE_at_96==True")
    )
    plottable_nonvrc01 = (
        combined_nonvrc01.groupby(["pubID", "trial"])
        .apply(_determine_metric)
        .stack()
        .reset_index()
        .rename({"level_2": "has_QE_at_96", 0: "freq"}, axis=1)
        .query("has_QE_at_96==True")
    )
    sns.boxplot(
        data=plottable_vrc01,
        x="trial",
        y="freq",
        linewidth=1,
        palette=pal,
        ax=axes[0],
        whis=[10, 90],
        order=["G001", "G002"],
        hue_order=["G001", "G002"],
        fliersize=0,
    )
    sns.stripplot(
        data=plottable_vrc01,
        x="trial",
        y="freq",
        hue="trial",
        linewidth=1,
        palette=pal,
        ax=axes[0],
        order=["G001", "G002"],
        hue_order=["G001", "G002"],
        size=4,
        edgecolor="black",
    )
    adjust_boxplot(axes[0])

    # plot 2 is nonvrc01 that have 5 len
    sns.boxplot(
        data=plottable_nonvrc01,
        x="trial",
        y="freq",
        linewidth=1,
        palette=pal,
        ax=axes[1],
        whis=[10, 90],
        order=["G001", "G002"],
        hue_order=["G001", "G002"],
        fliersize=0,
    )
    sns.stripplot(
        data=plottable_nonvrc01,
        x="trial",
        y="freq",
        hue="trial",
        linewidth=1,
        jitter=0.20,
        palette=pal,
        ax=axes[1],
        order=["G001", "G002"],
        hue_order=["G001", "G002"],
        size=4,
        edgecolor="black",
    )
    adjust_boxplot(axes[1])

    control_df = (
        oas_5_len.groupby(["oas_subject", "oas_author"])
        .apply(lambda x: x["has_QE_at_96"].value_counts(normalize=True))
        .to_frame()
        .reset_index()
        .query("level_2")
        .assign(dose_group="control")
    )

    # plot 3 is oas
    sns.boxplot(
        data=control_df,
        x="dose_group",
        y="has_QE_at_96",
        linewidth=1,
        color="#F2F0F2",
        ax=axes[2],
        whis=[10, 90],
        fliersize=0,
    )
    sns.stripplot(
        data=control_df,
        x="dose_group",
        y="has_QE_at_96",
        linewidth=1,
        color="#F2F0F2",
        jitter=0.15,
        ax=axes[2],
        size=4,
        edgecolor="black",
    )
    adjust_boxplot(axes[2])

    for ax_index, ax in enumerate(axes):
        ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        if ax_index == 0:
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
            label = "% BCRs with E/Q at\nLCDR3 position 96"
            ax.set_ylabel(label, size=10)
            ax.set_title("VRC01-Class")
        if ax_index == 1:
            ax.set_title("Non-VRC01-Class\nwith 5aa LCDR3")
        if ax_index > 0:
            ax.set_ylabel("")
            ax.set_yticklabels([])

        if ax_index < 2:
            ax.set_xticklabels(["G001", "G002"], rotation=0)
        if ax_index == 2:
            ax.set_xticklabels(["OAS"], rotation=0)
            ax.set_xticklabels(["OAS 5aa L3"])
            ax.set_title("Control")
        ax.set_ylim(0, 1.03)
        ax.set_xlabel("")

    axes[0].legend_.remove()  # type: ignore
    axes[1].legend_.remove()  # type: ignore
    sns.despine()
    axes[-1].set_axis_off()
    figure.savefig(outpath + ".png", dpi=300, transparent=True)


def plot_light_chain_dist_boost(data: Data, outpath: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    # apply_global_font_settings(9)
    g002_seqs = data.get_g002_sequences_boost()
    palette = data.get_week_palette()
    vrc01_ref_airr = data.get_vrc01_class_bnabs()

    vrc01_ref_airr_kappa_seqs = vrc01_ref_airr.query("locus_light == 'IGK'")["cdr3_aa_light"].to_list()
    vrc01_ref_airr_lambda_seqs = vrc01_ref_airr.query("locus_light == 'IGL'")["cdr3_aa_light"].to_list()

    def find_lowest_len(x, ref) -> int:
        return min(list(map(lambda y: distance(x, y), ref)))

    def report_lcdr3(df: pd.DataFrame):
        b = df["distance_to_known_lcdr3"].value_counts().sort_index()
        c = (b / b.sum()).cumsum()
        c = c.reindex([float(i) for i in range(6)]).ffill().fillna(0.0).transpose()
        return c

    combined = g002_seqs.reset_index(drop=True).query("weeks != '-5'")
    figure, axes = plt.subplots(
        1,
        3,
        figsize=(4.83, 3.5),
        gridspec_kw={
            "width_ratios": [1, 1, 0.2],
            # "hspace": 0.12,
            # "wspace": 0.1,
            # "bottom": 0.25,
            # "top": 0.84,
            # "left": 0.18,
            # "right": 0.9,
        },
        sharey=True,
    )

    class_df_kappa = combined.query("locus_light == 'IGK'").query("is_vrc01_class").copy()
    class_df_lambda = combined.query("locus_light == 'IGL'").query("is_vrc01_class").copy()
    class_df_kappa["distance_to_known_lcdr3_kappa"] = class_df_kappa["cdr3_aa_light"].apply(
        lambda x: find_lowest_len(x, vrc01_ref_airr_kappa_seqs)
    )
    class_df_lambda["distance_to_known_lcdr3_lambda"] = class_df_lambda["cdr3_aa_light"].apply(
        lambda x: find_lowest_len(x, vrc01_ref_airr_lambda_seqs)
    )
    combined_df = pd.concat([class_df_kappa, class_df_lambda]).reset_index(drop=True)
    combined_df["distance_to_known_lcdr3"] = combined_df[
        ["distance_to_known_lcdr3_kappa", "distance_to_known_lcdr3_lambda"]
    ].min(axis=1)
    plottable = (
        combined_df.query("locus_light=='IGK'")
        .groupby(["pubID", "pseudogroup", "weeks"])
        .apply(report_lcdr3)
        .stack()
        .reset_index()
        .rename(
            {0: "frequency"},
            axis=1,
        )
        .query("distance_to_known_lcdr3 == 0")
    )
    kappa_seq_df = combined_df.query("locus_light=='IGK'")
    kappa_df = plottable.copy(deep=True)

    sns.stripplot(
        data=plottable,
        x="pseudogroup",
        hue="weeks",
        y="frequency",
        palette=palette,
        edgecolor="black",
        linewidth=1,
        dodge=False,
        ax=axes[0],
        size=4,
    )
    sns.boxplot(
        data=plottable,
        x="pseudogroup",
        y="frequency",
        hue="weeks",
        palette=palette,
        linewidth=1,
        dodge=False,
        ax=axes[0],
        whis=[10, 90],
        fliersize=0,
    )
    axes[0].legend_.remove()  # type: ignore
    adjust_boxplot(axes[0])

    plottable = (
        combined_df.query("locus_light=='IGL'")
        .groupby(["pubID", "pseudogroup", "weeks"])
        .apply(report_lcdr3)
        .stack()
        .reset_index()
        .rename(
            {0: "frequency"},
            axis=1,
        )
        .query("distance_to_known_lcdr3 == 0")
    )
    lambda_seq_df = combined_df.query("locus_light=='IGL'")
    lambda_df = plottable.copy(deep=True)
    # from IPython import embed

    # embed()
    sns.stripplot(
        data=plottable,
        x="pseudogroup",
        hue="weeks",
        y="frequency",
        palette=palette,
        edgecolor="black",
        linewidth=1,
        order=range(1, 8),
        dodge=False,
        ax=axes[1],
        size=4,
    )
    sns.boxplot(
        data=plottable,
        x="pseudogroup",
        y="frequency",
        hue="weeks",
        palette=palette,
        order=range(1, 8),
        linewidth=1,
        dodge=False,
        ax=axes[1],
        whis=[10, 90],
        fliersize=0,
    )
    axes[1].legend_.remove()  # type: ignore
    adjust_boxplot(axes[1])
    for ax_index, ax in enumerate(axes):
        if ax_index == 1:
            ax.set_ylim(-0.05, 1.05)
            ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
            ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        if ax_index == 0:
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
            label = "% BCRs with VRC01-class\nbnAb LCDR3"
            ax.set_ylabel(label, size=10)
            ax.set_title("VRC01-Class Kappa")
        if ax_index == 1:
            ax.set_title("VRC01-Class Lambda")
        if ax_index > 0:
            ax.set_ylabel("")
        #     ax.set_yticklabels([])

        if ax_index in [0, 1]:
            ax.set_xticklabels(
                [
                    r"eOD",
                    r"core",
                    r"eOD$\rightarrow$eOD",
                    r"eOD$\rightarrow$core",
                    r"eOD$\rightarrow$eOD",
                    r"eOD$\rightarrow$core",
                    r"eOD$\rightarrow$eOD$\rightarrow$core",
                ],
                rotation=45,
                fontsize=8,
                ha="right",
                rotation_mode="anchor",
            )
        ax.set_xlabel("")
    pal = {
        8: "gold",
        16: "#E377C2",
        24: "#2078B4",
    }
    custom_lines = []
    for label in pal:
        color = pal[label]
        crude_label = label
        custom_lines.append(
            Patch(
                facecolor=color,
                edgecolor="black",
                linewidth=1,
                label=crude_label,
            )
        )
    axes[-1].legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=1,
        fontsize=7,
        bbox_to_anchor=(0.5, 1.05),
        labelspacing=0.1,
        title="Weeks",
        title_fontsize=8,
    )
    axes[-1].set_axis_off()
    sns.despine()
    plt.tight_layout()

    figure.savefig(outpath + ".png", dpi=300, transparent=False)
    return kappa_seq_df, kappa_df, lambda_seq_df, lambda_df


def plot_light_chain_dist(data: Data, outpath: str) -> None:
    g001_seqs = data.get_g001_sequences_prime()
    g002_seqs = data.get_g002_sequences_prime()
    palette = data.get_trial_palette()
    oas_5_len = data.get_oas_5_len()
    vrc01_ref_airr = data.get_vrc01_class_bnabs()

    vrc01_ref_airr_kappa_seqs = vrc01_ref_airr.query("locus_light == 'IGK'")["cdr3_aa_light"].to_list()
    vrc01_ref_airr_lambda_seqs = vrc01_ref_airr.query("locus_light == 'IGL'")["cdr3_aa_light"].to_list()

    def find_lowest_len(x, ref) -> int:
        return min(list(map(lambda y: distance(x, y), ref)))

    def report_lcdr3(df: pd.DataFrame):
        b = df["distance_to_known_lcdr3"].value_counts().sort_index()
        c = (b / b.sum()).cumsum()
        c = c.reindex([float(i) for i in range(6)]).ffill().fillna(0.0).transpose()
        return c

    combined = pd.concat([g001_seqs, g002_seqs]).reset_index(drop=True).query("weeks != '-5'")
    figure, axes = plt.subplots(
        1,
        4,
        figsize=(4.83, 2.5),
        gridspec_kw={
            "width_ratios": [1, 1, 0.5, 0.2],
            "hspace": 0.12,
            "wspace": 0.1,
            "bottom": 0.25,
            "top": 0.84,
            "left": 0.18,
            "right": 0.9,
        },
    )

    def plot_partial(is_vrc01_class, ax):
        if is_vrc01_class:
            class_df_kappa = combined.query("locus_light == 'IGK'").query("is_vrc01_class").copy()
            class_df_lambda = combined.query("locus_light == 'IGL'").query("is_vrc01_class").copy()
        else:
            class_df_kappa = (
                combined.query("locus_light == 'IGK'").query("is_vrc01_class==False").query("has_5_len").copy()
            )
            class_df_lambda = (
                combined.query("locus_light == 'IGL'").query("is_vrc01_class==False").query("has_5_len").copy()
            )

        class_df_kappa["distance_to_known_lcdr3_kappa"] = class_df_kappa["cdr3_aa_light"].apply(
            lambda x: find_lowest_len(x, vrc01_ref_airr_kappa_seqs)
        )
        class_df_lambda["distance_to_known_lcdr3_lambda"] = class_df_lambda["cdr3_aa_light"].apply(
            lambda x: find_lowest_len(x, vrc01_ref_airr_lambda_seqs)
        )
        combined_df = pd.concat([class_df_kappa, class_df_lambda]).reset_index(drop=True)
        combined_df["distance_to_known_lcdr3"] = combined_df[
            ["distance_to_known_lcdr3_kappa", "distance_to_known_lcdr3_lambda"]
        ].min(axis=1)
        plottable = (
            combined_df.groupby(["pubID", "trial", "locus_light"])
            .apply(report_lcdr3)
            .stack()
            .reset_index()
            .rename(
                {0: "frequency"},
                axis=1,
            )
            .query("distance_to_known_lcdr3 == 0")
        )
        plottable["x-axis"] = plottable["trial"] + plottable["locus_light"]

        sns.stripplot(
            data=plottable,
            x="x-axis",
            hue="trial",
            y="frequency",
            palette=palette,
            edgecolor="black",
            linewidth=1,
            dodge=False,
            order=["G001IGK", "G001IGL", "G002IGK", "G002IGL"],
            ax=ax,
            size=4,
        )
        sns.boxplot(
            data=plottable,
            x="x-axis",
            y="frequency",
            hue="trial",
            palette=palette,
            linewidth=1,
            dodge=False,
            ax=ax,
            order=["G001IGK", "G001IGL", "G002IGK", "G002IGL"],
            whis=[10, 90],
            fliersize=0,
        )
        ax.legend_.remove()  # type: ignore
        adjust_boxplot(ax)

    plot_partial(True, axes[0])
    plot_partial(False, axes[1])

    class_df_kappa = oas_5_len.query("locus == 'IGK'").copy()
    class_df_lambda = oas_5_len.query("locus == 'IGL'").copy()
    class_df_kappa["distance_to_known_lcdr3_kappa"] = class_df_kappa["cdr3_aa"].apply(
        lambda x: find_lowest_len(x, vrc01_ref_airr_kappa_seqs)
    )
    class_df_lambda["distance_to_known_lcdr3_lambda"] = class_df_lambda["cdr3_aa"].apply(
        lambda x: find_lowest_len(x, vrc01_ref_airr_lambda_seqs)
    )
    combined_df = pd.concat([class_df_kappa, class_df_lambda]).reset_index(drop=True)
    combined_df["distance_to_known_lcdr3"] = combined_df[
        ["distance_to_known_lcdr3_kappa", "distance_to_known_lcdr3_lambda"]
    ].min(axis=1)
    plottable = (
        combined_df.groupby(["oas_subject", "locus", "oas_author"])
        .apply(report_lcdr3)
        .stack()
        .reset_index()
        .rename({0: "frequency"}, axis=1)
        .query("distance_to_known_lcdr3 == 0")
    )
    sns.stripplot(
        data=plottable,
        x="locus",
        y="frequency",
        edgecolor="black",
        hue="locus",
        linewidth=1,
        dodge=False,
        palette={"IGK": "#F2F0F2", "IGL": "#F2F0F2"},
        order=["IGK", "IGL"],
        ax=axes[2],
        size=4,
    )
    sns.boxplot(
        data=plottable,
        x="locus",
        y="frequency",
        linewidth=1,
        dodge=False,
        ax=axes[2],
        order=["IGK", "IGL"],
        whis=[10, 90],
        fliersize=0,
    )
    axes[2].legend_.remove()  # type: ignore
    adjust_boxplot(axes[-1])

    for ax_index, ax in enumerate(axes):
        ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        if ax_index == 0:
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
            label = "% BCRs with VRC01-class\nbnAb LCDR3"
            ax.set_ylabel(label, size=10)
            ax.set_title("VRC01-Class")
        if ax_index == 1:
            ax.set_title("Non-VRC01-Class\nwith 5aa LCDR3")
        if ax_index > 0:
            ax.set_ylabel("")
            ax.set_yticklabels([])

        if ax_index < 2:
            ax.set_xticklabels(
                ["G001\nIGK", "G001\nIGL", "G002\nIGK", "G002\nIGL"],
                rotation=0,
                size=8,
            )
        if ax_index == 2:
            ax.set_title("Control")
            ax.set_xticklabels(["OAS\nIGK", "OAS\nIGL"], rotation=0, size=8)
            ax.set_xlabel("")
        ax.set_ylim(0, 1.03)
        ax.set_xlabel("")
    sns.despine()
    axes[-1].set_axis_off()
    figure.savefig(outpath + ".png", dpi=300, transparent=True)


def plot_has_100b_boost(data: Data, outpath: str) -> pd.DataFrame:
    """Does the antibody have a W100b?"""
    oas_vh12 = data.get_oas_vh12()
    g002_seq = data.get_g002_sequences_boost()
    g002_seq_ig = g002_seq_ig = g002_seq.query("top_c_call=='IGHG'")
    combined = pd.concat([g002_seq_ig]).reset_index(drop=True)

    combined["has_w100b"] = combined["junction_aa_heavy"].str[-6] == "W"

    # briney sotosoto_briney_vh12["has_w100b"] = soto_briney_vh12["cdr3_aa"].str[-5] == "W"
    oas_vh12["has_w100b"] = oas_vh12["cdr3_aa"].str[-5] == "W"

    combined_vrc01 = combined.query("is_vrc01_class==True").reset_index(drop=True)
    pal = data.get_week_palette()
    figure, axes = plt.subplots(
        1,
        3,
        figsize=(4.83, 3.5),
        gridspec_kw={
            "width_ratios": [1, 0.5, 0.2],
        },
    )

    # First is VRC01 class
    def _determine_metric(df: pd.DataFrame) -> pd.Series:
        return pd.Series(
            {
                True: df["has_w100b"].sum() / len(df),
                False: 1 - (df["has_w100b"].sum() / len(df)),
            }
        )

    plottable_vrc01 = (
        combined_vrc01.groupby(["pubID", "pseudogroup", "weeks"])
        .apply(_determine_metric)
        .stack()
        .reset_index()
        .rename({"level_3": "has_w100b", 0: "freq"}, axis=1)
        .query("has_w100b==True")
    )
    sns.boxplot(
        data=plottable_vrc01,
        x="pseudogroup",
        y="freq",
        hue="weeks",
        dodge=False,
        linewidth=1,
        palette=pal,
        ax=axes[0],
        whis=[10, 90],
        # order=["G001", "G002"],
        # hue_order=["G001", "G002"],
        fliersize=0,
    )
    sns.stripplot(
        data=plottable_vrc01,
        x="pseudogroup",
        y="freq",
        hue="weeks",
        dodge=False,
        linewidth=1,
        palette=pal,
        ax=axes[0],
        # order=["G001", "G002"],
        # hue_order=["G001", "G002"],
        size=4,
        edgecolor="black",
    )
    adjust_boxplot(axes[0])

    control_df = (
        oas_vh12.groupby(["oas_subject", "oas_author"])
        .apply(lambda x: x["has_w100b"].value_counts(normalize=True))
        .reset_index()
        .assign(dose_group="control")
    )

    # plot 3 is oas
    sns.boxplot(
        data=control_df,
        x="dose_group",
        y=True,
        linewidth=1,
        dodge=False,
        color="#F2F0F2",
        ax=axes[1],
        whis=[10, 90],
        fliersize=0,
    )
    sns.stripplot(
        data=control_df,
        x="dose_group",
        y=True,
        dodge=False,
        linewidth=1,
        color="#F2F0F2",
        jitter=0.15,
        ax=axes[1],
        size=4,
        edgecolor="black",
    )
    adjust_boxplot(axes[1])
    pal = {
        8: "gold",
        16: "#E377C2",
        24: "#2078B4",
    }
    custom_lines = []
    for label in pal:
        color = pal[label]
        crude_label = label
        custom_lines.append(
            Patch(
                facecolor=color,
                edgecolor="black",
                linewidth=1,
                label=crude_label,
            )
        )
    axes[-1].legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=1,
        fontsize=8,
        bbox_to_anchor=(0.5, 1.05),
        labelspacing=0.1,
        title="Weeks",
    )

    for ax_index, ax in enumerate(axes):
        ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        if ax_index == 0:
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
            label = "% BCRs with\n" + r"$\mathregular{Trp_{103-5}}$ in HCDR3"
            ax.set_ylabel(label, size=10)
            ax.set_title("VRC01-Class")
            ax.set_xticklabels(
                [
                    r"eOD",
                    r"core",
                    r"eOD$\rightarrow$eOD",
                    r"eOD$\rightarrow$core",
                    r"eOD$\rightarrow$eOD",
                    r"eOD$\rightarrow$core",
                    r"eOD$\rightarrow$eOD$\rightarrow$core",
                ],
                rotation=45,
                fontsize=8,
                ha="right",
                rotation_mode="anchor",
            )
        if ax_index > 0:
            ax.set_ylabel("")
            ax.set_yticklabels([])

        if ax_index == 1:
            ax.set_xticklabels(["OAS VH1-2"])
            ax.set_title("Control")
        ax.set_ylim(0, 1.03)
        ax.set_xlabel("")

    axes[0].legend_.remove()  # type: ignore
    sns.despine()
    axes[-1].set_axis_off()
    plt.tight_layout()
    figure.savefig(outpath + ".png", dpi=300, transparent=False)
    return combined_vrc01, plottable_vrc01


def plot_has_100b(data: Data, outpath: str) -> None:
    """Does the antibody have a W100b?"""
    oas_vh12 = data.get_oas_vh12()
    g001_seq = data.get_g001_sequences_prime()
    g002_seq = data.get_g002_sequences_prime()
    g002_seq_ig = g002_seq_ig = g002_seq.query("top_c_call=='IGHG'")
    combined = pd.concat([g001_seq, g002_seq_ig]).reset_index(drop=True)

    combined["has_w100b"] = combined["junction_aa_heavy"].str[-6] == "W"

    # briney sotosoto_briney_vh12["has_w100b"] = soto_briney_vh12["cdr3_aa"].str[-5] == "W"
    oas_vh12["has_w100b"] = oas_vh12["cdr3_aa"].str[-5] == "W"

    combined_vrc01 = combined.query("is_vrc01_class==True").reset_index(drop=True)
    pal = data.get_trial_palette()
    combined_nonvrc01 = combined.query("is_vrc01_class==False").query("`has_vh1-2`==True").reset_index(drop=True)
    figure, axes = plt.subplots(
        1,
        4,
        figsize=(4.83, 2.5),
        gridspec_kw={
            "width_ratios": [1, 1, 0.5, 0.2],
            "hspace": 0.12,
            "wspace": 0.1,
            "bottom": 0.25,
            "top": 0.85,
            "left": 0.18,
            "right": 0.9,
        },
    )

    # First is VRC01 class
    def _determine_metric(df: pd.DataFrame) -> pd.Series:
        return pd.Series(
            {
                True: df["has_w100b"].sum() / len(df),
                False: 1 - (df["has_w100b"].sum() / len(df)),
            }
        )

    plottable_vrc01 = (
        combined_vrc01.groupby(["pubID", "trial"])
        .apply(_determine_metric)
        .stack()
        .reset_index()
        .rename({"level_2": "has_w100b", 0: "freq"}, axis=1)
        .query("has_w100b==True")
    )
    plottable_nonvrc01 = (
        combined_nonvrc01.groupby(["pubID", "trial"])
        .apply(_determine_metric)
        .stack()
        .reset_index()
        .rename({"level_2": "has_w100b", 0: "freq"}, axis=1)
        .query("has_w100b==True")
    )
    sns.boxplot(
        data=plottable_vrc01,
        x="trial",
        y="freq",
        linewidth=1,
        palette=pal,
        ax=axes[0],
        whis=[10, 90],
        order=["G001", "G002"],
        hue_order=["G001", "G002"],
        fliersize=0,
    )
    sns.stripplot(
        data=plottable_vrc01,
        x="trial",
        y="freq",
        hue="trial",
        linewidth=1,
        palette=pal,
        ax=axes[0],
        order=["G001", "G002"],
        hue_order=["G001", "G002"],
        size=4,
        edgecolor="black",
    )
    adjust_boxplot(axes[0])

    # plot 2 is nonvrc01 that have 5 len
    sns.boxplot(
        data=plottable_nonvrc01,
        x="trial",
        y="freq",
        linewidth=1,
        palette=pal,
        ax=axes[1],
        whis=[10, 90],
        order=["G001", "G002"],
        hue_order=["G001", "G002"],
        fliersize=0,
    )
    sns.stripplot(
        data=plottable_nonvrc01,
        x="trial",
        y="freq",
        hue="trial",
        linewidth=1,
        jitter=0.20,
        palette=pal,
        ax=axes[1],
        order=["G001", "G002"],
        hue_order=["G001", "G002"],
        size=4,
        edgecolor="black",
    )
    adjust_boxplot(axes[1])

    control_df = (
        oas_vh12.groupby(["oas_subject", "oas_author"])
        .apply(lambda x: x["has_w100b"].value_counts(normalize=True))
        .reset_index()
        .assign(dose_group="control")
    )

    # plot 3 is oas
    sns.boxplot(
        data=control_df,
        x="dose_group",
        y=True,
        linewidth=1,
        color="#F2F0F2",
        ax=axes[2],
        whis=[10, 90],
        fliersize=0,
    )
    sns.stripplot(
        data=control_df,
        x="dose_group",
        y=True,
        linewidth=1,
        color="#F2F0F2",
        jitter=0.15,
        ax=axes[2],
        size=4,
        edgecolor="black",
    )
    adjust_boxplot(axes[2])

    for ax_index, ax in enumerate(axes):
        ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        if ax_index == 0:
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
            label = "% BCRs with\n" + r"$\mathregular{Trp_{103-5}}$ in HCDR3"
            ax.set_ylabel(label, size=10)
            ax.set_title("VRC01-Class")
        if ax_index == 1:
            ax.set_title("Non-VRC01-Class\nwith VH1-2")
        if ax_index > 0:
            ax.set_ylabel("")
            ax.set_yticklabels([])

        if ax_index < 2:
            ax.set_xticklabels(["G001", "G002"], rotation=0)
        if ax_index == 2:
            ax.set_xticklabels(["OAS VH1-2"])
            ax.set_title("Control")
        ax.set_ylim(0, 1.03)
        ax.set_xlabel("")

    axes[0].legend_.remove()  # type: ignore
    axes[1].legend_.remove()  # type: ignore
    sns.despine()
    axes[-1].set_axis_off()
    figure.savefig(outpath + ".png", dpi=300, transparent=True)


def plot_key_mutations(data: Data, outpath: str) -> None:
    """Plot Key Mutations"""
    g001_seqs = data.get_g001_sequences_prime()
    g002_seqs = data.get_g002_sequences_prime()
    vrc01_seqs = data.get_vrc01_class_bnabs()
    palette = data.get_trial_palette()
    vrc01_seqs = vrc01_seqs[vrc01_seqs["sequence_id"].isin(minimum_set)].copy()

    def _get_top_n_percet(df, percent=0.9):
        """Groupby funciton to get qunatile of key residues"""
        return pd.Series(
            {"residues": df["cottrell_focused_v_common_score"].quantile(percent, interpolation="midpoint")}
        )

    combined = pd.concat([g001_seqs, g002_seqs]).reset_index(drop=True).query("is_vrc01_class")
    combined["weeks"] = combined["weeks"].astype(int)
    combined = combined.query("weeks > 0")
    plottable = combined.groupby(["pubID", "trial", "weeks"]).apply(_get_top_n_percet).reset_index()

    figure, axes = plt.subplots(
        1,
        2,
        figsize=(6, 3),
        gridspec_kw={
            "width_ratios": [1, 0.4],
        },
    )
    sns.boxplot(
        data=plottable,
        x="weeks",
        y="residues",
        hue="trial",
        dodge=True,
        ax=axes[0],
        linewidth=1,
        whis=[10, 90],
        hue_order=["G001", "G002"],
        fliersize=0,
        palette=palette,
    )
    sns.stripplot(
        data=plottable,
        x="weeks",
        y="residues",
        hue="trial",
        dodge=True,
        ax=axes[0],
        linewidth=1,
        palette=palette,
        size=7,
        alpha=0.8,
        hue_order=["G001", "G002"],
        edgecolor="black",
    )
    adjust_boxplot(axes[0])
    axes[0].legend_.remove()  # type: ignore

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

    axes[0].set_ylabel(
        r"$90^{th}$" + "percentile\nnumber of key VRC01-class\nHC residues",
        labelpad=14,
    )
    axes[0].set_xlabel("Weeks Post Vaccination")
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
    custom_lines = []
    for x in palette:
        custom_lines.append(Patch(facecolor=palette[x], edgecolor="black", linewidth=1, label=x))
    axes[0].legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=2,
        bbox_to_anchor=(0.5, 1.2),
        labelspacing=0.1,
    )
    sns.despine()
    plt.tight_layout()
    figure.savefig(outpath + ".png", dpi=300, transparent=False)


def plot_key_mutations_boost(data: Data, outpath: str) -> None:
    """Plot Key Mutations"""
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

    plottable = (
        g002_seqs.query("pseudogroup!=1")
        .groupby(["pubID", "weeks", "pseudogroup"])
        .apply(_get_top_n_percet_cottrell)
        .reset_index()
    )
    plottable_hcdr2 = (
        g002_seqs.query("pseudogroup!=1")
        .groupby(["pubID", "weeks", "pseudogroup"])
        .apply(_get_top_n_percet_hcdr2)
        .reset_index()
    )

    figure, axes = plt.subplots(
        1,
        5,
        figsize=(8.5, 4.8),
        gridspec_kw={
            "width_ratios": [1, 0.4, 0.4, 1, 0.4],
            "bottom": 0.33,
            "top": 0.9,
            "right": 0.95,
        },
    )
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

    custom_lines = []
    pal = {
        16: "#E377C2",
        24: "#2078B4",
    }
    for x in pal:
        custom_lines.append(Patch(facecolor=pal[x], edgecolor="black", linewidth=1, label=x))

    axes[0].legend(
        custom_lines,
        ["wk " + i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=1,
        ncol=3,
        bbox_to_anchor=(0.5, -0.4),
        labelspacing=0.1,
    )
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
        ha="right",
        rotation_mode="anchor",
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

    # turn off the buffer axis
    axes[2].set_axis_off()

    # # hcdr2
    sns.boxplot(
        data=plottable_hcdr2,
        x="pseudogroup",
        y="residues",
        hue="weeks",
        dodge=False,
        ax=axes[3],
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
        ax=axes[3],
        linewidth=1,
        palette=palette,
        size=7,
        alpha=0.8,
        edgecolor="black",
    )
    adjust_boxplot(axes[3])
    axes[3].set_ylabel(
        "90th percentile\nnumber of key VRC01-class\nHCDR2 residues",
        labelpad=14,
    )
    axes[3].set_xlabel("")
    axes[3].yaxis.set_major_locator(mtick.MultipleLocator(1))
    axes[3].yaxis.set_minor_locator(mtick.MultipleLocator(0.5))
    axes[3].set_xticklabels(
        [
            r"eOD",
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
    axes[3].legend(
        custom_lines,
        ["wk " + i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=1,
        ncol=3,
        bbox_to_anchor=(0.5, -0.4),
        labelspacing=0.1,
    )

    # here is the control plot
    y = "num_hcdr2_mutations"
    sns.boxplot(
        data=vrc01_seqs,
        x="is_vrc01_class",
        y=y,
        color="#F2F0F2",
        linewidth=1,
        ax=axes[4],
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
        ax=axes[4],
    )
    adjust_boxplot(axes[4])

    axes[4].yaxis.set_major_locator(mtick.MultipleLocator(1))
    axes[4].yaxis.set_minor_locator(mtick.MultipleLocator(0.5))
    axes[4].set_xlabel("")
    axes[4].set_ylabel("")
    axes[4].set_xticklabels(["VRC01-class\nbnAbs"])
    axes[4].set_xlabel("")
    axes[4].set_title("Control")

    sns.despine()
    plt.tight_layout()
    figure.savefig(outpath + ".png", dpi=300, transparent=False)


def plot_seq_logos_vrc01_class(data: Data, outpath: str) -> None:
    figure, axes = plt.subplots(4, 2, gridspec_kw={"wspace": 0.1, "hspace": 0.1, "left": 0.2})
    g001_seqs = data.get_g001_sequences_prime()
    g002_seqs = data.get_g002_sequences_prime()
    human_naive = data.get_human_naive_5_len_df().query("is_vrc01_class")
    kappa_vrc01 = data.get_kappa_vrc01_select_df()
    lambda_vrc01 = data.get_lambda_vrc01_select_df()

    # bnabs
    kappa_vrc01 = kappa_vrc01["cdr3_aa_no_gaps"].to_list()
    lambda_vrc01 = lambda_vrc01["cdr3_aa_no_gaps"].to_list()

    # g001 kappa and lambda
    g001_kappa = g001_seqs.query("is_vrc01_class").query("locus_light=='IGK'")["cdr3_aa_light"].to_list()
    g001_lambda = g001_seqs.query("is_vrc01_class").query("locus_light=='IGL'")["cdr3_aa_light"].to_list()

    # g002 kappa and lambda
    g002_kappa = g002_seqs.query("is_vrc01_class").query("locus_light=='IGK'")["cdr3_aa_light"].to_list()
    g002_lambda = g002_seqs.query("is_vrc01_class").query("locus_light=='IGL'")["cdr3_aa_light"].to_list()

    # human naive
    human_naive_kappa = human_naive.query("locus_light=='IGK'")["cdr3_aa_light"].to_list()
    human_naive_lambda = human_naive.query("locus_light=='IGL'")["cdr3_aa_light"].to_list()

    plot_sequence_logo(
        axes[0][0],
        char_set=kappa_vrc01,
        ytitle="bnAbs",
        ytitle_font=12,
        show_left_spine=True,
        title="Kappa",
    )
    plot_sequence_logo(axes[0][1], char_set=lambda_vrc01, title="Lambda")
    plot_sequence_logo(
        axes[1][0],
        char_set=g001_kappa,
        ytitle="G001",
        ytitle_font=12,
        show_left_spine=True,
    )
    plot_sequence_logo(axes[1][1], char_set=g001_lambda)
    plot_sequence_logo(
        axes[2][0],
        char_set=g002_kappa,
        ytitle="G002",
        ytitle_font=12,
        show_left_spine=True,
    )
    plot_sequence_logo(axes[2][1], char_set=g002_lambda)
    plot_sequence_logo(
        axes[3][0],
        char_set=human_naive_kappa,
        ytitle="Human naive\nGT8 binders",
        ytitle_font=12,
        show_left_spine=True,
        xticks=list(range(93, 98)),
    )
    plot_sequence_logo(axes[3][1], char_set=human_naive_lambda, xticks=list(range(93, 98)))
    axes[3][0].annotate(
        text="Light Chain Position",
        xy=(1, -0.5),
        xytext=(1.05, -0.4),
        xycoords="axes fraction",
        ha="center",
        va="top",
        fontsize=14,
    )
    figure.suptitle("VRC01-Class", x=0.55, fontsize=16)
    figure.savefig(outpath + ".png", dpi=300, transparent=True)


def plot_seq_logos_nonvrc01_class(data: Data, outpath: str) -> None:
    figure, axes = plt.subplots(4, 2, gridspec_kw={"wspace": 0.1, "hspace": 0.1, "left": 0.2})
    g001_seqs = data.get_g001_sequences_prime()
    g002_seqs = data.get_g002_sequences_prime()
    oas_5 = data.get_oas_5_len()
    kappa_vrc01 = data.get_kappa_vrc01_select_df()
    lambda_vrc01 = data.get_lambda_vrc01_select_df()

    # bnabs
    kappa_vrc01 = kappa_vrc01["cdr3_aa_no_gaps"].to_list()
    lambda_vrc01 = lambda_vrc01["cdr3_aa_no_gaps"].to_list()

    # g001 kappa and lambda
    g001_kappa = (
        g001_seqs.query("is_vrc01_class==False")
        .query("has_5_len")
        .query("locus_light=='IGK'")["cdr3_aa_light"]
        .to_list()
    )
    g001_lambda = (
        g001_seqs.query("is_vrc01_class==False")
        .query("has_5_len")
        .query("locus_light=='IGL'")["cdr3_aa_light"]
        .to_list()
    )

    # g002 kappa and lambda
    g002_kappa = (
        g002_seqs.query("is_vrc01_class==False")
        .query("has_5_len")
        .query("locus_light=='IGK'")["cdr3_aa_light"]
        .to_list()
    )
    g002_lambda = (
        g002_seqs.query("is_vrc01_class==False")
        .query("has_5_len")
        .query("locus_light=='IGL'")["cdr3_aa_light"]
        .to_list()
    )

    # human naive
    oas_kappa = oas_5.query("locus=='IGK'")["cdr3_aa"].to_list()
    oas_lambda = oas_5.query("locus=='IGL'")["cdr3_aa"].to_list()

    plot_sequence_logo(
        axes[0][0],
        char_set=kappa_vrc01,
        ytitle="bnAbs",
        ytitle_font=12,
        show_left_spine=True,
        title="Kappa",
    )
    plot_sequence_logo(axes[0][1], char_set=lambda_vrc01, title="Lambda")
    plot_sequence_logo(
        axes[1][0],
        char_set=g001_kappa,
        ytitle="G001",
        ytitle_font=12,
        show_left_spine=True,
    )
    plot_sequence_logo(axes[1][1], char_set=g001_lambda)
    plot_sequence_logo(
        axes[2][0],
        char_set=g002_kappa,
        ytitle="G002",
        ytitle_font=12,
        show_left_spine=True,
    )
    plot_sequence_logo(axes[2][1], char_set=g002_lambda)
    plot_sequence_logo(
        axes[3][0],
        char_set=oas_kappa,
        ytitle="OAS 5aa L3",
        ytitle_font=12,
        show_left_spine=True,
        xticks=list(range(93, 98)),
    )
    plot_sequence_logo(axes[3][1], char_set=oas_lambda, xticks=list(range(93, 98)))
    axes[3][0].annotate(
        text="Light Chain Position",
        xy=(1, -0.5),
        xytext=(1.05, -0.4),
        xycoords="axes fraction",
        ha="center",
        va="top",
        fontsize=14,
    )
    figure.suptitle("Non-VRC01-class 5aa LCDR3s", x=0.55, fontsize=16)
    figure.savefig(outpath + ".png", dpi=300, transparent=True)


def plot_isotype_data(seq: pd.DataFrame) -> None:
    """Plot isotype distributions for all sequences and VRC01-class sequences.

    Args:
        seq: DataFrame of G002 sequences from primary immunization
        seq_boost: DataFrame of G002 sequences from boost immunization
    """

    def assign_top_c_call_subtype(df):
        df["top_c_call_subtype"] = df["c_call_heavy"].str.split("*").str.get(0).fillna("")
        return df

    seq = assign_top_c_call_subtype(seq)
    # seq_boost = assign_top_c_call_subtype(seq_boost)

    # anything where the top c call is empty, we want to set it to UNK
    seq.loc[seq[seq["top_c_call"] == ""].index, "top_c_call"] = "UNK"
    seq.loc[seq[seq["top_c_call"] == "~"].index, "top_c_call"] = "UNK"
    # seq_boost.loc[
    #     seq_boost[seq_boost["top_c_call"] == ""].index,
    #     "top_c_call",
    # ] = "UNK"

    seq_igg = seq[seq["top_c_call_subtype"].str.startswith("IGHG")]

    def get_plotable(df, groupby="weeks", col: str = "top_c_call"):
        return (
            df.sort_values(col)
            .groupby([groupby] if isinstance(groupby, str) else groupby)
            .apply(lambda x: x[col].value_counts(normalize=True, dropna=False))
            .reset_index()
            .rename({"level_1": "isotype", col: "frequency"}, axis=1)
        )

    # all sequences
    plottable_1 = get_plotable(seq)

    # only vrc01 class
    plottable_2 = get_plotable(seq.query("is_vrc01_class"))

    # all sequences boosted
    plottable_3 = get_plotable(seq_igg, col="top_c_call_subtype", groupby="weeks")

    # only boosted vrc01 class
    plottable_4 = get_plotable(
        seq_igg.query("is_vrc01_class"),
        col="top_c_call_subtype",
        groupby="weeks",
    )

    plottable_1.name = "all_sequences"
    plottable_2.name = "vrc01_class"
    plottable_3.name = "all_sequences_boosted"
    plottable_4.name = "vrc01_class_boosted"

    plottable_1.key = ["isotype", "frequency"]
    plottable_2.key = ["isotype", "frequency"]
    plottable_3.key = ["isotype", "frequency"]
    plottable_4.key = ["isotype", "frequency"]

    dfs = [plottable_1, plottable_2, plottable_3, plottable_4]

    fontsize = 10
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))

    axes[0][0].letter = "A"
    axes[0][1].letter = "B"
    axes[1][0].letter = "C"
    axes[1][1].letter = "D"

    colors = [
        "#91FCC0",
        "#6C65FF",
        "#E377C2",
        "#9567BD",
        "#17BFD0",
        "#FF7F0F",
        "#BDBD23",
        "white",
    ]
    isotype_color_mapping = {
        "IGHA": "#9567BD",
        "IGHD": "#91FCC0",
        "IGHE": "#BDBD23",
        "IGHG": "#6C65FF",
        "IGHM": "#E377C2",
        "UNK": "white",
    }
    isotype_subtype_color_mapping = {
        "IGHG1": "#6C65FF",
        "IGHG2": "#FF7F0F",
        "IGHG3": "#17BFD0",
        "IGHG4": "#E377C2",
    }
    # plot all sequences
    (
        plottable_1.pivot("weeks", "isotype", "frequency")
        .fillna(0)
        .plot(
            kind="bar",
            stacked=True,
            edgecolor="black",
            color=isotype_color_mapping,
            # color=["purple", "green", "gold", "blue", "red", "grey"],
            ax=axes[0][0],
        )
    )
    (
        plottable_3.pivot("weeks", "isotype", "frequency")
        .fillna(0)
        .plot(
            kind="bar",
            stacked=True,
            edgecolor="black",
            color=isotype_subtype_color_mapping,
            ax=axes[0][1],
        )
    )
    (
        plottable_2.pivot("weeks", "isotype", "frequency")
        .fillna(0)
        .plot(
            kind="bar",
            stacked=True,
            edgecolor="black",
            color=isotype_color_mapping,
            ax=axes[1][0],
        )
    )
    # plottable_4 = plottable_4.pivot("weeks", "isotype", "frequency")
    # from IPython import embed

    # embed()
    # plottable_4.loc[-5] = [np.nan] * len(plottable_4.columns)
    # plottable_4 = plottable_4.sort_index()
    (
        plottable_4.pivot("weeks", "isotype", "frequency")
        .fillna(0)
        .plot(
            kind="bar",
            stacked=True,
            edgecolor="black",
            color=isotype_subtype_color_mapping,
            ax=axes[1][1],
        )
    )

    for i, ax in enumerate(axes.flatten()):
        ax.legend(bbox_to_anchor=(1, 1), prop={"size": fontsize})
        if i in [0, 2]:
            ax.set_ylabel("Frequency", fontsize=fontsize)
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=fontsize)
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, symbol=None))
        else:
            ax.set_ylabel("")
            ax.set_yticklabels([])

        ax.set_xlabel("")
        ax.set_xticklabels(
            [f"wk{week}" for week in [-5, 4, 8, 16, 24]],
            rotation=0,
            # verticalalignment="top",
            # horizontalalignment="right",
            fontsize=fontsize,
        )

    axes[0][0].set_title("All Sequences", fontsize=fontsize)
    axes[1][0].set_title("VRC01-class", fontsize=fontsize)
    axes[0][1].set_title("All Sequences", fontsize=fontsize)
    axes[1][1].set_title("VRC01-class", fontsize=fontsize)

    for ax in axes.flatten():
        ax.text(-0.2, 1.1, ax.letter, transform=ax.transAxes, size=14, weight="bold")

    sns.despine()
    plt.tight_layout(w_pad=-2)

    return fig, dfs


def plot_isotype_data_pseudogroups(seq_boost: Data) -> tuple[plt.Figure, list[pd.DataFrame]]:
    """
    TODO:
    17. Add two additional panels to fig. S12 (currently showing Ig isotype distributions post-GT8),
    to make a 2x2 figure. The new panels should go to the right of the current two panels.
    The new panels should show IgG subclass distributions for the same samples indicated in the current panels.
    """

    def assign_top_c_call_subtype(df) -> pd.DataFrame:
        df["top_c_call_subtype"] = df["c_call_heavy"].str.split("*").str.get(0).fillna("")
        return df

    seq_boost = assign_top_c_call_subtype(seq_boost)

    # anything where the top c call is empty, we want to set it to UNK
    seq_boost.loc[
        seq_boost[seq_boost["top_c_call"] == ""].index,
        "top_c_call",
    ] = "UNK"

    seq_boost_igg = seq_boost[seq_boost["top_c_call_subtype"].str.startswith("IGHG")]

    isotype_subtype_color_mapping = {
        "IGHG1": "#6C65FF",
        "IGHG2": "#FF7F0F",
        "IGHG3": "#17BFD0",
        "IGHG4": "#E377C2",
    }

    def get_plotable(df, groupby="weeks", col: str = "top_c_call"):
        return (
            df.sort_values(col)
            .groupby([groupby] if isinstance(groupby, str) else groupby)
            .apply(lambda x: x[col].value_counts(normalize=True, dropna=False))
            .reset_index()
            .rename({"level_1": "isotype", col: "frequency"}, axis=1)
        )

    # all sequences boosted
    plottable_1 = get_plotable(seq_boost.query("is_vrc01_class"), groupby="pseudogroup")

    # only boosted vrc01 class
    plottable_2 = get_plotable(
        seq_boost_igg.query("is_vrc01_class"),
        groupby="pseudogroup",
        col="top_c_call_subtype",
    )
    fontsize = 10
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    axes[0].letter = "A"
    axes[1].letter = "B"

    plottable_1.name = "vrc01_class"
    plottable_2.name = "vrc01_class_igg"

    plottable_1.key = ["isotype", "frequency"]
    plottable_2.key = ["isotype", "frequency"]

    dfs = [plottable_1, plottable_2]

    colors = [
        "#91FCC0",
        "#6C65FF",
        "#E377C2",
        "#9567BD",
        "#17BFD0",
        "#FF7F0F",
        "#BDBD23",
        "white",
    ]
    isotype_color_mapping = {
        "IGHA": "#9567BD",
        "IGHD": "#91FCC0",
        "IGHE": "#BDBD23",
        "IGHG": "#6C65FF",
        "IGHM": "#E377C2",
        "UNK": "white",
    }
    # plot all sequences
    (
        plottable_1.pivot("pseudogroup", "isotype", "frequency")
        .fillna(0)
        .plot(
            kind="bar",
            stacked=True,
            edgecolor="black",
            color=isotype_color_mapping,
            # color=["purple", "green", "gold", "blue", "red", "grey"],
            ax=axes[0],
        )
    )

    (
        plottable_2.pivot("pseudogroup", "isotype", "frequency")
        .fillna(0)
        .plot(
            kind="bar",
            stacked=True,
            edgecolor="black",
            color=isotype_subtype_color_mapping,
            ax=axes[1],
        )
    )

    axes[0].set_xticklabels(
        [
            "eOD" + "\nwk8",
            "core" + "\nwk8",
            r"eOD$\rightarrow$eOD" + "\nwk16",
            r"eOD$\rightarrow$core" + "\nwk16",
            r"eOD$\rightarrow$eOD" + "\nwk24",
            r"eOD$\rightarrow$core" + "\nwk24",
            r"eOD$\rightarrow$eOD$\rightarrow$core" + "\nwk24",
        ],
        rotation=45,
        ha="right",
        fontsize=8,
        rotation_mode="anchor",
    )
    axes[1].set_xticklabels(
        [
            "eOD" + "\nwk8",
            "core" + "\nwk8",
            r"eOD$\rightarrow$eOD" + "\nwk16",
            r"eOD$\rightarrow$core" + "\nwk16",
            r"eOD$\rightarrow$eOD" + "\nwk24",
            r"eOD$\rightarrow$core" + "\nwk24",
            r"eOD$\rightarrow$eOD$\rightarrow$core" + "\nwk24",
        ],
        rotation=45,
        ha="right",
        fontsize=8,
        rotation_mode="anchor",
    )
    for i, ax in enumerate(axes.flatten()):
        ax.legend(bbox_to_anchor=(1.1, 1.05), prop={"size": 8})
        if i in [0, 2]:
            ax.set_ylabel("Frequency", fontsize=fontsize)
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, symbol=None))
        else:
            ax.set_ylabel("")
            ax.set_yticklabels([])
        ax.set_xlabel("")

    axes[0].set_title("VRC01-class", fontsize=fontsize)
    axes[1].set_title("VRC01-class", fontsize=fontsize)

    sns.despine()
    # plt.tight_layout(h_pad=1.5)

    for ax in axes.flatten():
        ax.text(-0.2, 1.1, ax.letter, transform=ax.transAxes, size=14, weight="bold")

    plt.tight_layout(w_pad=-2)
    # sns.despine()
    # plt.tight_layout()

    return fig, dfs
    # fig.savefig(str(outpath) + ".png", dpi=300, transparent=False)
