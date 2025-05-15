from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Levenshtein import distance
from matplotlib import pyplot as plt
from matplotlib import ticker as mtick
from matplotlib.patches import Patch

from g00x_figures.data import Data, Transforms
from g00x_figures.g00x_plot_templates.bar_plots import stacked_bar_plot_pivot_df
from g00x_figures.plot_helpers.font import apply_global_font_settings
from g00x_figures.plots import Plot


def single_group_freqs(df: pd.DataFrame, pivot_on: str, values: str = "v_call_top_light_gene") -> pd.DataFrame:
    """Only return frequency DF"""
    return (
        df.groupby(pivot_on)
        .apply(lambda x: x[values].value_counts(normalize=True))
        .reset_index()
        .pivot(index=pivot_on, columns=f"level_1", values=values)
    )


def find_lowest_len(x, ref) -> int:
    return min(list(map(lambda y: distance(x, y), ref)))


def report_lcdr3(df: pd.DataFrame):
    b = df["distance_to_known_lcdr3"].value_counts().sort_index()
    c = (b / b.sum()).cumsum()
    c = c.reindex([float(i) for i in range(6)]).ffill().fillna(0.0).transpose()
    return c


def determine_metric(df: pd.DataFrame, metric: str) -> pd.Series:
    return pd.Series({True: df[metric].sum() / len(df), False: 1 - (df[metric].sum() / len(df))})


def bnab_lcdr3(seq_df: pd.DataFrame, vrc01_ref_airr: pd.DataFrame, is_vrc01_class: bool) -> pd.DataFrame:
    vrc01_ref_airr_kappa_seqs = vrc01_ref_airr.query("locus_light == 'IGK'")["cdr3_aa_light"].to_list()
    vrc01_ref_airr_lambda_seqs = vrc01_ref_airr.query("locus_light == 'IGL'")["cdr3_aa_light"].to_list()

    if is_vrc01_class:
        class_df_kappa = seq_df.query("locus_light == 'IGK'").query("is_vrc01_class").copy()
        class_df_lambda = seq_df.query("locus_light == 'IGL'").query("is_vrc01_class").copy()
    else:
        class_df_kappa = seq_df.query("locus_light == 'IGK'").query("is_vrc01_class==False").query("has_5_len").copy()
        class_df_lambda = seq_df.query("locus_light == 'IGL'").query("is_vrc01_class==False").query("has_5_len").copy()

    class_df_kappa["distance_to_known_lcdr3_kappa"] = class_df_kappa["cdr3_aa_light"].apply(
        lambda x: find_lowest_len(x, vrc01_ref_airr_kappa_seqs)
    )
    class_df_lambda["distance_to_known_lcdr3_lambda"] = class_df_lambda["cdr3_aa_light"].apply(
        lambda x: find_lowest_len(x, vrc01_ref_airr_lambda_seqs)
    )
    seq_df = pd.concat([class_df_kappa, class_df_lambda]).reset_index(drop=True)
    seq_df["distance_to_known_lcdr3"] = seq_df[["distance_to_known_lcdr3_kappa", "distance_to_known_lcdr3_lambda"]].min(
        axis=1
    )
    return (
        seq_df.groupby(["pseudogroup", "is_cp"])
        .apply(report_lcdr3)
        .stack()
        .reset_index()
        .rename(
            {0: "frequency"},
            axis=1,
        )
        .query("distance_to_known_lcdr3 == 0")
    )


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


def get_light_chain_usage_df(df: pd.DataFrame, pivot_on: str, is_cp: bool, normalize: bool = True) -> pd.DataFrame:
    return (
        df.query("is_vrc01_class==True")
        .query("is_cp==@is_cp")
        .groupby(pivot_on)
        .apply(lambda x: x["v_call_top_light_gene"].value_counts(normalize=normalize))
        .reset_index()
        .pivot(index=pivot_on, columns="level_1", values="v_call_top_light_gene")
    )


def spr_properties(data: Data) -> tuple[plt.Figure, list[pd.DataFrame]]:
    """Generate combined light chain and SPR properties plots."""
    plot = Plot()
    dfs = []
    apply_global_font_settings(10)
    palette = {False: "#91FCC0", True: "#6C65FF"}

    # Common variables
    light_genes = ["IGKV1-33", "IGKV3-20", "IGKV1-5", "IGKV3-15", "IGLV2-14", "IGLV2-23", "IGLV2-11", "Other"]
    psuedogroups = [
        r"eOD" + "\nwk8",
        r"eOD$\rightarrow$eOD" + "\nwk16",
        r"eOD$\rightarrow$core" + "\nwk 16",
        r"eOD$\rightarrow$eOD" + "\nwk 24",
        r"eOD$\rightarrow$core" + "\nwk 24",
        r"eOD$\rightarrow$eOD$\rightarrow$core" + "\nwk 24",
    ]

    # Data preparation
    g002_spr_boost_df = data.get_g002_spr_df_boost()
    g002_spr_boost_df = g002_spr_boost_df.query("top_c_call=='IGHG'").fillna(0)
    g002_spr_boost_df.loc[
        g002_spr_boost_df["v_call_top_light_gene"] == "IGKV3D-15", "v_call_top_light_gene"
    ] = "IGKV3-15"
    g002_spr_boost_df.loc[
        ~g002_spr_boost_df["v_call_top_light_gene"].isin(light_genes), "v_call_top_light_gene"
    ] = "Other"

    # Create figure with 5 subplots (2 for light chain, 3 for SPR)
    fig = plt.figure(figsize=(10, 12), dpi=500)
    gs = fig.add_gridspec(5, 3, width_ratios=[1, 1, 0.33], height_ratios=[2, 0.5, 1, 1, 1])
    axes = [
        fig.add_subplot(gs[0, 0]),  # Light chain random
        fig.add_subplot(gs[0, 1]),  # Light chain selected
        fig.add_subplot(gs[0, 2]),  # Reference data
        fig.add_subplot(gs[2, :]),  # LCDR3
        fig.add_subplot(gs[3, :]),  # QE at 96
        fig.add_subplot(gs[4, :]),  # W100b
    ]
    axes[0].letter = "A"
    axes[1].letter = "A"
    axes[2].letter = "A"
    axes[3].letter = "B"
    axes[4].letter = "C"
    axes[5].letter = "D"

    # Light chain color mapping
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
    colors = [explicit_mapping[i] for i in light_genes]

    # Plot light chain usage (random and selected)
    df = g002_spr_boost_df.copy(deep=True).drop_duplicates(["Ligand"])

    # Random sequences
    df_random = get_light_chain_usage_df(df, ["pseudogroup"], False)
    for c in light_genes:
        if c not in df_random.columns:
            df_random[c] = 0
    _df = df_random.loc[:, light_genes].copy(deep=True)
    plot.stacked_bar_plot_pivot_df(
        pivot_df=_df,
        title="Random",
        ylabel="% BCRs with VRC01-class\nbnAb $\mathregular{V_{K/L}}$",
        xticklabels=psuedogroups,
        xticklabel_rotation=45,
        yaxis_mtick_major=0.2,
        yaxis_mtick_minor=0.1,
        colors=colors,
        yaxis_mtick_normalize=False,
        fig=fig,
        ax=axes[0],
    )
    dfs.append(annotate_df(_df.reset_index(), axes[0], list(_df.columns), "random"))

    # Selected sequences
    df = g002_spr_boost_df.copy(deep=True).drop_duplicates(["Ligand"])
    _df = get_light_chain_usage_df(df, "pseudogroup", True).loc[:, light_genes].copy(deep=True)
    plot.stacked_bar_plot_pivot_df(
        pivot_df=_df,
        title="Selected",
        xticklabels=psuedogroups,
        remove_yticklabels=True,
        xticklabel_rotation=45,
        yaxis_mtick_major=0.2,
        yaxis_mtick_minor=0.1,
        colors=colors,
        fig=fig,
        ax=axes[1],
    )
    dfs.append(annotate_df(_df.reset_index(), axes[0], list(_df.columns), "selected"))

    # Reference data
    plot.plot_stacked_bar_dekosky_pivot_replicate_donors_by_light_value_counts_df(light_genes, ax=axes[2])

    # SPR properties
    g002_spr_boost_df["cdr3_aa_light"] = g002_spr_boost_df["cdr3_aa_light"].fillna("")
    g002_spr_boost_df = g002_spr_boost_df.query("Analyte=='core-Hx_r4.0D_TH6_g28v2_pHLsecAvi'")
    g002_spr_boost_df["has_QE_at_96"] = g002_spr_boost_df["cdr3_aa_light"].apply(
        lambda x: True if x[-2] in ["Q", "E"] else False
    )
    g002_spr_boost_df["has_w100b"] = g002_spr_boost_df["junction_aa_heavy"].str[-6] == "W"

    # Plot LCDR3, QE at 96, and W100b
    vrc01_ref_airr = data.get_vrc01_class_bnabs()
    df = bnab_lcdr3(g002_spr_boost_df.copy(deep=True), vrc01_ref_airr, is_vrc01_class=True)
    sns.barplot(
        data=df,
        y="frequency",
        x="pseudogroup",
        hue="is_cp",
        dodge=True,
        palette=palette,
        linewidth=1,
        edgecolor="black",
        ax=axes[3],
    )
    dfs.append(annotate_df(df, axes[3], ["frequency", "is_cp"], "bnab_lcdr3"))

    # QE at 96
    df = (
        g002_spr_boost_df.query("is_vrc01_class==True")
        .groupby(by=["pseudogroup", "is_cp"])
        .apply(determine_metric, "has_QE_at_96")
        .stack()
        .reset_index()
        .rename({"level_2": "has_QE_at_96", 0: "freq"}, axis=1)
        .query("has_QE_at_96==True")
    )
    sns.barplot(
        data=df,
        y="freq",
        x="pseudogroup",
        hue="is_cp",
        dodge=True,
        palette=palette,
        linewidth=1,
        edgecolor="black",
        ax=axes[4],
    )
    dfs.append(annotate_df(df, axes[4], ["freq", "is_cp"], "qe_at_96"))

    # W100b
    df = (
        g002_spr_boost_df.query("is_vrc01_class==True")
        .groupby(["pseudogroup", "is_cp"])
        .apply(determine_metric, "has_w100b")
        .stack()
        .reset_index()
        .rename({"level_2": "has_w100b", 0: "freq"}, axis=1)
        .query("has_w100b==True")
    )
    sns.barplot(
        data=df,
        y="freq",
        x="pseudogroup",
        hue="is_cp",
        dodge=True,
        palette=palette,
        linewidth=1,
        edgecolor="black",
        ax=axes[5],
    )
    dfs.append(annotate_df(df, axes[5], ["freq", "is_cp"], "w100b"))

    # Format SPR plots
    for ax in axes[3:]:
        ax.set_xlabel("")
        ax.set_ylim(0, 1)
        ax.get_legend().remove()
        ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0, symbol=None))

    axes[3].set_xticklabels([])
    axes[4].set_xticklabels([])
    axes[5].set_xticklabels(psuedogroups, rotation=45, ha="right", rotation_mode="anchor")

    # Set labels
    axes[3].set_ylabel("% BCRs with VRC01-class\nbnAb LCDR3")
    axes[4].set_ylabel("% BCRs with E/Q at\nLCDR3 position 96")
    axes[5].set_ylabel("% BCRs with\n" + r"$\mathregular{Trp_{103-5}}$ in HCDR3")

    # Add legends
    legend_pallete = {
        **{k[2:].replace("V", ""): explicit_mapping[k] for k in light_genes[:-1]},
        **{"Other": explicit_mapping["Other"]},
    }
    plot.side_legend([axes[0], axes[1], axes[2]], palette=legend_pallete, fontsize=10, bbox_to_anchor=(1.5, 1))
    plot.bottom_legend_random_selected(axes[5], distance=-0.75)

    sns.despine()
    # plt.tight_layout()

    for i in [3, 4, 5]:
        ax = axes[i]
        ax.text(-0.13, 1.1, ax.letter, transform=ax.transAxes, size=14, weight="bold")

    axes[0].text(-0.3, 1.1, axes[0].letter, transform=axes[0].transAxes, size=14, weight="bold")

    return fig, dfs
