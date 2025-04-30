from pprint import pprint as pp
from typing import Any, Dict, List, Literal, Tuple

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
import patchworklib as pw
import seaborn as sns
from Levenshtein import distance
from matplotlib import colormaps
from matplotlib.axes import Axes
from matplotlib.ticker import AutoMinorLocator

from g00x_figures.data import Data, calculate_resonse
from g00x_figures.plot_helpers.boxplot import adjust_boxplot, format_y_axis
from g00x_figures.plot_helpers.font import apply_global_font_settings
from g00x_figures.plot_helpers.legend import plot_legend
from g00x_figures.plots import Plot

data = Data()
palette = data.get_trial_g001_g002_g003_palette()
pw.overwrite_axisgrid()
apply_global_font_settings()


def get_g00x_spr(is_vrc01_class: bool = True) -> pd.DataFrame:
    g001_spr = data.get_g001_spr_df()
    g001_spr = g001_spr[~g001_spr["cdr3_aa_light"].isna()]
    g002_spr = data.get_g002_spr_df_prime()
    g003_spr = data.get_g003_spr_df_prime()
    df = pd.concat([g001_spr, g002_spr, g003_spr])
    df = df.query("weeks in [4, 8, 10, 16]")
    # df["hcdr3_len"] = df["hcdr3_len"].astype(int)
    df = df.query("is_vrc01_class==@is_vrc01_class")
    df = df.query('top_c_call == "IGHG"')
    df = df.query("is_cp == False")
    df = df.query("is_igl == False")
    return df


def get_g00x_prime(is_vrc01_class: bool = True) -> pd.DataFrame:
    df = data.get_g00x_sequences_prime().copy(deep=True)
    df = df.query("weeks in [4, 8, 10, 16]")
    df = df.query("is_vrc01_class==@is_vrc01_class")
    df = df.query('top_c_call == "IGHG"')
    # df = df.query("is_cp == False")
    # df = df.query("is_igl == False")
    return df


def get_g00x_all(is_vrc01_class: bool = True) -> pd.DataFrame:
    cols = [
        "pubID",
        "weeks",
        "trial",
        "is_spr",
        "hcdr3_len",
        "cdr3_aa_light",
        "v_mutation_aa_heavy",
        "v_mutation_aa_light",
        "cottrell_focused_v_common_score",
        "junction_aa_heavy",
        "locus_light",
    ]
    spr = get_g00x_spr(is_vrc01_class=is_vrc01_class)
    spr["is_spr"] = True
    study = get_g00x_prime(is_vrc01_class=is_vrc01_class)
    study["is_spr"] = False
    df = pd.concat([spr, study])[cols]
    df["has_QE_at_96"] = df["cdr3_aa_light"].apply(lambda x: True if x[-2] in ["Q", "E"] else False)
    df["has_w100b"] = df["junction_aa_heavy"].str[-6] == "W"
    df["hcdr3_len"] = df["hcdr3_len"].astype(int)
    df["cottrell_focused_v_common_score"] = df["cottrell_focused_v_common_score"].astype(np.float32)
    return df


def append_spr_fields_from_spr_df_column(df: pd.DataFrame, column: str = "Ligand") -> pd.DataFrame:
    """Pull needed "airr" fields from complex name give to SPR data"""
    # df[column] = df[column].apply(lambda x: x.split("_")[0])
    # map timepoint (str), garunteed at index 1, to weeks (int)
    df["weeks"] = (
        df.Ligand.str.split("_")
        .str[1]
        .map(
            {
                "V91": -5,
                "V101": 0,
                "V108": 1,
                "V115": 2,
                "V122": 3,
                "V201": 8,
                "V208": 9,
                "V215": 10,
                "V222": 11,
                "V229": 12,
                "V257": 16,
                "V292": 21,
            }
        )
    )
    df["is_igl"] = df["Ligand"].apply(lambda x: True if "_iGL_" in x else False)
    df["is_vrc01_class"] = df["Ligand"].apply(lambda x: True if "_VRC01" in x else False)
    df["trial"] = df.Ligand.str.split("-").str[0].str.split("_").str[0]
    df["ptid"] = df.Ligand.str.split("_").str[0]
    return df


def modify_spr_rows(
    df: pd.DataFrame,
    lt: float = 5e-5,
    trial: str = "G003",
    kd_lookup: Literal["KD_fix", "Kon", "Koff"] = "KD_fix",
) -> pd.DataFrame:
    if "trial" not in df.columns:
        df["trial"] = trial
    """Filter for rows with valid SPR data"""
    df["Analyte"] = df.Analyte.replace(
        {
            "core-Hx_r4.0D_TH6_g28v2_pHLsecAvi": "core-g28v2",
            "core-Hx_r4.0D_TH6_g28v2_KO11b_pHLsecAvi": "core-KO11b",
            "191084_SOSIP_MD39_N276D_mC": "191084-N276D",
            "eOD-GT8.1_His-Avi_mC": "eOD-GT8",
        }
    )
    df = df[df.Ligand.str.startswith(trial)]
    df.loc[df[df["estimated"] == True].index, kd_lookup] = lt
    df.loc[df[df[kd_lookup] > lt].index, kd_lookup] = lt
    return df


def filter_spr_rows(
    df: pd.DataFrame,
    is_igl: bool | None = False,
    is_vrc01_class: bool | None = True,
    is_cp: bool | None = False,
) -> pd.DataFrame:
    df = df[~df.Ligand.isna()]
    if "is_vrc01_class" in df.columns and is_vrc01_class is not None:
        df = df.query("is_vrc01_class == @is_vrc01_class")
    if "is_igl" in df.columns and is_igl is not None:
        df = df.query("is_igl == @is_igl")
    if "is_cp" in df.columns and is_cp is not None:
        df = df.query("is_cp == @is_cp")
    df = df.sort_values("Chi2").groupby(["Ligand", "Analyte"]).head(1)
    return df


def plot_violin(
    ax: pw.Brick,
    df: pd.DataFrame,
    x: Literal["weeks", "pseudogroup"],
    y: Literal[
        "hcdr3_len",
        "v_mutation_aa_heavy",
        "v_mutation_aa_light",
        "v_mutation_heavy",
        "v_mutation_light",
    ],
    outline: str,
    hue: Literal["trial", "weeks", "is_spr"],
    palette: dict,
    sharey: bool = True,
    sharex: bool = True,
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
        palette=palette,
        # dodge=dodge,
        inner="quartile",
        saturation=1,
        scale="area",
        hue_order=palette.keys(),
        # cut=0,
        gap=2,
        linewidth=2,
        edgecolor="black",
        split=True,
    )
    ax.legend().remove()
    if sharey:
        ax.set_ylabel(None)
    if sharex:
        ax.set_xlabel(None)
        ax.set_xticklabels([])
    sns.despine(ax=ax)
    if df[y].dtype == "float64":
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    for i, violin in enumerate(ax.collections):
        # Check if the current violin corresponds to the second hue level
        if i % 2 == 1:
            # Modify the linestyle
            # r, g, b = hex_to_rgb(palette[True])
            # violin.set_linestyle("--")
            violin.set_edgecolor(outline)
    # violin.set_edgecolor("black")
    # for line in violin.get_paths():
    #     line.set_facecolor("black")\
    return ax


def plot_bar(
    ax: pw.Brick,
    df: pd.DataFrame,
    x: Literal["weeks", "pseudogroup"],
    y: Literal[
        "hcdr3_len",
        "v_mutation_aa_heavy",
        "v_mutation_aa_light",
        "v_mutation_heavy",
        "v_mutation_light",
    ],
    hue: Literal["trial", "weeks", "is_spr"],
    palette: dict,
    sharey: bool = True,
    sharex: bool = True,
    *args,
    **kwargs,
) -> None:
    """Plot bar plot.

    >>> plot_bar(ax=ax df=g00x_seq, x=weeks, y=v_mutation_heavy, hue=trial)
    >>> plot_bar(ax=ax df=g002_seq, x=pseudogroup, y=v_mutation_heavy, hue=weeks)

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
    sns.barplot(
        ax=ax,
        data=df,
        x=x,
        y=y,
        hue=hue,
        hue_order=palette.keys(),
        palette=palette,
        ci=None,
    )
    ax.legend().remove()
    if sharey:
        ax.set_ylabel(None)
    if sharex:
        ax.set_xlabel(None)
        ax.set_xticklabels([])
    sns.despine(ax=ax)
    if df[y].dtype == "float64":
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    return ax


def _determine_metric(df: pd.DataFrame, metric: str) -> pd.Series:
    return pd.Series(
        {
            True: df[metric].sum() / len(df),
            False: 1 - (df[metric].sum() / len(df)),
        }
    )


def find_lowest_len(x, ref) -> int:
    return min(list(map(lambda y: distance(x, y), ref)))


def report_lcdr3(df: pd.DataFrame):
    b = df["distance_to_known_lcdr3"].value_counts().sort_index()
    c = (b / b.sum()).cumsum()
    c = c.reindex([float(i) for i in range(6)]).ffill().fillna(0.0).transpose()
    return c


def get_plottable() -> pd.DataFrame:
    vrc01_ref_airr = data.get_vrc01_class_bnabs()

    vrc01_ref_airr_kappa_seqs = vrc01_ref_airr.query("locus_light == 'IGK'")["cdr3_aa_light"].to_list()
    vrc01_ref_airr_lambda_seqs = vrc01_ref_airr.query("locus_light == 'IGL'")["cdr3_aa_light"].to_list()

    class_df_kappa = get_g00x_all().query("locus_light == 'IGK'")
    class_df_lambda = get_g00x_all().query("locus_light == 'IGL'")
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
        combined_df.groupby(["pubID", "weeks", "trial", "is_spr", "locus_light"])
        .apply(report_lcdr3)
        .stack()
        .reset_index()
        .rename(
            {0: "Kappa"},
            axis=1,
        )
        .query("distance_to_known_lcdr3 == 0")
    )
    plottable["x-axis"] = plottable["trial"] + plottable["locus_light"]
    return plottable


def plot_qc():
    dfs = []
    apply_global_font_settings(fontsize=14)

    palette = data.get_trial_g001_g002_g003_palette()
    plottable = get_plottable()

    def _determine_metric(df: pd.DataFrame, metric: str) -> pd.Series:
        return pd.Series(
            {
                True: df[metric].sum() / len(df),
                False: 1 - (df[metric].sum() / len(df)),
            }
        )

    b = (
        get_g00x_all()
        .groupby(
            [
                "pubID",
                "weeks",
                "trial",
                "is_spr",
            ]
        )
        .apply(_determine_metric, metric="has_QE_at_96")
        .stack()
        .reset_index()
        .rename({"level_4": "has_QE_at_96", 0: "QE_at_96_freq"}, axis=1)
        .query("has_QE_at_96==True")
    )
    d = (
        get_g00x_all()
        .groupby(
            [
                "pubID",
                "weeks",
                "trial",
                "is_spr",
            ]
        )
        .apply(_determine_metric, metric="has_w100b")
        .stack()
        .reset_index()
        .rename({"level_4": "has_w100b", 0: "has_w100b_freq"}, axis=1)
        .query("has_w100b==True")
    )
    d.head()

    def _determine_metric(df: pd.DataFrame, metric: str) -> pd.Series:
        return pd.Series(
            {
                True: df["has_QE_at_96"].sum() / len(df),
                False: 1 - (df["has_QE_at_96"].sum() / len(df)),
            }
        )

    g00x_all = get_g00x_all(is_vrc01_class=True)

    figsize = (4, 2)
    ylabel_map = {
        "cottrell_focused_v_common_score": "number of VRC01-class\nHC residues",
        "QE_at_96_freq": "% BCRs with E/Q at\nLCDR3 position 96",
        "Kappa": "% BCRs with VRC01-class\nbnAb LCDR3",
        "has_w100b_freq": "% BCRs with\n" + r"$\mathregular{Trp_{103-5}}$ in HCDR3",
        "hcdr3_len": "HCDR3 Length",
        "v_mutation_aa_heavy": r"$\mathregular{V_H}$" + f" gene\n%mutation (aa)",
        "v_mutation_aa_light": r"$\mathregular{V_{K/L}}$" + f" gene\n%mutation (aa)",
    }
    darker_palette = {
        "G001": "#0800B2",
        "G002": "#05C157",
        "G003": "#8F1D6C",
    }
    trials = ["G001", "G002", "G003"]
    # the number of VRC01-class HC residues
    # Willis et al. submitted fig. S31 B
    # Willis et al. submitted fig. S31 C
    # Willis et al. submitted fig. S31 D

    # palette = {True: palette[trial], False: "gold"}

    y_df = [
        (
            g00x_all,
            "hcdr3_len",
            plot_violin,
        ),
        (g00x_all, "v_mutation_aa_heavy", plot_violin),
        (g00x_all, "v_mutation_aa_light", plot_violin),
        (g00x_all, "cottrell_focused_v_common_score", plot_violin),
        (b, "QE_at_96_freq", plot_bar),
        (
            plottable[plottable["locus_light"] == "IGK"],
            "Kappa",
            plot_bar,
        ),
        (d, "has_w100b_freq", plot_bar),
    ]
    cols = []
    for (j, trial), letter in zip(enumerate(trials), "ABCDEFG"):
        axes = []
        for i, (df, y, func) in enumerate(
            y_df,
        ):
            ax = pw.Brick(figsize=figsize)
            ax.letter = letter
            _df = df.query(f"trial == '{trial}'").copy(deep=True)
            g = func(
                ax=ax,
                df=_df,
                x="weeks",
                y=y,
                hue="is_spr",
                palette={True: palette[trial], False: "#bf9b30"},
                outline=darker_palette[trial],
                sharex=False if i == len(y_df) - 1 else True,
                sharey=False if j == 0 else True,
            )
            _df.name = f"{trial}_{y}"
            _df.key = [y, "is_spr"]
            dfs.append(_df)
            if j == 0:
                ax.set_ylabel(ylabel_map.get(y, y))
            if y == "hcdr3_len":
                ax.set_ylim(0, 30)
            elif y == "v_mutation_aa_heavy":
                ax.set_ylim(0, 0.2)
            elif y == "v_mutation_aa_light":
                ax.set_ylim(0, 0.2)
            elif y == "cottrell_focused_v_common_score":
                ax.set_ylim(-1, 6)
            else:
                ax.set_ylim(0, 1)
            ax.set_title(" ", fontsize=16, fontweight="bold", x=-0.2, y=1.1)
            axes.append(g)
        handles, labels = axes[-1].get_legend_handles_labels()
        # if isinstance(func, type(plot_violin)):
        #     for handle in handles[1:]:
        #         handle.set_linewidth(2)
        #         handle.set_edgecolor(darker_palette[trial])
        #         handle.set_facecolor("#bf9b30")
        axes[-1].legend(
            handles,
            ["SPR Sampled", "Study Population"],
            # title="",
            loc="lower center",
            bbox_to_anchor=(0.5, -0.65),
            frameon=False,
        )
        if j == 0:
            axes[0].set_title("A", fontsize=16, fontweight="bold", x=-0.2, y=1.1)
            axes[1].set_title("B", fontsize=16, fontweight="bold", x=-0.2, y=1.1)
            axes[2].set_title("C", fontsize=16, fontweight="bold", x=-0.2, y=1.1)
            axes[3].set_title("D", fontsize=16, fontweight="bold", x=-0.2, y=1.1)
            axes[4].set_title("E", fontsize=16, fontweight="bold", x=-0.2, y=1.1)
            axes[5].set_title("F", fontsize=16, fontweight="bold", x=-0.2, y=1.1)
            axes[6].set_title("G", fontsize=16, fontweight="bold", x=-0.2, y=1.1)

        c = pw.stack(axes, margin=0.15, operator="/")
        c.set_suptitle(f"{trial} VRC01", fontsize=16, y=0.98)
        cols.append(c)

    g = pw.stack(cols, margin=-0.1, operator="|")
    return g, dfs
    # TODO:
    # - A: G001 SPR
    #
    g.savefig("Sup/Leggat-S29-plus.png", dpi=500)
