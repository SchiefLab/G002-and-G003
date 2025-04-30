from pathlib import Path
from typing import Any

import matplotlib.ticker as mtick
import numpy as np
import pandas as pd
import patchworklib as pw
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

from g00x_figures.box_and_scatter.flow_frequencies import (
    MinorSymLogLocator,
    adjust_boxplot,
)
from g00x_figures.data import Data
from g00x_figures.plot_helpers.font import apply_global_font_settings
from g00x_figures.plot_helpers.legend import plot_legend
from g00x_figures.plots import Plot

lt_mapping = {
    "191084_SOSIP_MD39_N276D_mC": 5e-5,
    "eOD-GT8.1_His-Avi_mC": 1e-4,
    "eOD-GT8.1_KO11_pHLsecAvi": 1e-4,
    "core-Hx_r4.0D_TH6_g28v2_pHLsecAvi": 1e-4,
    "core-Hx_r4.0D_TH6_g28v2_KO11b_pHLsecAvi": 1e-4,
    "core-Hx_r4.0_TH6_g28v2_N276_pHLsecAvi": 1e-4,
}

fixed_locator_lookup = {
    1e-13: "100 fM",
    1e-12: "1 pM",
    1e-11: "10 pM",
    1e-10: "100 pM",
    1e-09: "1 nM",
    1e-08: "10 nM",
    1e-07: "100 nM",
    1e-06: "1 $\mathregular{\mu}$M",  # noqa: W605
    1e-05: "10 $\mathregular{\mu}$M",  # noqa: W605
    0.0001: "100$\mathregular{\mu}$M",  # noqa: W605
    # 0.0001: "$\mathregular{\geq}$100$\mathregular{\mu}$M",  # noqa: W605
    0.001: "1 mM",
    0.01: "10 mM",
    0.1: "100 mM",  # 1e1
    0: "1 M",  # 1
    1: "10 M",  # 10
    2: "100 M",
    10: "1 G",
}


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


def plot_spr_core_candidates(
    data: Data,
    outpath: Path,
    metric_outdir: Path,
    lim: float = 5e-5,
    is_vrc01_class: bool = True,
    is_cp: bool = True,
    ylabel: str = "$\mathregular{K_D}$",
    analyte_mapping: dict[str, str] | None = None,
    use_geomean: bool = True,
    median_scale: float = 1,
) -> None:
    # apply_global_font_settings()
    # Analyte mapping is a dictionary of the form {analyte: publication_name}
    # Can only pick 4 analytes at a time
    analytes = analyte_mapping or {
        "191084_SOSIP_MD39_N276D_mC": "191084-N276D",
        "1HD2_B4_S62_MD39v2_L14_N276Q_mC2": "1HD2-N276Q",
        "001428_MD39_L14_T278M_m2": "001428-T278M",
        "BG505_SOSIP_MD39_N276Q_mC": "BG505-N276Q",
        "235_47_RnS_2G_L14_T278M_mC2": "235-T278M",
        "BG505_MD39.3_cd4bsHxB2_M278_mC2": "BG505-cd4bsHxB2-T278M",
        # "001428_MD39_L14_m2": "001428",
        # "BG505_SOSIP_MD39_m": "BG505",
        # "235_47_RnS_2G_L14_T278M_m2": "235-T278M",
        # "core-Hx_r4.0D_TH6_g28v2_pHLsecAvi": "core-g28v2",
        # "core-Hx_r4.0_TH6_g28v2_N276_pHLsecAvi": "core-N276",
        # "BG505_MD39.3_cd4bsHxB2_M278_mC2": "BG505-HXB2CD4bs-T278M",
        # "BG505_MD39.3_cd4bsHxB2_mC2": "BG505-HXB2CD4bs",
    }
    # fig, axes = plt.subplots(
    #     6,
    #     2,
    #     figsize=(8.5, 11),
    #     height_ratios=[0.5, 2, 0.5, 2, 0.5, 2],
    #     sharex=False,
    #     sharey="row",
    #     gridspec_kw={"hspace": 1},
    # )
    # axes = axes.flatten()
    tablesize = (4, 1)
    figsize = (4, 2)
    axes = [
        pw.Brick(figsize=tablesize),
        pw.Brick(figsize=tablesize),
        pw.Brick(figsize=figsize),
        pw.Brick(figsize=figsize),
        pw.Brick(figsize=tablesize),
        pw.Brick(figsize=tablesize),
        pw.Brick(figsize=figsize),
        pw.Brick(figsize=figsize),
        pw.Brick(figsize=tablesize),
        pw.Brick(figsize=tablesize),
        pw.Brick(figsize=figsize),
        pw.Brick(figsize=figsize),
    ]

    for i, analyte, letter in zip([0, 1, 4, 5, 8, 9], analytes, ["A", "B", "C", "D", "E", "F"]):
        ax = axes[i]
        analyte_pub_name = analytes[analyte]
        plot_df = data.get_g002_spr_df_boost().query("Analyte==@analyte").query(f"is_cp=={is_cp}")
        plot_df = plot_df.query(f"is_vrc01_class=={is_vrc01_class}")
        plot_df.loc[plot_df[plot_df["estimated"] == True].index, "KD_fix"] = lim
        plot_df.loc[plot_df[plot_df["KD_fix"] > lim].index, "KD_fix"] = lim
        plot_df = plot_df.sort_values("Chi2").groupby(["Ligand", "Analyte"]).head(1)
        if i in [0, 4, 8]:
            remove_label = False
        else:
            remove_label = True
        try:
            table_ax, table_df = Plot().plot_dyn_table(
                plot_df,
                groupby=["pseudogroup"],
                ax=ax,
                remove_label=remove_label,
                lim=lim,
                loc=None,
                median_scale=median_scale,
                geomean=use_geomean,
                # loc="upper",
            )
            table_df = table_df.T.reset_index()
            table_df = data.populate_psname_spr(table_df)
            table_df.to_csv(
                metric_outdir / f"fig{letter}-table-spr-{analyte}.csv",
                index=False,
            )
        except:
            print(f"failed to plot table {analyte}")
    for i, analyte, letter in zip([2, 3, 6, 7, 10, 11], analytes, ["A", "B", "C", "D", "E", "F"]):
        ax = axes[i]
        analyte_pub_name = analytes[analyte]
        plot_df = data.get_g002_spr_df_boost().query("Analyte==@analyte").query(f"is_cp=={is_cp}")
        plot_df = plot_df.query("is_vrc01_class==@is_vrc01_class")
        plot_df.loc[plot_df[plot_df["estimated"] == True].index, "KD_fix"] = lim
        plot_df.loc[plot_df[plot_df["KD_fix"] > lim].index, "KD_fix"] = lim
        plot_df = plot_df.sort_values("Chi2").groupby(["Ligand", "Analyte"]).head(1)

        data_for_strip = plot_df.copy(deep=True)
        data_for_box = plot_df.copy(deep=True)

        # if plot_lt:
        #     # any data less than
        #     di = plot_df[plot_df["KD_fix"] < lt].index

        #     # any data less than, just make it that lt number
        #     data_for_strip.loc[di, ""] = lt

        #     # but for box, just remove it
        #     data_for_box.loc[di, g002_label] = np.nan

        sns.stripplot(
            x="pseudogroup",
            y="KD_fix",
            hue="weeks",
            data=plot_df.query("is_cp==@is_cp"),
            dodge=False,
            edgecolor="black",
            linewidth=1,
            # jitter=0.5,
            order=range(1, 7),
            palette=data.get_week_palette(),
            ax=ax,
        )
        sns.boxplot(
            x="pseudogroup",
            y="KD_fix",
            hue="weeks",
            data=plot_df.query("is_cp==@is_cp"),
            dodge=False,
            fliersize=0,
            whis=[10, 90],
            order=range(1, 7),
            palette=data.get_week_palette(),
            ax=ax,
        )
        # data.populate_psname_spr(
        plot_df[
            [
                "pubID",
                "pseudogroup",
                "is_cp",
                "weeks",
                "KD_fix",
                "Kon",
                "Koff",
                "Ligand",
                "Analyte",
            ]
        ].copy(deep=True).to_csv(
            metric_outdir / f"fig{letter}-spr-{analyte}.csv",
            index=False,
        )
        ax.set_xlabel("")
        ax.set_ylim(1e-10, 0.00015)
        ax.set(yscale="log")
        ax.set_ylabel("")
        ax.set_xlabel("")
        ax.set_xticklabels([])
        ax.set_title(analyte_pub_name)
        ax.legend().remove()
        ax.axhline(y=lim, linestyle="--", linewidth=1, color="black", zorder=-10)
        adjust_boxplot(ax)

    axes[2].set_ylabel(ylabel, size=14)
    axes[6].set_ylabel(ylabel, size=14)
    axes[10].set_ylabel(ylabel, size=14)
    for i in [2, 3, 6, 7, 10, 11]:
        # axes[4].set_ylabel(ylabel, size=14)
        axes[i].yaxis.set_major_locator(mtick.LogLocator(base=10, numticks=9))
        # custom locator, shouldnt be linear anywhere
        axes[i].yaxis.set_minor_locator(MinorSymLogLocator(linthresh=1))
        axes[i].yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: fixed_locator_lookup[x]))
    fontsize = 9
    for i in [2, 3, 6, 7, 10, 11]:
        axes[i].set_xticklabels(
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
            fontsize=fontsize,
        )

    pal = {
        8: "gold",
        16: "#E377C2",
        24: "#2078B4",
    }
    custom_lines = []
    for x in pal:
        custom_lines.append(Patch(facecolor=pal[x], edgecolor="black", linewidth=1, label=x))

    for i in [6, 7]:
        axes[i].legend(
            custom_lines,
            ["wk " + i._label for i in custom_lines],
            loc="upper center",
            frameon=False,
            handlelength=0.8,
            ncol=3,
            bbox_to_anchor=(0.5, -0.5),
            labelspacing=0.1,
            fontsize=12,
        )
    for ax in axes:
        sns.despine(ax=ax)
    for i, letter in zip([0, 1, 4, 5], "ABCD"):
        if i % 2 == 0:
            x = -0.2
        else:
            x = -0.1
        axes[i].set_title(letter, fontsize=18, fontweight="bold", x=x, y=0.9)
    # axes[5].set_axis_off()
    plt.tight_layout()
    # c1 = pw.stack(
    #     [axes[0], axes[2], axes[4], axes[6], axes[8], axes[10]],
    #     margin=0,
    #     operator="/",
    # )
    # c2 = pw.stack(
    #     [axes[1], axes[3], axes[5], axes[7], axes[9], axes[11]],
    #     margin=0,
    #     operator="/",
    # )
    c1 = pw.stack(
        [axes[0], axes[2], axes[4], axes[6]],
        margin=0,
        operator="/",
    )
    c2 = pw.stack(
        [axes[1], axes[3], axes[5], axes[7]],
        margin=0,
        operator="/",
    )
    g = pw.stack([c1, c2], margin=0.25, operator="|")
    g.savefig(outpath, dpi=300)
    # plt.savefig(outpath, dpi=300)


def plot_spr_prime(
    data: Data,
    analyte: str,
    is_vrc01_class: bool = True,
    letter: str = "A",
    igl_pad: int = 7,
    remove_igl: bool = False,
    mature_pad: int = 7,
    use_geomean: bool = True,
    median_scale: float = 1,
    sci_notation: bool = True,
) -> tuple[pw.Bricks, list[pd.DataFrame]]:
    # fig23 in paper
    dfs = []
    apply_global_font_settings(8)
    trial_palette = data.get_trial_g001_g002_g003_palette()
    tablesize = (4, 1)
    figsize = (4, 2)
    axes = [
        pw.Brick(figsize=(1, 1)),
        pw.Brick(figsize=tablesize),
        pw.Brick(figsize=(1, 2)),
        pw.Brick(figsize=figsize),
    ]

    g001_spr = data.get_g001_spr_df()
    g001_spr["trial"] = "G001"
    g002_spr = data.get_g002_spr_df_prime()
    g002_spr["trial"] = "G002"
    g003_spr = data.get_g003_spr_df_prime()
    g003_spr["trial"] = "G003"

    lt = lt_mapping[analyte]

    # estimated existing only in G002/3
    g002_spr.loc[g002_spr[g002_spr["estimated"] == True].index, "KD_fix"] = lt
    g002_spr.loc[
        g002_spr[(g002_spr["estimated"] == True) & ((~g002_spr["KD_fix_iGL"].isna()))].index,
        "KD_fix_iGL",
    ] = lt

    g003_spr.loc[g003_spr[g003_spr["estimated"] == True].index, "KD_fix"] = lt
    g003_spr.loc[
        g003_spr[(g003_spr["estimated"] == True) & ((~g003_spr["KD_fix_iGL"].isna()))].index,
        "KD_fix_iGL",
    ] = lt

    def get_filtered_df(g001_spr, g002_spr, g003_spr):
        combined_df = pd.concat([g001_spr, g002_spr, g003_spr])
        combined_df = combined_df.query("weeks in [4, 8, 10, 16, 21]")
        combined_df.loc[combined_df[combined_df["KD_fix"] > lt].index, "KD_fix"] = lt
        combined_df.loc[combined_df[combined_df["KD_fix_iGL"] > lt].index, "KD_fix_iGL"] = lt
        combined_df = combined_df.query(f"is_vrc01_class=={is_vrc01_class}").query("is_cp==False")
        combined_df = combined_df.query(f"Analyte=='{analyte}'")
        combined_df = combined_df.sort_values("Chi2").groupby(["Ligand", "Analyte", "trial"]).head(1)
        combined_df = combined_df.drop_duplicates(["Ligand"])
        return combined_df

    mature_df = get_filtered_df(
        g001_spr.query("is_igl==False"),
        g002_spr,
        g003_spr,
    )
    igl_df = get_filtered_df(g001_spr.query("is_igl"), g002_spr, g003_spr)
    igl_df["weeks"] = "iGL"
    if remove_igl:
        axes[1].set_title(letter, fontsize=12, fontweight="bold", x=-0.2, y=0.8)
    else:
        axes[0].set_title(letter, fontsize=12, fontweight="bold", x=-0.8, y=0.8)

    # mature_df = mature_df.query("top_c_call=='IGHG'")
    igl_ax = axes[2]
    mature_ax = axes[3]

    mature_df.name = letter + "_mature_" + analyte
    igl_df.name = letter + "_iGL_" + analyte

    mature_df.key = ["KD_fix", "cellid", "Ligand"]
    igl_df.key = ["KD_fix_iGL"]

    dfs.append(mature_df)
    dfs.append(igl_df)

    # from IPython import embed

    # embed()

    # igl_df = combined_df.copy(deep=True)
    # igl_df["weeks"] = "iGL"
    table_ax, table_df = Plot().plot_dyn_table(
        igl_df,
        groupby=["weeks", "trial"],
        ax=axes[0],
        lim=lt,
        geomean=use_geomean,
        median_scale=median_scale,
        sci_notation=sci_notation,
        KD_fix="KD_fix_iGL",
        fontsize=12,
        fill_na=True,
        loc=None,
        pad=igl_pad,
    )
    table_df = table_df.T.reset_index()
    # Convert empty strings to NaN
    table_df.replace("-", np.nan, inplace=True)
    table_df.replace(" ", np.nan, inplace=True)
    # Convert numeric columns to float
    table_df = table_df.apply(pd.to_numeric, errors="ignore")
    table_df.name = letter + "_header_iGL_" + analyte
    table_df.key = list(table_df.columns)
    dfs.append(table_df)

    sns.stripplot(
        x="weeks",
        y="KD_fix_iGL",
        hue="trial",
        data=igl_df,
        ax=igl_ax,
        dodge=True,
        palette=trial_palette,
        edgecolor="black",
        hue_order=["G001", "G002", "G003"],
        linewidth=1,
    )

    sns.boxplot(
        x="weeks",
        y="KD_fix_iGL",
        hue="trial",
        data=igl_df,
        ax=igl_ax,
        dodge=True,
        fliersize=0,
        whis=[10, 90],
        hue_order=["G001", "G002", "G003"],
        palette=trial_palette,
    )

    igl_ax.set_yscale("log")
    igl_ax.legend_.set_visible(False)
    antigen = analyte
    if antigen == "core-Hx_r4.0D_TH6_g28v2_pHLsecAvi":
        antigen = "core-g28v2\n"
    if antigen == "eOD-GT8.1_His-Avi_mC":
        antigen = "eOD-GT8\n"
    if is_vrc01_class:
        vrc01 = "VRC01-class"
    else:
        vrc01 = "Non-VRC01-class"
    igl_ax.set_ylabel("$\mathregular{K_D}$\n" + antigen + vrc01)
    igl_ax.set_xlabel("")
    igl_ax.axhline(y=lt, linestyle="--", linewidth=1, color="black", zorder=-10)
    adjust_boxplot(igl_ax)

    # from IPython import embed

    # embed()

    table_ax, table_df = Plot().plot_dyn_table(
        mature_df,
        groupby=["weeks", "trial"],
        ax=axes[1],
        remove_label=False if remove_igl else True,
        geomean=use_geomean,
        median_scale=median_scale,
        sci_notation=sci_notation,
        lim=lt,
        fill_na=True,
        fontsize=12,
        loc=None,
        pad=mature_pad,
    )
    table_df = table_df.T.reset_index()
    # Convert empty strings to NaN
    table_df.replace("-", np.nan, inplace=True)
    table_df.replace(" ", np.nan, inplace=True)
    # Convert numeric columns to float
    table_df = table_df.apply(pd.to_numeric, errors="ignore")

    table_df.name = letter + "_header_mature_" + analyte
    table_df.key = list(table_df.columns)
    dfs.append(table_df)

    # ax.set_label("")
    sns.stripplot(
        x="weeks",
        y="KD_fix",
        hue="trial",
        data=mature_df,
        ax=mature_ax,
        dodge=True,
        palette=trial_palette,
        edgecolor="black",
        hue_order=["G001", "G002", "G003"],
        linewidth=1,
    )
    sns.boxplot(
        x="weeks",
        y="KD_fix",
        hue="trial",
        data=mature_df,
        ax=mature_ax,
        dodge=True,
        fliersize=0,
        whis=[10, 90],
        hue_order=["G001", "G002", "G003"],
        palette=trial_palette,
    )
    mature_ax.set_yscale("log")
    mature_ax.set_xlabel("Weeks Post Vaccine")
    if remove_igl:
        mature_ax.set_ylabel("$\mathregular{K_D}$\n" + antigen + vrc01)
    else:
        mature_ax.set_ylabel("")
        mature_ax.tick_params(axis="y", labelleft=False)
    # mature_ax.set_xticklabels([4, 8, 16])
    mature_ax.legend_.set_visible(False)
    mature_ax.axhline(y=lt, linestyle="--", linewidth=1, color="black", zorder=-10)

    mature_ax.set_ylim(1e-12, 0.00015)
    igl_ax.set_ylim(1e-12, 0.00015)

    adjust_boxplot(mature_ax)
    plot_legend(
        ax=mature_ax,
        pallete=trial_palette,
        fontsize=8,
        bbox_to_anchor=(0.5, -0.35),
    )
    for ax in [igl_ax, mature_ax]:
        ax.yaxis.set_major_locator(mtick.LogLocator(base=10, numticks=9))
        # custom locator, shouldnt be linear anywhere
        ax.yaxis.set_minor_locator(MinorSymLogLocator(linthresh=1))
        ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: fixed_locator_lookup[x]))
    for ax in axes:
        sns.despine(ax=ax)
    if remove_igl:
        c1 = pw.stack(
            [axes[1], pw.spacer(axes[1], 0.025)],
            margin=0.2,
            operator="|",
        )
        g = pw.stack(
            [c1, axes[3]],
            margin=0.1,
            operator="/",
        )
        return g, dfs
    c1 = pw.stack(
        [axes[0], pw.spacer(axes[1], 0.02), axes[1], pw.spacer(axes[1], 0.03)],
        margin=0.2,
        operator="|",
    )
    c2 = pw.stack(
        [axes[2], axes[3]],
        margin=0.1,
        operator="|",
    )
    g = pw.stack([c1, c2], margin=0, operator="/")
    return g, dfs


def plot_core_spr(
    data: Data,
    outpath: Path,
    metric_outdir: Path,
    is_vrc01_class: bool = True,
    kd_lookup: str = "KD_fix",
    use_geomean: bool = True,
    median_scale: float = 1,
    sci_notation: bool = False,
) -> None:
    lt = lt_mapping["core-Hx_r4.0D_TH6_g28v2_pHLsecAvi"]
    spr_df = data.get_g002_spr_df_boost()
    spr_df["ptid"] = spr_df["pubID"]
    plot = Plot()
    spr_df = spr_df.query(f"is_vrc01_class=={is_vrc01_class}").query("top_c_call=='IGHG'")
    spr_df.loc[spr_df[spr_df["estimated"] == True].index, kd_lookup] = lt
    spr_df.loc[spr_df[spr_df[kd_lookup] > lt].index, kd_lookup] = lt
    spr_df = spr_df.sort_values("Chi2").groupby(["Ligand", "Analyte"]).head(1)
    spr_df = spr_df.query("Analyte=='core-Hx_r4.0D_TH6_g28v2_pHLsecAvi'")
    fig, axes = plt.subplots(
        4,
        1,
        figsize=(9.5, 10),
        gridspec_kw={
            # # "width_ratios": [0.5],
            "height_ratios": [0.05, 1, 0.05, 1],
            "hspace": 0.9,
            "wspace": 0.1,
            # "bottom": 0.,
            "top": 0.85,
            "left": 0.18,
            "right": 0.9,
        },
    )
    table_ax, table_df = plot.plot_dyn_table(
        spr_df,
        ax=axes[0],
        groupby=["pseudogroup", "is_cp"],
        lim=lt,
        KD_fix=kd_lookup,
        geomean=use_geomean,
        median_scale=median_scale,
        sci_notation=sci_notation,
    )
    # table_df = table_df.reset_index(drop=False)
    table_df = table_df.T.reset_index()
    table_df = data.populate_psname_spr(table_df)
    # table_df = table_df.sort_values(["is_cp", "pseudogroup", "KD_fix"])
    table_df.to_csv(metric_outdir / "figA-table-spr-g28v2.csv", index=False)
    # plot.plot_dyn_table(analyte="core-Hx_r4.0_TH6_g28v2_N276_pHLsecAvi", ax=axes[2], is_vrc01_class=is_vrc01_class)

    pal = {True: "#6C65FF", False: "#91FCC0"}

    metric_df = spr_df[
        [
            "ptid",
            "pseudogroup",
            "is_cp",
            "weeks",
            "KD_fix",
            # "Kon",
            # "Koff",
            "Ligand",
            "Analyte",
            # "mutations_heavy_and_light_count",
            # "cottrell_focused_v_common_score",
        ]
    ]
    metric_df = data.populate_psname_spr(metric_df)
    metric_df.to_csv(metric_outdir / f"figA-spr-g28v2.csv", index=False)

    sns.stripplot(
        data=spr_df,
        x="pseudogroup",
        y=kd_lookup,
        dodge=True,
        hue="is_cp",
        edgecolor="black",
        linewidth=1,
        palette=pal,
        ax=axes[1],
    )
    ax = sns.boxplot(
        data=spr_df,
        x="pseudogroup",
        y=kd_lookup,
        dodge=True,
        hue="is_cp",
        fliersize=0,
        palette=pal,
        whis=[10, 90],
        ax=axes[1],
    )
    ax.set_yscale("log")
    ax.legend().set_visible(False)
    ax.yaxis.set_major_locator(mtick.LogLocator(base=10, numticks=9))
    # custom locator, shouldnt be linear anywhere
    ax.yaxis.set_minor_locator(MinorSymLogLocator(linthresh=1))
    ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: fixed_locator_lookup[x]))
    ax.set_ylabel("$\mathregular{K_D}$ for core-g28v2", size=14)
    ax.set_xlabel("")
    # ax[0].set_title("Core")
    adjust_boxplot(ax)
    ax.set_xticklabels([])
    ax.set_xticklabels(
        [
            r"eOD" + "\nwk8",
            r"eOD$\rightarrow$eOD" + "\nwk16",
            r"eOD$\rightarrow$core" + "\nwk 16",
            r"eOD$\rightarrow$eOD" + "\nwk 24",
            r"eOD$\rightarrow$core" + "\nwk 24",
            r"eOD$\rightarrow$eOD$\rightarrow$core" + "\nwk 24",
        ],
        rotation=0,
        verticalalignment="top",
        horizontalalignment="center",
    )
    ax.axhline(y=lt, linestyle="--", linewidth=1, color="black", zorder=-10)
    custom_lines = []
    for x in {
        "Random": "#91FCC0",
        "Selected": "#6C65FF",
    }:
        custom_lines.append(
            Patch(
                facecolor={"Random": "#91FCC0", "Selected": "#6C65FF"}[x],
                edgecolor="black",
                linewidth=1,
                label=x,
            )
        )
    ax.legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=3,
        bbox_to_anchor=(0.5, -0.2),
        labelspacing=0.1,
    )
    ax.set_ylim(1e-11, 0.00015)

    lt = lt_mapping["core-Hx_r4.0_TH6_g28v2_N276_pHLsecAvi"]
    spr_df = data.get_g002_spr_df_boost()
    spr_df["ptid"] = spr_df["pubID"]
    spr_df = spr_df.query(f"is_vrc01_class=={is_vrc01_class}").query("top_c_call=='IGHG'")
    spr_df = spr_df.query("Analyte=='core-Hx_r4.0_TH6_g28v2_N276_pHLsecAvi'")
    spr_df.loc[spr_df[spr_df["estimated"] == True].index, kd_lookup] = lt
    spr_df.loc[spr_df[spr_df[kd_lookup] > lt].index, kd_lookup] = lt
    spr_df = spr_df.sort_values("Chi2").groupby(["Ligand", "Analyte"]).head(1)
    table_ax, table_df = plot.plot_dyn_table(
        spr_df,
        ax=axes[2],
        groupby=["pseudogroup", "is_cp"],
        lim=lt,
        geomean=use_geomean,
        sci_notation=sci_notation,
        median_scale=median_scale,
    )
    table_df = table_df.T.reset_index()
    table_df = data.populate_psname_spr(table_df)
    table_df.to_csv(metric_outdir / "figB-table-spr-N276.csv", index=False)

    metric_df = spr_df[
        [
            "ptid",
            "pseudogroup",
            "is_cp",
            "weeks",
            "KD_fix",
            # "Kon",
            # "Koff",
            "Ligand",
            "Analyte",
            # "mutations_heavy_and_light_count",
            # "cottrell_focused_v_common_score",
        ]
    ]
    metric_df = data.populate_psname_spr(metric_df)
    metric_df.to_csv(metric_outdir / "figB-spr-N276.csv", index=False)
    sns.stripplot(
        data=spr_df,
        x="pseudogroup",
        y=kd_lookup,
        dodge=True,
        hue="is_cp",
        edgecolor="black",
        linewidth=1,
        palette=pal,
        ax=axes[3],
    )
    ax = sns.boxplot(
        data=spr_df,
        x="pseudogroup",
        y=kd_lookup,
        dodge=True,
        hue="is_cp",
        fliersize=0,
        palette=pal,
        whis=[10, 90],
        ax=axes[3],
    )
    ax.set_yscale("log")
    ax.set_ylim(1e-11, 0.00015)
    ax.legend().set_visible(False)
    ax.yaxis.set_major_locator(mtick.LogLocator(base=10, numticks=9))
    # custom locator, shouldnt be linear anywhere
    ax.yaxis.set_minor_locator(MinorSymLogLocator(linthresh=1))
    ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: fixed_locator_lookup[x]))
    ax.set_ylabel("$\mathregular{K_D}$ for core-N276", size=14)
    ax.set_xlabel("")
    adjust_boxplot(ax)
    ax.set_xticklabels(
        [
            r"eOD" + "\nwk8",
            r"eOD$\rightarrow$eOD" + "\nwk16",
            r"eOD$\rightarrow$core" + "\nwk 16",
            r"eOD$\rightarrow$eOD" + "\nwk 24",
            r"eOD$\rightarrow$core" + "\nwk 24",
            r"eOD$\rightarrow$eOD$\rightarrow$core" + "\nwk 24",
        ],
        rotation=0,
        verticalalignment="top",
        horizontalalignment="center",
    )

    custom_lines = []
    for x in {
        "Random": "#91FCC0",
        "Selected": "#6C65FF",
    }:
        custom_lines.append(
            Patch(
                facecolor={"Random": "#91FCC0", "Selected": "#6C65FF"}[x],
                edgecolor="black",
                linewidth=1,
                label=x,
            )
        )
    ax.legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=3,
        bbox_to_anchor=(0.5, -0.2),
        labelspacing=0.1,
    )
    ax.axhline(y=lt, linestyle="--", linewidth=1, color="black", zorder=-10)
    sns.despine()
    plt.tight_layout()
    fig.savefig(outpath, dpi=300)


def plot_core_spr_kon_off(
    data: Data,
    analyte_long_name: str,
    analyte_short_name: str,
    lt: float = 5e-5,
) -> tuple[plt.Figure, list[pd.DataFrame]]:
    lt = lt_mapping[analyte_long_name]
    dfs = []
    # SPR BOOST Load
    spr_df = data.get_g002_spr_df_boost()
    # Plot utility
    plot = Plot()
    pal = {True: "#6C65FF", False: "#91FCC0"}

    # Filters
    spr_df = spr_df.query(f"is_vrc01_class")
    spr_df = spr_df.query("top_c_call=='IGHG'")
    spr_df = spr_df.query(f"Analyte=='{analyte_long_name}'")

    # First filter for valid KD_fix then use Koff or Kon
    kd_lookup = "KD_fix"
    spr_df.loc[spr_df[spr_df["estimated"] == True].index, kd_lookup] = np.nan
    spr_df.loc[spr_df[spr_df[kd_lookup] > lt].index, kd_lookup] = np.nan
    spr_df = spr_df[~spr_df[kd_lookup].isna()]
    spr_df = spr_df.sort_values("Chi2").groupby(["Ligand", "Analyte"]).head(1)
    kd_lookup = "Koff"

    fig, axes = plt.subplots(
        4,
        1,
        figsize=(9.5, 10),
        gridspec_kw={
            # "width_ratios": [0.5],
            "height_ratios": [0.1, 1, 0.1, 1],
            "hspace": 0.9,
            "wspace": 0.1,
            # "bottom": 0.,
            "top": 0.85,
            "left": 0.18,
            "right": 0.9,
        },
    )

    axes[0].letter = "A"
    axes[1].letter = "A"
    axes[2].letter = "B"
    axes[3].letter = "B"

    ax = axes[0]
    _df = spr_df.copy(deep=True)
    # breakpoint()
    table_ax, table_df = plot.plot_dyn_table(
        spr_df=_df,
        ax=ax,
        groupby=["pseudogroup", "is_cp"],
        lim=lt,
        KD_fix=kd_lookup,
        add_kd_lt=False,
        median_name="Median $\mathregular{K_{off}}$ (s$^{-1}$)",
        sci_notation=True,
        median_scale=1,
    )
    table_df = table_df.T.reset_index()
    table_df = data.populate_psname_spr(table_df)
    dfs.append(
        annotate_df(
            table_df.reset_index(), ax=ax, key=list(table_df.columns), label=f"header-Koff-{analyte_short_name}"
        )
    )
    print(kd_lookup)
    ax = axes[1]
    _df = spr_df.copy(deep=True)
    sns.stripplot(
        data=_df,
        x="pseudogroup",
        y=kd_lookup,
        dodge=True,
        hue="is_cp",
        edgecolor="black",
        linewidth=1,
        palette=pal,
        ax=ax,
    )
    sns.boxplot(
        data=_df,
        x="pseudogroup",
        y=kd_lookup,
        dodge=True,
        hue="is_cp",
        fliersize=0,
        palette=pal,
        whis=[10, 90],
        ax=ax,
    )
    dfs.append(annotate_df(_df, ax=ax, key=["is_cp", kd_lookup], label=f"Koff-{analyte_short_name}"))
    ax.set_yscale("log")
    ax.legend().set_visible(False)
    ax.set_ylabel(
        "$\mathregular{K_{off}}$ for " + analyte_short_name + " (s$^{-1}$)",
        size=13,
    )
    ax.set_xlabel("")
    adjust_boxplot(ax)
    ax.set_xticklabels([])
    ax.set_xticklabels(
        [
            r"eOD" + "\nwk8",
            r"eOD$\rightarrow$eOD" + "\nwk16",
            r"eOD$\rightarrow$core" + "\nwk 16",
            r"eOD$\rightarrow$eOD" + "\nwk 24",
            r"eOD$\rightarrow$core" + "\nwk 24",
            r"eOD$\rightarrow$eOD$\rightarrow$core" + "\nwk 24",
        ],
        rotation=0,
        verticalalignment="top",
        horizontalalignment="center",
    )
    # ax.axhline(y=5e-5, linestyle="--", linewidth=1, color="black", zorder=-10)
    custom_lines = []
    for x in {
        "Random": "#91FCC0",
        "Selected": "#6C65FF",
    }:
        custom_lines.append(
            Patch(
                facecolor={"Random": "#91FCC0", "Selected": "#6C65FF"}[x],
                edgecolor="black",
                linewidth=1,
                label=x,
            )
        )
    ax.legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=3,
        bbox_to_anchor=(0.5, -0.2),
        labelspacing=0.1,
    )
    # ax.set_ylim(1e-11, 0.00015)

    spr_df = data.get_g002_spr_df_boost()

    # Filters
    spr_df = spr_df.query(f"is_vrc01_class")
    spr_df = spr_df.query("top_c_call=='IGHG'")
    spr_df = spr_df.query(f"Analyte=='{analyte_long_name}'")

    # Transformations
    kd_lookup = "KD_fix"
    spr_df.loc[spr_df[spr_df["estimated"] == True].index, kd_lookup] = np.nan
    spr_df.loc[spr_df[spr_df[kd_lookup] > lt].index, kd_lookup] = np.nan
    spr_df = spr_df[~spr_df[kd_lookup].isna()]
    spr_df = spr_df.sort_values("Chi2").groupby(["Ligand", "Analyte"]).head(1)
    kd_lookup = "Kon"

    ax = axes[2]
    _df = spr_df.copy(deep=True)
    table_ax, table_df = plot.plot_dyn_table(
        spr_df=_df,
        ax=ax,
        lim=lt,
        groupby=["pseudogroup", "is_cp"],
        KD_fix=kd_lookup,
        add_kd_lt=False,
        median_name="Median $\mathregular{K_{on}}$ (M$^{-1}$s$^{-1}$)",
        median_scale=1,
        sci_notation=True,
    )
    table_df = table_df.T.reset_index()
    table_df = data.populate_psname_spr(table_df)
    dfs.append(
        annotate_df(table_df.reset_index(), ax=ax, key=list(table_df.columns), label=f"header-Kon-{analyte_short_name}")
    )

    ax = axes[3]
    _df = spr_df.copy(deep=True)
    sns.stripplot(
        data=_df,
        x="pseudogroup",
        y=kd_lookup,
        dodge=True,
        hue="is_cp",
        edgecolor="black",
        linewidth=1,
        palette=pal,
        ax=ax,
    )
    ax = sns.boxplot(
        data=_df,
        x="pseudogroup",
        y=kd_lookup,
        dodge=True,
        hue="is_cp",
        fliersize=0,
        palette=pal,
        whis=[10, 90],
        ax=ax,
    )
    dfs.append(annotate_df(_df, ax=ax, key=["is_cp", kd_lookup], label=f"Kon-{analyte_short_name}"))

    ax.set_yscale("log")
    # ax.set_ylim(1e-11, 0.00015)
    ax.legend().set_visible(False)
    # ax.yaxis.set_major_locator(mtick.LogLocator(base=10, numticks=9))
    # custom locator, shouldnt be linear anywhere
    # ax.yaxis.set_minor_locator(MinorSymLogLocator(linthresh=1))
    # ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: fixed_locator_lookup[x]))
    ax.set_ylabel(
        "$\mathregular{K_{on}}$ for " + analyte_short_name + " (M$^{-1}$s$^{-1}$)",
        size=13,
    )
    ax.set_xlabel("")
    adjust_boxplot(ax)
    ax.set_xticklabels(
        [
            r"eOD" + "\nwk8",
            r"eOD$\rightarrow$eOD" + "\nwk16",
            r"eOD$\rightarrow$core" + "\nwk 16",
            r"eOD$\rightarrow$eOD" + "\nwk 24",
            r"eOD$\rightarrow$core" + "\nwk 24",
            r"eOD$\rightarrow$eOD$\rightarrow$core" + "\nwk 24",
        ],
        rotation=0,
        verticalalignment="top",
        horizontalalignment="center",
    )

    custom_lines = []
    for x in {
        "Random": "#91FCC0",
        "Selected": "#6C65FF",
    }:
        custom_lines.append(
            Patch(
                facecolor={"Random": "#91FCC0", "Selected": "#6C65FF"}[x],
                edgecolor="black",
                linewidth=1,
                label=x,
            )
        )
    ax.legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=3,
        bbox_to_anchor=(0.5, -0.2),
        labelspacing=0.1,
    )

    sns.despine()
    plt.tight_layout()

    for i in [0, 2]:
        ax = axes[i]
        ax.text(-0.2, 1.1, ax.letter, transform=ax.transAxes, size=14, weight="bold")

    return fig, dfs


def plot_shm_spr(data: Data) -> tuple[pw.Bricks, list[pd.DataFrame]]:
    """Plot supplementary figure 32 A-D in a single 2x2 figure.

    Args:
        data: Data object containing required data
        outpath_base: Base output path for figures
        outpath_suffix: Output file suffix (default: 'png')
        dpi: DPI for output figures (default: 300)
    """

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

    dfs = []
    plot = Plot()
    g002_spr_df = data.get_g002_spr_df_boost()
    palette = {False: "#91FCC0", True: "#6C65FF"}

    # Common filtering
    g002_spr_df = g002_spr_df.query("top_c_call=='IGHG'").query("is_vrc01_class==True")
    g002_spr_df = g002_spr_df.query("Analyte=='core-Hx_r4.0D_TH6_g28v2_pHLsecAvi'")

    psuedogroups = [
        r"eOD" + "\nwk8",
        r"eOD$\rightarrow$eOD" + "\nwk16",
        r"eOD$\rightarrow$core" + "\nwk 16",
        r"eOD$\rightarrow$eOD" + "\nwk 24",
        r"eOD$\rightarrow$core" + "\nwk 24",
        r"eOD$\rightarrow$eOD$\rightarrow$core" + "\nwk 24",
    ]

    # Create 2x2 figure
    # fig, axes = plt.subplots(2, 2, figsize=(15, 15))
    # axes = axes.flatten()
    figsize = (3, 3)
    axes = [
        pw.Brick(figsize=figsize),
        pw.Brick(figsize=figsize),
        pw.Brick(figsize=figsize),
        pw.Brick(figsize=figsize),
    ]
    axes[0].letter = "A"
    axes[1].letter = "B"
    axes[2].letter = "C"
    axes[3].letter = "D"

    # Plot A - VH mutation
    ax = axes[0]
    ax, df = plot.plot_stripbox(
        df=g002_spr_df.copy(deep=True),
        x="pseudogroup",
        y="v_mutation_aa_heavy",
        hue="is_cp",
        xlabel="",
        ylabel=r"$\mathregular{V_H}$" + f" gene\n% mutation (aa)",
        remove_legend=True,
        remove_xticklabels=True,
        yaxis_mtick_major=0.02,
        yaxis_mtick_minor=0.01,
        palette=palette,
        tilt=True,
        ax=ax,
    )
    ax.set_ylim(0, 0.24)
    df = annotate_df(df, ax, ["v_mutation_aa_heavy", "is_cp"], "VH_mutation")
    dfs.append(df)
    adjust_boxplot(ax)

    # Plot B - VK/L mutation
    ax = axes[1]
    ax, df = plot.plot_stripbox(
        df=g002_spr_df.copy(deep=True),
        x="pseudogroup",
        y="v_mutation_aa_light",
        hue="is_cp",
        xlabel="",
        ylabel=r"$\mathregular{V_{K/L}}$" + f" gene\n%mutation (aa)",
        remove_legend=True,
        remove_xticklabels=True,
        yaxis_mtick_major=0.02,
        yaxis_mtick_minor=0.01,
        palette=palette,
        tilt=True,
        ax=ax,
    )
    ax.set_ylim(0, 0.24)
    df = annotate_df(df, ax, ["v_mutation_aa_light", "is_cp"], "VKL_mutation")
    dfs.append(df)
    adjust_boxplot(ax)

    # Plot C - Key VRC01-class HC residues
    ax = axes[2]
    ax, df = plot.plot_stripbox(
        df=g002_spr_df.copy(deep=True),
        x="pseudogroup",
        y="cottrell_focused_v_common_score",
        hue="is_cp",
        xlabel="",
        ylabel="Number of key VRC01-class\nHC residues",
        xticklabels=psuedogroups,
        yaxis_mtick_major=2,
        yaxis_mtick_minor=1,
        palette=palette,
        yaxis_mtick_normalize=True,
        ax=ax,
    )
    ax.set_ylim(0, 10)
    df = annotate_df(df, ax, ["cottrell_focused_v_common_score", "is_cp"], "VRC01_class_HC_residues")
    dfs.append(df)
    adjust_boxplot(ax)

    # Plot D - Key VRC01-class HCDR2 residues
    ax = axes[3]
    ax, df = plot.plot_stripbox(
        df=g002_spr_df.copy(deep=True),
        x="pseudogroup",
        y="num_hcdr2_mutations",
        hue="is_cp",
        xlabel="",
        ylabel="Number of key VRC01-class\nHCDR2 residues",
        xticklabels=psuedogroups,
        yaxis_mtick_major=2,
        yaxis_mtick_minor=1,
        palette=palette,
        yaxis_mtick_normalize=True,
        ax=ax,
    )
    ax.set_ylim(0, 10)
    df = annotate_df(df, ax, ["num_hcdr2_mutations", "is_cp"], "VRC01_class_HCDR2_residues")
    dfs.append(df)
    adjust_boxplot(ax)

    for ax in axes:
        ax.text(-0.2, 1.1, ax.letter, transform=ax.transAxes, size=14, weight="bold")

    palette = {"Random": "#91FCC0", "Selected": "#6C65FF"}
    plot.bottom_legend(ax=axes[2], palette=palette)
    plot.bottom_legend(ax=axes[3], palette=palette)

    g = pw.stack(
        [
            pw.stack([axes[0], axes[2]], margin=0.1, operator="/"),
            pw.stack([axes[1], axes[3]], margin=0.1, operator="/"),
        ],
        margin=0.1,
        operator="|",
    )

    return g, dfs
