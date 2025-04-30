from pathlib import Path

import pandas as pd
import patchworklib as pw
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import ticker as mtick
from matplotlib.patches import Patch

from g00x_figures.box_and_scatter.flow_frequencies import (
    MinorSymLogLocator,
    adjust_boxplot,
)
from g00x_figures.data import Data
from g00x_figures.plot_helpers.font import apply_global_font_settings


def plot_polyclonality(data: Data, outpath: str, metric_outdir: Path, use_last_week: bool = False) -> None:
    apply_global_font_settings(fontsize=18)
    # figure, axes = plt.subplots(3, 3, figsize=(8.5, 10))
    axes = [
        [
            pw.Brick(figsize=(4, 4)),
            pw.Brick(figsize=(4, 4)),
            pw.Brick(figsize=(4, 4)),
        ],
        [
            pw.Brick(figsize=(4, 4)),
            pw.Brick(figsize=(4, 4)),
            pw.Brick(figsize=(4, 4)),
        ],
        [
            pw.Brick(figsize=(4, 1)),
            # pw.Brick(figsize=(4, 4)),
            # pw.Brick(figsize=(4, 4)),
        ],
    ]

    trial_palette = {"G001": "#6C65FF", "G002": "#91FCC0", "G003": "#E377C2"}
    if use_last_week:
        weeks = [4, 8, 10, 16, 21, 24]
    else:
        weeks = [4, 8, 16]
    # histogram
    clustered_g001 = (
        data.get_clustered_g001_seqs().query(f"is_vrc01_class==True")
        # data.get_g001_cluster_ses().query(f"is_vrc01_class==True")
        # .query("weeks_post!='V02'")
        .query(f"weeks_post in {weeks}")
    )
    # breakpoint()
    clustered_g002 = (
        data.get_g002_sequences_prime(use_cluster_file=True, use_cluster_with_wk24=use_last_week)
        .query(f"is_vrc01_class==True")
        .query(f"weeks in {weeks}")
        .query("top_c_call=='IGHG'")
    )
    clustered_g002["weeks"] = clustered_g002["weeks"].replace({24: 21})
    clustered_g002["ptid"] = clustered_g002["pubID"]

    clustered_g003 = (
        data.get_g003_sequences_prime(use_cluster_file=True, use_cluster_with_wk21=use_last_week)
        .query(f"is_vrc01_class==True")
        .query(f"weeks in {weeks}")
        .query("top_c_call=='IGHG'")
    )
    clustered_g003["ptid"] = clustered_g003["pubID"]

    # plottable
    plottable_df_g001 = clustered_g001.groupby(["cluster"], as_index=False).apply(
        lambda x: pd.Series({"num_clusters": len(x)})
    )
    plottable_df_g002 = clustered_g002.groupby(["cluster"], as_index=False).apply(
        lambda x: pd.Series({"num_clusters": len(x)})
    )
    plottable_df_g003 = clustered_g003.groupby(["cluster"], as_index=False).apply(
        lambda x: pd.Series({"num_clusters": len(x)})
    )
    # first axis
    ax = axes[0][0]
    sns.histplot(
        data=plottable_df_g001,
        x="num_clusters",
        edgecolor="black",
        linewidth=1,
        binwidth=1,
        ax=ax,
        color=trial_palette["G001"],
    )
    ax.set_xlim(1, 30)
    ax.set_yscale("log")
    ax.set_xlabel("Cluster Size")
    # second axis
    ax = axes[0][1]
    sns.histplot(
        data=plottable_df_g002,
        x="num_clusters",
        edgecolor="black",
        linewidth=1,
        binwidth=1,
        ax=ax,
        color=trial_palette["G002"],
    )
    ax.set_xlim(1, 30)
    ax.set_yscale("log")
    ax.set_xlabel("Cluster Size")
    ax.set_ylabel("")

    # G003 Count
    ax = axes[0][2]
    sns.histplot(
        data=plottable_df_g003,
        x="num_clusters",
        edgecolor="black",
        linewidth=1,
        binwidth=1,
        ax=ax,
        color=trial_palette["G003"],
    )
    ax.set_xlim(1, 30)
    ax.set_yscale("log")
    ax.set_xlabel("Cluster Size")
    ax.set_ylabel("")

    # clones vs seqs
    plottable_df_g001 = clustered_g001.groupby(["pubid", "Dose_Group", "vaccine_group"], as_index=False).apply(
        lambda x: pd.Series({"num_clusters": len(x["cluster"].unique()), "num_seqs": len(x)})
    )
    plottable_df_g002 = clustered_g002.groupby(["ptid", "group"], as_index=False).apply(
        lambda x: pd.Series({"num_clusters": len(x["cluster"].unique()), "num_seqs": len(x)})
    )
    plottable_df_g003 = clustered_g003.groupby(["ptid", "group"], as_index=False).apply(
        lambda x: pd.Series({"num_clusters": len(x["cluster"].unique()), "num_seqs": len(x)})
    )
    # plottable_df_g001
    # from IPython import embed

    # embed()

    # third axis
    ax = axes[1][0]
    sns.scatterplot(
        data=plottable_df_g001,
        x="num_seqs",
        y="num_clusters",
        edgecolor="black",
        s=40,
        color=trial_palette["G001"],
        linewidth=1,
        ax=ax,
    )
    X = range(-1, 101)
    Y = range(-1, 101)
    ax.plot(
        X,
        Y,
        linewidth=2,
        linestyle="--",
        zorder=0,
        markeredgecolor="black",
        color=trial_palette["G001"],
    )
    ax.set_xlabel("Number of Sequences")
    ax.set_ylabel("Number of Clusters")
    # print(plottable_df_g002.groupby("num_clusters")["num_seqs"].median())

    # fourth axis
    ax = axes[1][1]
    sns.scatterplot(
        data=plottable_df_g002,
        x="num_seqs",
        y="num_clusters",
        edgecolor="black",
        color=trial_palette["G002"],
        s=40,
        linewidth=1,
        ax=ax,
    )
    X = range(-10, 1500)
    Y = range(-10, 1500)
    ax.plot(
        X,
        Y,
        linewidth=2,
        linestyle="--",
        zorder=0,
        markeredgecolor="black",
        color=trial_palette["G002"],
    )
    ax.set_xlabel("Number of Sequences")
    ax.set_ylabel("")
    if use_last_week:
        ax.set_xlim(xmax=2500)
        ax.set_ylim(ymax=2500)

    # G003 num of clusters
    ax = axes[1][2]
    sns.scatterplot(
        data=plottable_df_g003,
        x="num_seqs",
        y="num_clusters",
        edgecolor="black",
        color=trial_palette["G003"],
        s=40,
        linewidth=1,
        ax=ax,
    )
    X = range(-10, 2000)
    Y = range(-10, 2000)
    ax.plot(
        X,
        Y,
        linewidth=2,
        linestyle="--",
        zorder=0,
        markeredgecolor="black",
        color=trial_palette["G003"],
    )
    ax.set_xlabel("Number of Sequences")
    ax.set_ylabel("")
    ax.set_xlim(xmax=2000)
    ax.set_ylim(ymax=2000)

    # clonality
    plottable_df_g001 = (
        clustered_g001.groupby(["pubid", "weeks_post", "timepoint"])
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
    plottable_df_g001["ptid"] = plottable_df_g001["pubid"]
    plottable_df_g001["weeks"] = plottable_df_g001["weeks_post"].astype(int)

    clustered_g002["timepoint"] = clustered_g002["visit_id"]
    plottable_df_g002 = (
        clustered_g002.groupby(["ptid", "weeks", "timepoint"])
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
    plottable_df_g003 = (
        clustered_g003.groupby(["ptid", "weeks", "timepoint"])
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
    plottable_df_g001["trial"] = "G001"
    plottable_df_g002["trial"] = "G002"
    plottable_df_g003["trial"] = "G003"
    # breakpoint()
    plottable = pd.concat([plottable_df_g001, plottable_df_g002, plottable_df_g003])
    plottable[["trial", "ptid", "weeks", "timepoint", "clonality"]].to_csv(
        metric_outdir / "figC_polyclonality.csv", index=False
    )
    unique_timepoints = plottable.groupby(["trial", "weeks"])["timepoint"].unique().reset_index()
    unique_timepoints.to_csv(metric_outdir / "unique_timepoints.csv", index=False)
    # from IPython import embed

    # embed()

    # fifth axis
    ax = axes[2][0]
    sns.stripplot(
        data=plottable,
        x="weeks",
        y="clonality",
        linewidth=1,
        edgecolor="black",
        hue="trial",
        dodge=True,
        # size=7,
        s=7,
        # hue_order=["G001", "G002", "G003"],
        palette=trial_palette,
        # color=trial_palette["G001"],
        jitter=0.15,
        ax=ax,
    )
    sns.boxplot(
        data=plottable,
        x="weeks",
        y="clonality",
        fliersize=0,
        linewidth=1,
        hue="trial",
        # dodge=False,
        # color=trial_palette["G001"],
        # hue_order=["G001", "G002", "G003"],
        palette=trial_palette,
        whis=[10, 90],
        ax=ax,
    )

    adjust_boxplot(ax)
    ax.set_ylabel("Fraction Unique Clones")
    ax.set_xlabel("Weeks Post-Vaccination")
    ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
    ax.set_ylim(0, 1.05)

    custom_lines = []
    # pal = {"G001": "#6C65FF", "G002": "#91FCC0"}
    for x in trial_palette:
        custom_lines.append(
            Patch(
                facecolor=trial_palette[x],
                edgecolor="black",
                linewidth=1,
                label=x,
            )
        )
    axes[2][0].legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=3,
        bbox_to_anchor=(0.5, -0.17),
        labelspacing=0.1,
    )

    for ax in [ax for row in axes for ax in row]:
        sns.despine(ax=ax)
    g1, g2, g3 = axes[0]
    g4, g5, g6 = axes[1]
    g7 = axes[2][0]
    if use_last_week:
        g7.set_xticklabels(
            ["4", "8", "10", "16", "24/21"],
        )
    else:
        g7.set_xticklabels(
            ["4", "8", "16"],
        )

    g1.set_title("A", fontsize=20, fontweight="bold", x=-0.2, y=1)
    g4.set_title("B", fontsize=20, fontweight="bold", x=-0.2, y=1)
    g7.set_title("C", fontsize=20, fontweight="bold", x=-0.06, y=1)
    # g7 = pw.spacer(g7, 0.25) | g7 | pw.spacer(g7, 0.25)
    g = (g1 | g2 | g3) / (g4 | g5 | g6) / g7
    # plt.tight_layout()
    g.savefig(outpath, dpi=300)


def plot_multi_find_clonality(
    data: Data,
    outpath: str,
    use_last_week: bool = False,
) -> None:
    # g001_seqs = data.get_g001_cluster_ses().query(f"is_vrc01_class")
    # .query("weeks_post != 'V02'")
    # g002_seqs =  pd.read_feather(data.paths.g002_sequences).query("is_vrc01_class").query("weeks != -5")
    if use_last_week:
        weeks = [4, 8, 10, 16, 21, 24]
    else:
        weeks = [4, 8, 16]
    g001_seqs = data.get_clustered_g001_seqs().query(f"is_vrc01_class==True").query(f"weeks_post in {weeks}")
    g002_seqs = (
        data.get_g002_sequences_prime(use_cluster_file=True, use_cluster_with_wk24=use_last_week)
        .query(f"is_vrc01_class==True")
        .query(f"weeks in {weeks}")
        .query("top_c_call=='IGHG'")
    )
    g003_seqs = (
        data.get_g003_sequences_prime(use_cluster_file=True, use_cluster_with_wk21=use_last_week)
        .query(f"is_vrc01_class==True")
        .query(f"weeks in {weeks}")
        .query("top_c_call=='IGHG'")
    )

    cluster_summary_g003 = (
        g003_seqs.groupby("cluster")
        .apply(
            lambda x: pd.Series(
                {
                    "size": len(x),
                    "subjects": list(x["ptid"].unique()),
                    "num_subjects": len(list(x["ptid"].unique())),
                    "timepoints": list(x["weeks"].unique()),
                    "num_timeponts": len(list(x["weeks"].unique())),
                    "cdr3_aa_heavy": x["cdr3_aa_heavy"].to_list(),
                }
            )
        )
        .sort_values("num_subjects")[::-1]
    )

    cluster_summary_g002 = (
        g002_seqs.groupby("cluster")
        .apply(
            lambda x: pd.Series(
                {
                    "size": len(x),
                    "subjects": list(x["ptid"].unique()),
                    "num_subjects": len(list(x["ptid"].unique())),
                    "timepoints": list(x["weeks"].unique()),
                    "num_timeponts": len(list(x["weeks"].unique())),
                    "cdr3_aa_heavy": x["cdr3_aa_heavy"].to_list(),
                }
            )
        )
        .sort_values("num_subjects")[::-1]
    )

    cluster_summary_g001 = (
        g001_seqs.groupby("cluster")
        .apply(
            lambda x: pd.Series(
                {
                    "size": len(x),
                    "subjects": list(x["pubid"].unique()),
                    "num_subjects": len(list(x["pubid"].unique())),
                    "timepoints": list(x["weeks_post"].unique()),
                    "num_timeponts": len(list(x["weeks_post"].unique())),
                    "cdr3_aa_heavy": x["cdr3_aa_heavy"].to_list(),
                }
            )
        )
        .sort_values("num_subjects")[::-1]
    )
    plottable_df = []

    # G003
    number_singletons = len(cluster_summary_g003.query("num_timeponts==1 and num_subjects==1").index.unique())
    number_timepoints = len(cluster_summary_g003.query("num_timeponts > 1").index.unique())
    number_donors = len(cluster_summary_g003.query("num_subjects > 1").index.unique())
    plottable_df += [
        {
            "count": number_singletons,
            "label": "cluster_single_timepoint",
            "trial": "G003",
        },
        {
            "count": number_timepoints,
            "label": "cluster_more_one_timepoint",
            "trial": "G003",
        },
        {
            "count": number_donors,
            "label": "cluster_multiple_subjects",
            "trial": "G003",
        },
    ]

    # G002
    number_singletons = len(cluster_summary_g002.query("num_timeponts==1 and num_subjects==1").index.unique())
    number_timepoints = len(cluster_summary_g002.query("num_timeponts > 1").index.unique())
    number_donors = len(cluster_summary_g002.query("num_subjects > 1").index.unique())
    plottable_df += [
        {
            "count": number_singletons,
            "label": "cluster_single_timepoint",
            "trial": "G002",
        },
        {
            "count": number_timepoints,
            "label": "cluster_more_one_timepoint",
            "trial": "G002",
        },
        {
            "count": number_donors,
            "label": "cluster_multiple_subjects",
            "trial": "G002",
        },
    ]
    # from IPython import embed

    # embed()
    # G001
    number_singletons = len(cluster_summary_g001.query("num_timeponts==1 and num_subjects==1").index.unique())
    number_timepoints = len(cluster_summary_g001.query("num_timeponts > 1").index.unique())
    number_donors = len(cluster_summary_g001.query("num_subjects > 1").index.unique())
    plottable_df += [
        {
            "count": number_singletons,
            "label": "cluster_single_timepoint",
            "trial": "G001",
        },
        {
            "count": number_timepoints,
            "label": "cluster_more_one_timepoint",
            "trial": "G001",
        },
        {
            "count": number_donors,
            "label": "cluster_multiple_subjects",
            "trial": "G001",
        },
    ]
    plottable_df = pd.DataFrame(plottable_df)
    print(plottable_df)

    fig, ax = plt.subplots(1, 1, figsize=(4, 5))
    ax = sns.barplot(
        data=plottable_df,
        y="count",
        x="label",
        hue="trial",
        edgecolor="black",
        hue_order=["G001", "G002", "G003"],
        palette={"G001": "#6C65FF", "G002": "#91FCC0", "G003": "#E377C2"},
    )
    ax.set_ylim(1, 10e4)
    ax.set_yscale("symlog", linthresh=10)
    ax.set_xticklabels(
        [
            "Single Donor\n and Timepoint",
            "Multiple\nTimepoints",
            "Multiple\nDonors",
        ],
        rotation=90,
        ha="right",
        va="center",
        rotation_mode="anchor",
    )
    ax.set_ylabel("Number of Clusters", labelpad=0)
    ax.set_xlabel("")
    ax.yaxis.set_minor_locator(MinorSymLogLocator(10))
    custom_lines = []
    pal = {"G001": "#6C65FF", "G002": "#91FCC0", "G003": "#E377C2"}
    for x in pal:
        custom_lines.append(Patch(facecolor=pal[x], edgecolor="black", linewidth=1, label=x))
    ax.legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=3,
        bbox_to_anchor=(0.5, -0.5),
        labelspacing=0.1,
    )
    sns.despine()
    plt.tight_layout()
    fig.savefig(outpath, dpi=300)
    return plottable_df


def plot_boost_clonality(
    data: Data,
    img_outdir: Path,
    metric_outdir: Path,
    fig: str,
    y_axis: str = "total_unique",
) -> None:
    apply_global_font_settings()
    clustered_g002 = (
        data.get_g002_sequences_boost(use_cluster_file=True)
        .query(f"is_vrc01_class==True")
        .query("weeks != -5")
        .query("top_c_call=='IGHG'")
    )

    # First calculate median cluster size per participant
    cluster_sizes = (
        clustered_g002.groupby(["pseudogroup", "pubID", "weeks", "cluster"]).size().reset_index(name="cluster_size")
    )
    cluster_sizes.to_csv(metric_outdir / f"{fig}_cluster_report.csv", index=False)
    median_cluster_sizes = (
        cluster_sizes.groupby(["pseudogroup", "pubID", "weeks"])["cluster_size"]
        .median()
        .reset_index(name="median_cluster_size")
    )

    plottable_df_g002 = (
        clustered_g002.groupby(["pseudogroup", "pubID", "weeks"])
        .apply(
            lambda x: pd.Series(
                {
                    "clonality": len(x["cluster"].unique()) / len(x),
                    "total_unique": len(x["cluster"].unique()),
                    "total_non_unique": len(x) - len(x["cluster"].unique()),
                    "num_seqs": len(x),
                    "num_clusters": len(x["cluster"].unique()),
                }
            )
        )
        .reset_index()
    )

    # Merge median cluster size into plottable_df_g002
    plottable_df_g002 = plottable_df_g002.merge(median_cluster_sizes, on=["pseudogroup", "pubID", "weeks"], how="left")
    # from IPython import embed; embed()
    ax = pw.Brick(figsize=(6, 4))
    plottable_df_g002 = data.populate_psname(plottable_df_g002)
    plottable_df_g002[["psname", "pseudogroup", "pubID", "weeks", y_axis]].sort_values(["pseudogroup"]).to_csv(
        metric_outdir / f"{fig}_boost_clonality_{y_axis}.csv", index=False
    )
    median = plottable_df_g002.groupby(["psname", "pseudogroup"])[y_axis].median().round(3)
    median = median.reset_index()
    median.sort_values(["pseudogroup"]).to_csv(metric_outdir / f"{fig}_medians_{y_axis}.csv", index=False)

    sns.boxplot(
        data=plottable_df_g002,
        x="pseudogroup",
        y=y_axis,
        dodge=False,
        fliersize=0,
        linewidth=1,
        palette=data.get_week_palette(),
        hue="weeks",
        whis=[10, 90],
        ax=ax,
    )
    sns.stripplot(
        data=plottable_df_g002,
        x="pseudogroup",
        y=y_axis,
        linewidth=1,
        dodge=False,
        hue="weeks",
        edgecolor="black",
        palette=data.get_week_palette(),
        size=7,
        jitter=0.1,
        ax=ax,
    )
    adjust_boxplot(ax)
    ax.get_legend().remove()
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
        ha="right",
        rotation_mode="anchor",
    )

    ax.set_xlabel("")
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
    ax.legend(
        custom_lines,
        [f"wk {i._label}" for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=1,
        fontsize=18,
        bbox_to_anchor=(1, -0.15),
        labelspacing=0.1,
        # title="Weeks",
        # title_fontsize=18,
    )
    if y_axis == "clonality":
        ax.set_ylabel("Fraction Unique VRC01-class Clones", labelpad=10)
    elif y_axis == "total_unique":
        ax.set_ylabel("Total Unique VRC01-class Clones")
        ax.set_yscale("symlog", linthresh=10)
        ax.yaxis.set_minor_locator(MinorSymLogLocator(10))
        # ax.minorticks_on()
        ax.set_ylim(0, 1000)
    elif y_axis == "total_non_unique":
        ax.set_yscale("symlog", linthresh=10)
        ax.set_ylabel("Total Non-Unique VRC01-class Clones")
        ax.yaxis.set_minor_locator(MinorSymLogLocator(10))
        ax.set_ylim(0, 1000)
    elif y_axis == "num_seqs":
        ax.set_yscale("symlog", linthresh=10)
        ax.set_ylabel("Total VRC01-class Sequences")
        ax.yaxis.set_minor_locator(MinorSymLogLocator(10))
        ax.set_ylim(0, 1000)

    sns.despine(ax=ax)
    return ax
