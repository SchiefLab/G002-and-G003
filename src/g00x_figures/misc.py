# from pathlib import Path

# import pandas as pd
# import seaborn as sns
# from Levenshtein import distance
# from matplotlib import pyplot as plt
# from matplotlib import ticker as mtick
# from matplotlib.patches import Patch

# from g00x_figures.box_and_scatter.flow_frequencies import adjust_boxplot
# from g00x_figures.data import Data


# def plot_methodology_comparision(data: Data, output: str):
#     df = data.get_g001_methodology()
#     ax = sns.stripplot(
#         data=df,
#         x="method",
#         y="Percent of IgG+ B cells detected as VRC01-class (missing seq to 0)",
#         edgecolor="black",
#         order=["Sanger", "10X"],
#         linewidth=1,
#     )
#     ax = sns.boxplot(
#         data=df,
#         x="method",
#         y="Percent of IgG+ B cells detected as VRC01-class (missing seq to 0)",
#         fliersize=0,
#         order=["Sanger", "10X"],
#     )
#     plt.savefig(output + ".png", dpi=300)


from pathlib import Path
from typing import List

import pandas as pd
import patchworklib as pw
import seaborn as sns
from Levenshtein import distance
from matplotlib import pyplot as plt
from matplotlib import ticker as mtick
from matplotlib.patches import Patch

from g00x_figures.box_and_scatter.flow_frequencies import adjust_boxplot
from g00x_figures.data import Data
from g00x_figures.flow_frequencies import generic_plot_box_and_whisker
from g00x_figures.plot_helpers.font import apply_global_font_settings
from g00x_figures.plot_helpers.legend import plot_legend

apply_global_font_settings()
METHODOLOGY_pubIDs = [
    "PubID_023",
    "PubID_005",
    "PubID_051",
    "PubID_187",
    "PubID_153",
    "PubID_079",
    "PubID_047",
    "PubID_154",
    "PubID_110",
    "PubID_046",
]


def get_g003_methodology(data: Data):
    # TODO: add path to data
    g003_method_flow = data.get_g003_methology_flow_and_seq()
    g003_method_flow["trial"] = "G003"

    # g003_method_flow = g003_method_flow.query("run_purpose == 'G001Sort'")

    column_g003_2_g001_mapping = {
        "Percent antigen-specific among IgG^{+}": "percentage_of_antigen++",
        # "Percent epitope-specific (KO^-Ag^{++}) among IgD^-": "percentage_of_epitope_specific",
        "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}": "percentage_of_epitope_specific",
        "pubID": "PTID",
        "Percent of VRC01-class sequences among IgG": "Percent of IgG+ B cells detected as VRC01-class (missing seq to 0)",
        # "v_mutation_aa_heavy": "v_mutation_aa_heavy",
    }
    for column_g003, column_g001 in column_g003_2_g001_mapping.items():
        g003_method_flow[column_g001] = g003_method_flow[column_g003]

    g003_method_flow["method"] = g003_method_flow["run_purpose"].apply(lambda x: "KWTRP" if x == "G001Sort" else "NAC")
    return g003_method_flow


def get_g003_seq_methods(data: Data):
    # TODO: add path to data
    g003_seq_methods = data.get_g003_methodology_seq()

    # g003_seq_methods = g003_seq_methods[g003_seq_methods.run_purpose == "G001Sort"]
    # g003_seq_methods = g003_seq_methods[
    #     g003_seq_methods["PTID"].isin(g003_seq_methods.PTID.unique())
    # ]
    g003_seq_methods = g003_seq_methods.query("timepoint=='V08'")
    # g003_seq_methods["method"] = "10X_G003"
    g003_seq_methods["method"] = g003_seq_methods["run_purpose"].apply(lambda x: "KWTRP" if x == "G001Sort" else "NAC")
    g003_seq_methods.PTID.unique()
    return g003_seq_methods


def get_g00x_flow_and_seq(data: Data):
    # successful sequences when using 10x G002 methodology
    g002_METHODOLOGY_pubIDs = [
        "PubID_023",
        "PubID_005",
        "PubID_051",
        "PubID_187",
        "PubID_153",
        "PubID_079",
        "PubID_047",
        "PubID_154",
        "PubID_110",
        "PubID_046",
    ]
    # successful sequences when using 10x G003 methodology
    g003_METHODOLOGY_pubIDs = [
        # "PubID_079",  # also did not make it
        "PubID_187",
        "PubID_153",
        "PubID_046",
    ]
    mg001_flow_and_seq = data.get_g001_flow_and_seq_prime(allow_week_10=True)
    mg001_flow_and_seq = mg001_flow_and_seq.query("pubID in @g002_METHODOLOGY_pubIDs")

    mg001_flow_and_seq = mg001_flow_and_seq.query("weeks == 10")
    mg001_flow_and_seq["method"] = "G001"
    mg001_flow_and_seq["trial"] = "G001"

    mg002_flow_and_seq = data.get_g002_methodology_flow_and_seq()
    mg002_flow_and_seq["method"] = "G002"
    mg002_flow_and_seq["trial"] = "G002"

    gb_columns: list[str] = [
        "pubID",
        "Visit",
        "Probeset",
        "SampleType",
    ]
    isotype = "IGHG"
    column = "top_c_call"
    df = data.get_g001_10x_sequences()
    df["is_vrc01_class"] = df["is_vrc01_class"].astype(bool)
    g_df = df.query(f"{column}=='{isotype}'").groupby(gb_columns)
    break_down = (
        g_df["is_vrc01_class"]
        .value_counts(normalize=True)
        .reset_index(name="Percent of IGHG sequences that are VRC01-class")
        .query("is_vrc01_class")
    )

    _df = break_down.set_index("pubID")[["Percent of IGHG sequences that are VRC01-class"]]

    mg002_flow_and_seq = pd.concat([_df, mg002_flow_and_seq.set_index("pubID")], axis=1).reset_index()
    mg002_flow_and_seq["Percent of VRC01-class sequences among IgG"] = (
        mg002_flow_and_seq["percent_ep_among_igg"]
        * mg002_flow_and_seq["Percent of IGHG sequences that are VRC01-class"]
    )

    mg003_flow_and_seq = data.get_g003_methodology_flow_and_seq()
    mg003_flow_and_seq["method"] = "G003"
    mg003_flow_and_seq["trial"] = "G003"

    mg001_2_flow_and_seq = pd.concat(
        [
            mg001_flow_and_seq.query("pubID in @g002_METHODOLOGY_pubIDs"),
            mg002_flow_and_seq.query("pubID in @g002_METHODOLOGY_pubIDs"),
        ]
    ).reset_index(drop=True)
    mg001_2_3_flow_and_seq = pd.concat(
        [
            mg001_flow_and_seq,
            mg002_flow_and_seq,
            mg003_flow_and_seq,
            # mg001_flow_and_seq.query("pubID in @g003_METHODOLOGY_pubIDs"),
            # mg002_flow_and_seq.query("pubID in @g003_METHODOLOGY_pubIDs"),
            # mg003_flow_and_seq.query("pubID in @g003_METHODOLOGY_pubIDs"),
        ]
    ).reset_index(drop=True)
    mg001_2_3_flow_and_seq = mg001_2_3_flow_and_seq[
        mg001_2_3_flow_and_seq.pubID.isin(mg003_flow_and_seq.pubID.unique())
    ]

    mapping = {
        "Percent of IgG+ B cells that are GT8++ (without regard to KO binding status)": "Percent antigen-specific among IgG^{+}",
        "Percent of IgG+ B cells that are epitope-specific (KO-GT8++)": "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}",
        "Percent of GT8++IgG+ B cells that are KO-": "Percent IgG^{+}KO^- among Ag^{++}",
        "Number of epitope-specific (KO-GT8++) sequenced IgG BCRs that are VRC01-class": "Number of IGHG sequences that are VRC01-class",
        "Percent of IgG+ B cells detected as VRC01-class (missing seq to 0)": "Percent of VRC01-class sequences among IgG",
        "Response (missing seq to 0)": "Response x",
        "Percent of epitope-specific (KO-GT8++) sequenced IgG BCRs that are VRC01-class": "Percent of IGHG sequences that are VRC01-class",
        "Percent of GT8++ IgG+ B cells detected as VRC01-class (missing seq to 0)": "Percent VRC01-class among eOD-specific IgG+ memory BCR sequences",
    }
    for k, v in mapping.items():
        if k not in mg001_flow_and_seq.columns:
            mg001_2_3_flow_and_seq[v] = mg001_2_3_flow_and_seq[k]
        if k not in mg001_2_flow_and_seq.columns:
            mg001_2_flow_and_seq[v] = mg001_2_flow_and_seq[k]
        if v not in mg001_2_3_flow_and_seq.columns:
            mg001_2_3_flow_and_seq[v] = mg001_2_3_flow_and_seq[k]
        if v not in mg001_2_flow_and_seq.columns:
            mg001_2_flow_and_seq[v] = mg001_2_flow_and_seq[k]

    a = "Percent of IGHG sequences that are VRC01-class"
    b = "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}"
    c = "Percent antigen-specific among IgG^{+}"
    d = "Percent VRC01-class among eOD-specific IgG+ memory BCR sequences"

    mg001_2_3_flow_and_seq[d] = mg001_2_3_flow_and_seq[a] * mg001_2_3_flow_and_seq[b] / mg001_2_3_flow_and_seq[c]
    mg001_2_flow_and_seq[d] = mg001_2_flow_and_seq[a] * mg001_2_flow_and_seq[b] / mg001_2_flow_and_seq[c]

    return mg001_2_flow_and_seq, mg001_2_3_flow_and_seq


def plot_g003_methodology_comparision(data: Data, axes: List[plt.Axes], sup_metric_outdir: Path, img_outdir: Path):
    g003_method_flow = get_g003_methodology(data=data)
    method_df = data.get_g001_methodology()
    ptid_hcl = [
        # "PubID_079",
        "PubID_187",
        "PubID_153",
        "PubID_046",
    ]
    drop_ptids = set(method_df.PTID.unique()) - set(g003_method_flow.PTID.unique())
    method_df = pd.concat(
        [
            method_df[~method_df.PTID.isin(drop_ptids)],
            g003_method_flow[
                [
                    "run_purpose",
                    "PTID",
                    "percentage_of_antigen++",
                    "method",
                    "percentage_of_epitope_specific",
                    "Percent of IgG+ B cells detected as VRC01-class (missing seq to 0)",
                    "num_vrc01_class",
                    "num_not_vrc01_class",
                    "IgD-IgM-IgG+ B cells",
                    "IgM-/IgG+/KO-/Antigen++ B cells",
                    "Number of IGHG sequences that are VRC01-class",
                ]
            ],
        ],
        axis=0,
    )

    pal = {
        "Sanger": "#6C65FF",
        "10X": "#91FCC0",
        "KWTRP": "#E377C2",
    }
    df = method_df.copy(deep=True).query('method!="NAC"')
    order = ["Sanger", "10X", "KWTRP"]
    df["trial"] = df["method"].replace({"Sanger": "G001", "10X": "G002", "KWTRP": "G003"})

    metrics = [
        (
            "percentage_of_antigen++",
            "% GT8++\namong IgG+ B cells",
            "figE percent GT8++ among IgG+ B cells",
            df,
        ),
        (
            "percentage_of_epitope_specific",
            "% CD4bs-specific\namong IgG+ B cells",
            "figF percent CD4bs-specific among IgG+ B cells",
            df,
        ),
        (
            "Percent of IgG+ B cells detected as VRC01-class (missing seq to 0)",
            "Percent of IgG+ B cells\ndetected as VRC01-class",
            "figG Percent of IgG+ B cells detected as VRC01-class",
            df[df.PTID.isin(ptid_hcl)],
        ),
    ]

    for i, (metric, ylabel, yname, df) in enumerate(metrics):
        ax = axes[i]
        sns.stripplot(
            data=df,
            x="method",
            y=metric,
            edgecolor="black",
            order=order,
            linewidth=1,
            s=7,
            palette=pal,
            ax=ax,
        )
        sns.boxplot(
            data=df,
            x="method",
            y=metric,
            fliersize=0,
            order=order,
            palette=pal,
            ax=ax,
        )
        adjust_boxplot(ax)
        ax.set_xticklabels([])
        ax.set_xlabel("")
        ax.set_ylabel(ylabel)
        df[["PTID", "trial", "method", metric]].to_csv(
            sup_metric_outdir / f"{yname.replace(' ', '_')}.csv", index=False
        )

    g003_seq_methods = get_g003_seq_methods(data=data)
    ptids = ptid_hcl
    sequences_g001_10x = data.get_g001_10x_sequences()
    sequences_g001_10x["method"] = "10X"
    # sequence_g002_sanger = pd.read_feather(data.paths.g001_sequences)
    sequence_g002_sanger = data.get_g001_seqs()
    sequence_g002_sanger = sequence_g002_sanger[sequence_g002_sanger["PTID"].isin(ptids)]
    sequence_g002_sanger = sequence_g002_sanger.query("Timepoint=='V08'")
    sequence_g002_sanger["method"] = "Sanger"
    combined = pd.concat(
        [
            sequence_g002_sanger,
            sequences_g001_10x,
            g003_seq_methods,
        ]
    ).reset_index(drop=True)
    combined = combined.query("is_vrc01_class == True")
    combined = combined.query("top_c_call=='IGHG'")
    seq_df = combined.groupby(["PTID", "method"])["v_mutation_aa_heavy"].median().to_frame().reset_index()

    seq_df = seq_df[seq_df.PTID.isin(ptid_hcl)].query('method!="NAC"')

    ax = axes[3]
    sns.stripplot(
        data=seq_df,
        x="method",
        y="v_mutation_aa_heavy",
        edgecolor="black",
        # order=order,
        linewidth=1,
        s=7,
        palette=pal,
        ax=ax,
    )
    sns.boxplot(
        data=seq_df,
        x="method",
        y="v_mutation_aa_heavy",
        fliersize=0,
        # order=order,
        palette=pal,
        ax=ax,
    )
    adjust_boxplot(ax)
    ax.set_xticklabels(["G001 Method", "G002 Method", "G003 KWTRP"])
    ax.set_xlabel("")
    ylabel = r"$\mathregular{V_H}$" + f" gene\n%mutation (aa)"
    yname = "H V Heavy Gene Mutation AA".replace(" ", "_")
    ax.set_ylabel(ylabel)
    seq_df["trial"] = seq_df["method"].replace({"Sanger": "G001", "10X": "G002", "KWTRP": "G003"})
    seq_df[
        [
            "PTID",
            "trial",
            "method",
            "v_mutation_aa_heavy",
        ]
    ].to_csv(sup_metric_outdir / f"{yname}.csv", index=False)

    for ax, title in zip(axes, ["E", "F", "G", "H"]):
        ax.set_title(title, fontsize=16, fontweight="bold", x=-0.1, y=1.1)


def plot_methodology_comparision(data: Data, metric_outdir: Path, img_outdir: Path):
    apply_global_font_settings(14)

    df = data.get_g001_methodology()
    df["trial"] = df["method"].replace({"Sanger": "G001", "10X": "G002", "KWTRP": "G003"})

    pal = {"Sanger": "#6C65FF", "10X": "#91FCC0"}
    fig, axes = plt.subplots(4, 2, figsize=(8.5, 10))
    # change plot with patchworklib
    axes = [
        [pw.Brick(figsize=(4, 2)), pw.Brick(figsize=(4, 2))],
        [pw.Brick(figsize=(4, 2)), pw.Brick(figsize=(4, 2))],
        [pw.Brick(figsize=(4, 2)), pw.Brick(figsize=(4, 2))],
        [pw.Brick(figsize=(4, 2)), pw.Brick(figsize=(4, 2))],
        [pw.Brick(figsize=(4, 2)), pw.Brick(figsize=(4, 2))],
    ]
    faxes = [ax for row in axes for ax in row]

    metrics = [
        (
            "percentage_of_antigen++",
            "% GT8++\namong IgG+ B cells",
            "figA percent GT8++ among IgG+ B cells",
        ),
        (
            "percentage_of_epitope_specific",
            "% CD4bs-specific\namong IgG+ B cells",
            "figB percent CD4bs-specific among IgG+ B cells",
        ),
        (
            "percent_iggko_among_ag",
            "% of GT8++ IgG+ B cells that are KO-",
            "figC percent CD4bs-specific among IgG+ B cells",
        ),
        (
            "Percent of IgG+ B cells detected as VRC01-class (missing seq to 0)",
            "Percent of IgG+ B cells\ndetected as VRC01-class",
            "figD Percent of IgG+ B cells detected as VRC01-class",
        ),
    ]

    for i, (metric, ylabel, yname) in enumerate(metrics):
        ax = axes[i][0]
        sns.stripplot(
            data=df,
            x="method",
            y=metric,
            edgecolor="black",
            order=["Sanger", "10X"],
            linewidth=1,
            palette=pal,
            s=7,
            ax=ax,
        )
        sns.boxplot(
            data=df,
            x="method",
            y=metric,
            fliersize=0,
            order=["Sanger", "10X"],
            palette=pal,
            ax=ax,
        )
        adjust_boxplot(ax)
        ax.set_xticklabels([])
        ax.set_ylabel(ylabel)
        ax.set_xlabel("")
        df[["PTID", "trial", "method", metric]].to_csv(metric_outdir / f"{yname.replace(' ', '_')}.csv", index=False)

    ptids = df["PTID"]
    sequences_g001_10x = data.get_g001_10x_sequences()
    sequences_g001_10x["method"] = "10X"
    sequence_g002_sanger = data.get_g001_seqs()
    sequence_g002_sanger = sequence_g002_sanger[sequence_g002_sanger["PTID"].isin(ptids)]
    sequence_g002_sanger = sequence_g002_sanger.query("Timepoint=='V08'")
    sequence_g002_sanger["method"] = "Sanger"
    combined = pd.concat(
        [
            sequence_g002_sanger,
            sequences_g001_10x,
        ]
    ).reset_index(drop=True)

    combined = combined.query("is_vrc01_class == True")
    combined = combined.query("top_c_call=='IGHG'")

    seq_df = combined.groupby(["PTID", "method"])["v_mutation_aa_heavy"].median().to_frame().reset_index()
    ax = axes[3][0]
    sns.stripplot(
        data=seq_df,
        x="method",
        y="v_mutation_aa_heavy",
        edgecolor="black",
        order=["Sanger", "10X"],
        linewidth=1,
        s=7,
        palette=pal,
        ax=ax,
    )
    sns.boxplot(
        data=seq_df,
        x="method",
        y="v_mutation_aa_heavy",
        fliersize=0,
        order=["Sanger", "10X"],
        palette=pal,
        ax=ax,
    )
    adjust_boxplot(ax)
    ax.set_xticklabels(["G001 Method", "G002 Method"])
    label = r"$\mathregular{V_H}$" + f" gene\n%mutation (aa)"
    ax.set_ylabel(label)
    ax.set_xlabel("")

    yname = "figD V Heavy Gene Mutation AA".replace(" ", "_")
    seq_df["trial"] = seq_df["method"].replace({"Sanger": "G001", "10X": "G002", "KWTRP": "G003"})
    seq_df[
        [
            "PTID",
            "trial",
            "method",
            "v_mutation_aa_heavy",
        ]
    ].to_csv(metric_outdir / f"{yname}.csv", index=False)

    for ax, title in zip([row[0] for row in axes], ["A", "B", "C", "D"]):
        ax.set_title(title, fontsize=18, fontweight="bold", x=-0.1, y=1.1)

    plot_g003_methodology_comparision(
        data=data,
        axes=[row[1] for row in axes],
        sup_metric_outdir=metric_outdir,
        img_outdir=img_outdir,
    )

    for ax in [row[1] for row in axes]:
        ax.set_ylabel("")

    pal = {
        "G001": "#6C65FF",
        "G002": "#91FCC0",
    }
    plot_legend(
        pallete=pal,
        ax=axes[3][0],
        bbox_to_anchor=(0.5, -0.45),
        fontsize=14,
    )

    pal = {
        "G001": "#6C65FF",
        "G002": "#91FCC0",
        "G003": "#E377C2",
    }
    plot_legend(
        pallete=pal,
        ax=axes[3][1],
        bbox_to_anchor=(0.5, -0.45),
        fontsize=14,
    )

    for ax in faxes:
        sns.despine(ax=ax)

    # sns.despine()
    plt.tight_layout()

    for ax in faxes:
        ax.set_ylim(bottom=0.0)
        ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: f"{int(x)}" if x == 0 else round(x, 2)))

    for ax in faxes:
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        # ax.tick_params(axis="y", which="minor", length=4)

    g1, g2, g3, g4, g5 = [row[0] for row in axes]
    g6, g7, g8, g9, g10 = [row[1] for row in axes]

    c1 = pw.stack([g1, g2, g3, g4, g5], margin=0.1, operator="/")
    c2 = pw.stack([g6, g7, g8, g9, g10], margin=0.1, operator="/")

    g = pw.hstack(c1, c2, margin=0.1)

    g.savefig(img_outdir / "methodology_comparison.png", dpi=300)


def plot_g003_methodology_comparision2(data: Data, axes: List[plt.Axes], sup_metric_outdir: Path, img_outdir: Path):
    ptid_hcl = [
        # "PubID_079",
        "PubID_187",
        "PubID_153",
        "PubID_046",
    ]
    _, mg001_2_3_flow_and_seq = get_g00x_flow_and_seq(data)

    df = mg001_2_3_flow_and_seq
    pal = {
        "Sanger": "#6C65FF",
        "10X": "#91FCC0",
        "KWTRP": "#E377C2",
        "G001": "#6C65FF",
        "G002": "#91FCC0",
        "G003": "#E377C2",
    }
    # df = method_df.copy(deep=True).query('method!="NAC"')
    order = ["Sanger", "10X", "KWTRP"]
    # order = ["G001", "G002", "G003"]
    # df["trial"] = df["method"].replace(
    #     {"Sanger": "G001", "10X": "G002", "KWTRP": "G003"}
    # )

    metrics = [
        (
            # "percentage_of_antigen++",
            "Percent antigen-specific among IgG^{+}",
            "% GT8++\namong IgG+ B cells",
            "figF percent GT8++ among IgG+ B cells",
            df,
        ),
        (
            # "percentage_of_epitope_specific",
            "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}",
            "% CD4bs-specific\namong IgG+ B cells",
            "figG percent CD4bs-specific among IgG+ B cells",
            df,
        ),
        # (
        #     "Percent IgG^{+}KO^- among Ag^{++}",
        #     "% of GT8$^{++}$IgG$^{+}$\nB cells that are KO$^{-}$",
        #     "figI Percent of GT8++IgG+ B cells that are KO-",
        #     # df[df.pubID.isin(ptid_hcl)],
        #     df,
        # ),
        (
            # "Percent of IgG+ B cells detected as VRC01-class (missing seq to 0)",
            "Percent of VRC01-class sequences among IgG",
            "Percent of IgG+ B cells\ndetected as VRC01-class",
            "figH Percent of IgG+ B cells detected as VRC01-class",
            df[df.pubID.isin(ptid_hcl)],
        ),
    ]

    for i, (metric, ylabel, yname, df) in enumerate(metrics):
        ax = axes[i]
        sns.stripplot(
            data=df,
            x="method",
            y=metric,
            edgecolor="black",
            # order=order,
            linewidth=1,
            s=7,
            palette=pal,
            ax=ax,
        )
        sns.boxplot(
            data=df,
            x="method",
            y=metric,
            fliersize=0,
            # order=order,
            palette=pal,
            ax=ax,
        )
        adjust_boxplot(ax)
        ax.set_xticklabels([])
        ax.set_xlabel("")
        ax.set_ylabel(ylabel)
        df[["pubID", "trial", "method", metric]].to_csv(
            sup_metric_outdir / f"{yname.replace(' ', '_')}.csv", index=False
        )

    g003_seq_methods = get_g003_seq_methods(data=data).query('method!="NAC"')
    g003_seq_methods["method"] = "G003"
    ptids = df["PTID"]
    sequences_g001_10x = data.get_g001_10x_sequences()
    # sequences_g001_10x["PTID"] = sequences_g001_10x["pubID"]
    sequences_g001_10x["method"] = "G002"
    # sequence_g002_sanger = pd.read_feather(data.paths.g001_sequences)
    sequence_g002_sanger = data.get_g001_seqs()
    # sequence_g002_sanger["PTID"] = sequence_g002_sanger["pubID"]
    sequence_g002_sanger = sequence_g002_sanger[sequence_g002_sanger["pubID"].isin(ptid_hcl)]
    sequence_g002_sanger = sequence_g002_sanger.query("Timepoint=='V08'")
    sequence_g002_sanger["method"] = "G001"
    # breakpoint()
    combined = pd.concat(
        [
            sequence_g002_sanger,
            sequences_g001_10x,
            g003_seq_methods,
        ]
    ).reset_index(drop=True)
    combined = combined.query("is_vrc01_class == True")
    combined = combined.query("top_c_call=='IGHG'")

    seq_df = combined.groupby(["PTID", "method"])["v_mutation_aa_heavy"].median().to_frame().reset_index()
    seq_df = seq_df[seq_df.PTID.isin(ptid_hcl)]  # .query('method!="NAC"')

    ax = axes[3]
    sns.stripplot(
        data=seq_df,
        x="method",
        y="v_mutation_aa_heavy",
        edgecolor="black",
        # order=order,
        linewidth=1,
        s=7,
        palette=pal,
        ax=ax,
    )
    sns.boxplot(
        data=seq_df,
        x="method",
        y="v_mutation_aa_heavy",
        fliersize=0,
        # order=order,
        palette=pal,
        ax=ax,
    )
    adjust_boxplot(ax)
    ax.set_xticklabels([])
    ylabel = r"$\mathregular{V_H}$" + f" gene\n%mutation (aa)"
    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    yname = "figI V Heavy Gene Mutation AA".replace(" ", "_")
    seq_df["trial"] = seq_df["method"].replace({"Sanger": "G001", "10X": "G002", "KWTRP": "G003"})
    seq_df[
        [
            "PTID",
            "trial",
            "method",
            "v_mutation_aa_heavy",
        ]
    ].to_csv(sup_metric_outdir / f"{yname}.csv", index=False)

    seq_df = combined.groupby(["PTID", "method"])["v_mutation_aa_light"].median().to_frame().reset_index()
    seq_df = seq_df[seq_df.PTID.isin(ptid_hcl)]  # .query('method!="NAC"')
    ax = axes[4]
    sns.stripplot(
        data=seq_df,
        x="method",
        y="v_mutation_aa_light",
        edgecolor="black",
        # order=order,
        linewidth=1,
        s=7,
        palette=pal,
        ax=ax,
    )
    sns.boxplot(
        data=seq_df,
        x="method",
        y="v_mutation_aa_light",
        fliersize=0,
        # order=order,
        palette=pal,
        ax=ax,
    )
    adjust_boxplot(ax)
    ax.set_xticklabels(["G001 Method", "G002 Method", "G003 KWTRP"])
    ax.set_xlabel("")
    ylabel = r"$\mathregular{V_{K/L}}$" + " gene\n% mutation (aa)"
    yname = "figJ V Light Gene Mutation AA".replace(" ", "_")
    ax.set_ylabel(ylabel)
    seq_df["trial"] = seq_df["method"].replace({"Sanger": "G001", "10X": "G002", "KWTRP": "G003"})
    seq_df[
        [
            "PTID",
            "trial",
            "method",
            "v_mutation_aa_light",
        ]
    ].to_csv(sup_metric_outdir / f"{yname}.csv", index=False)

    for ax, title in zip(axes, ["F", "G", "H", "I", "J"]):
        ax.set_title(title, fontsize=16, fontweight="bold", x=-0.1, y=1.1)


def plot_methodology_comparision2(data: Data, metric_outdir: Path, img_outdir: Path):
    apply_global_font_settings(14)
    df, _ = get_g00x_flow_and_seq(data)
    df["trial"] = df["method"].replace({"Sanger": "G001", "10X": "G002", "KWTRP": "G003"})

    pal = {
        "Sanger": "#6C65FF",
        "10X": "#91FCC0",
        "KWTRP": "#E377C2",
        "G001": "#6C65FF",
        "G002": "#91FCC0",
        "G003": "#E377C2",
    }
    fig, axes = plt.subplots(4, 2, figsize=(8.5, 10))
    # change plot with patchworklib
    axes = [
        [pw.Brick(figsize=(4, 2)), pw.Brick(figsize=(4, 2))],
        [pw.Brick(figsize=(4, 2)), pw.Brick(figsize=(4, 2))],
        # [pw.Brick(figsize=(4, 2)), pw.Brick(figsize=(4, 2))],
        [pw.Brick(figsize=(4, 2)), pw.Brick(figsize=(4, 2))],
        [pw.Brick(figsize=(4, 2)), pw.Brick(figsize=(4, 2))],
        [pw.Brick(figsize=(4, 2)), pw.Brick(figsize=(4, 2))],
    ]
    faxes = [ax for row in axes for ax in row]

    metrics = [
        (
            # "percentage_of_antigen++",
            "Percent antigen-specific among IgG^{+}",
            "% GT8++\namong IgG+ B cells",
            "figA percent GT8++ among IgG+ B cells",
        ),
        (
            # "percentage_of_epitope_specific",
            "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}",
            "% CD4bs-specific\namong IgG+ B cells",
            "figB percent CD4bs-specific among IgG+ B cells",
        ),
        # (
        #     "Percent IgG^{+}KO^- among Ag^{++}",
        #     "% of GT8$^{++}$IgG$^{+}$\nB cells that are KO$^{-}$",
        #     "figC Percent of GT8++IgG+ B cells that are KO-",
        # ),
        (
            "Percent of VRC01-class sequences among IgG",
            "Percent of IgG+ B cells\ndetected as VRC01-class",
            "figC Percent of IgG+ B cells detected as VRC01-class",
        ),
    ]

    for i, (metric, ylabel, yname) in enumerate(metrics):
        ax = axes[i][0]
        sns.stripplot(
            data=df,
            x="method",
            y=metric,
            edgecolor="black",
            order=["G001", "G002"],
            linewidth=1,
            palette=pal,
            s=7,
            ax=ax,
        )
        sns.boxplot(
            data=df,
            x="method",
            y=metric,
            fliersize=0,
            order=["G001", "G002"],
            palette=pal,
            ax=ax,
        )
        adjust_boxplot(ax)
        ax.set_xticklabels([])
        ax.set_ylabel(ylabel)
        ax.set_xlabel("")
        df[["pubID", "trial", "method", metric]].to_csv(metric_outdir / f"{yname.replace(' ', '_')}.csv", index=False)
    # breakpoint()

    ptids = df["PTID"]
    sequences_g001_10x = data.get_g001_10x_sequences()
    sequences_g001_10x["method"] = "G002"
    sequence_g002_sanger = data.get_g001_seqs()
    sequence_g002_sanger = sequence_g002_sanger[sequence_g002_sanger["PTID"].isin(sequences_g001_10x.pubID.unique())]
    sequence_g002_sanger = sequence_g002_sanger.query("Timepoint=='V08'")
    sequence_g002_sanger["method"] = "G001"
    combined = pd.concat(
        [
            sequence_g002_sanger,
            sequences_g001_10x,
        ]
    ).reset_index(drop=True)
    combined = combined.query("is_vrc01_class == True")
    combined = combined.query("top_c_call=='IGHG'")
    seq_df = combined.groupby(["PTID", "method"])["v_mutation_aa_heavy"].median().to_frame().reset_index()

    ax = axes[3][0]
    sns.stripplot(
        data=seq_df,
        x="method",
        y="v_mutation_aa_heavy",
        edgecolor="black",
        # order=order,
        linewidth=1,
        s=7,
        palette=pal,
        ax=ax,
    )
    sns.boxplot(
        data=seq_df,
        x="method",
        y="v_mutation_aa_heavy",
        fliersize=0,
        # order=order,
        palette=pal,
        ax=ax,
    )
    adjust_boxplot(ax)
    ax.set_xticklabels([])
    ylabel = r"$\mathregular{V_H}$" + f" gene\n%mutation (aa)"
    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    yname = "figD V Heavy Gene Mutation AA".replace(" ", "_")
    seq_df["trial"] = seq_df["method"].replace({"Sanger": "G001", "10X": "G002", "KWTRP": "G003"})
    seq_df[
        [
            "PTID",
            "trial",
            "method",
            "v_mutation_aa_heavy",
        ]
    ].to_csv(metric_outdir / f"{yname}.csv", index=False)

    seq_df = combined.groupby(["PTID", "method"])["v_mutation_aa_light"].median().to_frame().reset_index()
    ax = axes[4][0]
    sns.stripplot(
        data=seq_df,
        x="method",
        y="v_mutation_aa_light",
        # order=["G001", "G002"],
        edgecolor="black",
        linewidth=1,
        s=7,
        palette=pal,
        ax=ax,
    )
    sns.boxplot(
        data=seq_df,
        x="method",
        y="v_mutation_aa_light",
        fliersize=0,
        # order=["G001", "G002"],
        palette=pal,
        ax=ax,
    )
    adjust_boxplot(ax)
    ax.set_xticklabels(["G001 Method", "G002 Method"])
    label = r"$\mathregular{V_{K/L}}$" + " gene\n% mutation (aa)"
    ax.set_ylabel(label)
    ax.set_xlabel("")

    yname = "figE V Light Gene Mutation AA".replace(" ", "_")
    seq_df["trial"] = seq_df["method"].replace({"Sanger": "G001", "10X": "G002", "KWTRP": "G003"})
    seq_df[
        [
            "PTID",
            "trial",
            "method",
            "v_mutation_aa_light",
        ]
    ].to_csv(metric_outdir / f"{yname}.csv", index=False)

    for ax, title in zip([row[0] for row in axes], ["A", "B", "C", "D", "E"]):
        ax.set_title(title, fontsize=18, fontweight="bold", x=-0.1, y=1.1)

    plot_g003_methodology_comparision2(
        data=data,
        axes=[row[1] for row in axes],
        sup_metric_outdir=metric_outdir,
        img_outdir=img_outdir,
    )

    for ax in [row[1] for row in axes]:
        ax.set_ylabel("")

    pal = {
        "G001": "#6C65FF",
        "G002": "#91FCC0",
    }
    plot_legend(
        pallete=pal,
        ax=axes[4][0],
        bbox_to_anchor=(0.5, -0.45),
        fontsize=14,
    )

    pal = {
        "G001": "#6C65FF",
        "G002": "#91FCC0",
        "G003": "#E377C2",
    }
    plot_legend(
        pallete=pal,
        ax=axes[4][1],
        bbox_to_anchor=(0.5, -0.45),
        fontsize=14,
    )

    for ax in faxes:
        sns.despine(ax=ax)

    # sns.despine()
    plt.tight_layout()

    for ax in faxes:
        ax.set_ylim(bottom=0.0)
        ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: f"{int(x)}" if x == 0 else round(x, 2)))

    for ax in faxes:
        ax.yaxis.set_minor_locator(mtick.AutoMinorLocator())
        # ax.tick_params(axis="y", which="minor", length=4)

    g1, g2, g4, g5, g6 = [row[0] for row in axes]
    g7, g8, g10, g11, g12 = [row[1] for row in axes]

    c1 = pw.stack([g1, g2, g4, g5, g6], margin=0.1, operator="/")
    c2 = pw.stack([g7, g8, g10, g11, g12], margin=0.1, operator="/")

    g = pw.hstack(c1, c2, margin=0.1)

    g.savefig(img_outdir / "methodology_comparison.png", dpi=300)
