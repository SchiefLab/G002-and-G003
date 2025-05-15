"""This is our main entry point"""

# %%
import functools
import logging
import math
import subprocess
from pathlib import Path

import click
import matplotlib.pyplot as plt
import numpy as np
import patchworklib as pw
import seaborn as sns
from matplotlib import ticker as mtick
from matplotlib.patches import Patch

from g00x_figures.alleles import plot_allele_table, plot_allele_table_group4
from g00x_figures.box_and_scatter.flow_frequencies import (  # plot_flow_frequencies,
    MinorSymLogLocator,
    adjust_boxplot,
    plot_boost_frequences,
    plot_cp_frequency,
    plot_group_4_responders,
    plot_pre_post_frequency,
    plot_response_summary,
)
from g00x_figures.counts import FlowCytometryPlot
from g00x_figures.data import Data, Transforms
from g00x_figures.features import (  # plot_key_mutations_boost,
    plot_has_100b,
    plot_has_100b_boost,
    plot_isotype_data,
    plot_isotype_data_pseudogroups,
    plot_key_mutations,
    plot_light_chain_dist,
    plot_light_chain_dist_boost,
    plot_light_chain_usage,
    plot_light_chain_usage_boost,
    plot_qe_on_light_chain,
    plot_qe_on_light_chain_boost,
    plot_seq_logos_nonvrc01_class,
    plot_seq_logos_vrc01_class,
)
from g00x_figures.flow_frequencies import plot_flow_frequencies
from g00x_figures.misc import (
    plot_methodology_comparision,
    plot_methodology_comparision2,
)
from g00x_figures.mutations import (
    plot_key_mutations_boost,
    plot_somatic_mutation_aa_frequencies,
    plot_somatic_mutation_aa_frequencies_boost,
    plot_somatic_mutation_nt_frequencies,
    plot_v_mutations,
    run_90_percentile_hc_residues,
)
from g00x_figures.plot_helpers.font import apply_global_font_settings
from g00x_figures.plots import Plot
from g00x_figures.polyclonality import (
    plot_boost_clonality,
    plot_multi_find_clonality,
    plot_polyclonality,
)
from g00x_figures.properties import spr_properties
from g00x_figures.qc import plot_qc
from g00x_figures.responders import run_percent_igg_responders
from g00x_figures.spr import (
    plot_core_spr,
    plot_core_spr_kon_off,
    plot_shm_spr,
    plot_spr_core_candidates,
    plot_spr_prime,
)


@click.group(invoke_without_command=True)
@click.pass_context
@click.option(
    "--outdir",
    "-o",
    default=None,
    type=click.Path(file_okay=True, resolve_path=True),
)
@click.option(
    "--fig",
    "-f",
    default=None,
    type=str,
    help="Figure number",
    show_default=True,
)
@click.option(
    "--is-main",
    "-m",
    is_flag=True,
)
@click.option(
    "--use-geomean",
    "-g",
    is_flag=True,
    help="Add geometric mean to header and metric files",
)
@click.option(
    "--median-scale",
    "-s",
    type=float,
    help="Scale factor for median values",
    default=1e9,
)
def cli(
    ctx: click.Context,
    outdir: Path,
    fig: str,
    is_main: bool,
    use_geomean,
    median_scale,
) -> None:
    """Plot figures."""
    ctx.obj = {}
    ctx.obj = {
        "data": Data(),
        "transforms": Transforms(),
        "output": Path(__file__).parent.parent.parent,
    }

    data = ctx.obj["data"]
    ctx.obj["outdir"] = outdir if outdir else data.paths.figure_outdir

    if is_main:
        dest = "Main"
    else:
        dest = "Sup"

    # IMG default to file specific: /figN/ not needed for most figures
    ctx.obj["img_outdir"] = ctx.obj["outdir"] / f"{dest}"
    # METRIC default to folder specific
    ctx.obj["metric_outdir"] = ctx.obj["outdir"] / f"{dest}-Metrics/{fig}"
    if fig is not None:
        ctx.obj["img_outdir"].mkdir(parents=True, exist_ok=True)
        ctx.obj["metric_outdir"].mkdir(parents=True, exist_ok=True)

    ctx.obj["fig"] = fig
    ctx.obj["is_main"] = is_main
    ctx.obj["use_geomean"] = use_geomean
    ctx.obj["median_scale"] = median_scale

    for key in ["fig", "img_outdir", "metric_outdir"]:
        value = ctx.obj[key]
        if fig is not None:
            logging.info(f"{key}: {value}")


@cli.command("flow-freq")
@click.pass_context
def fig2(ctx: click.Context) -> None:
    """Flow frequencies."""
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    use_alt = False
    img_outdir = Path(outdir) / f"{fig}"
    metric_outdir = Path(outdir) / "metrics"
    img_outdir.mkdir(parents=True, exist_ok=True)
    metric_outdir.mkdir(parents=True, exist_ok=True)

    plot_flow_frequencies(
        data=data,
        img_outdir=img_outdir,
        metric_outdir=metric_outdir,
        use_alt=use_alt,
    )


# @figures.command("fig8")
@cli.command("fig8")
@click.pass_context
def fig8(ctx: click.Context, *args, **kwargs) -> None:
    breakpoint()
    data = ctx.obj["data"]
    img_outdir = ctx.obj["img_outdir"]
    img_outpath = img_outdir / f"{ctx.obj['fig']}.png"
    metric_outdir = ctx.obj["metric_outdir"]
    breakpoint()
    plot_cp_frequency(data, img_outpath=img_outpath, metric_outdir=metric_outdir)


@cli.command("prime-mut")
@click.pass_context
@click.option(
    "--aa",
    "-a",
    is_flag=True,
)
@click.option(
    "--method",
    type=str,
    help="Method to np.qunatile",
    default="nearest",
)
@click.option(
    "--no-panel-e",
    is_flag=True,
    help="Do not plot panel E",
)
def prime_mutations(ctx: click.Context, aa: bool, method: str, no_panel_e: bool) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outdir = Path(outdir)
    metric_outdir = Path(outdir) / "metrics"

    name = ctx.command.name

    logging.info(f"Figure {fig}")
    logging.info(f"Plotting {name}")
    logging.info(f"Output directory: {outdir}")

    (g1, g2, g3, g4) = plot_v_mutations(
        data=data,
        src_type="prime",
        seq_type="aa" if aa else "nt",
        mectric_outdir=metric_outdir,
    )
    c1 = pw.stack([g1, g3], operator="/", margin=0)
    c2 = pw.stack([g2, g4], operator="/", margin=0)
    b1 = c1 | c2

    if no_panel_e:
        g = b1
    else:
        g5, df = run_90_percentile_hc_residues(data=data, metric_outdir=metric_outdir, method=method)
        figE_medians = df.groupby(["trial", "pseudogroup", "weeks"])[f"residues"].median()
        figE_medians.to_csv(metric_outdir / "figE_medians.csv", index=True)

        b2 = pw.spacer(g5, 0.05) | g5 | pw.spacer(g5, 0.15)
        g = b1 / b2
    g.savefig(img_outdir / f"{fig}.png", dpi=700)


@cli.command("boost-freq")
@click.pass_context
def boost_freq(ctx: click.Context) -> None:
    """Figure 4."""
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    outpath = Path(outdir)
    logging.info(f"Figure {fig}")
    # TODO: hard coded for now in path.
    plot_boost_frequences(data=data, outpath=outpath, fig_num=fig)


@cli.command("boost-mut-aa")
@click.pass_context
@click.option(
    "--method",
    type=str,
    help="Method to np.qunatile",
    default="nearest",
)
def fig5(ctx: click.Context, method: str) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    is_main = ctx.obj["is_main"]
    outdir = ctx.obj["outdir"]
    logging.info(f"Figure {fig}")
    if is_main:
        dest = "Main"
    else:
        dest = "Sup"
    img_outdir = Path(outdir) / f"{dest}"
    metric_outdir = Path(outdir) / f"{dest}-Metrics/{fig}"
    img_outdir.mkdir(parents=True, exist_ok=True)
    metric_outdir.mkdir(parents=True, exist_ok=True)

    (g1, g2, g3, g4) = plot_v_mutations(
        data=data,
        src_type="boost",
        seq_type="aa",
        mectric_outdir=metric_outdir,
    )
    g5, g6, g7, g8 = plot_key_mutations_boost(data=data, metric_outdir=metric_outdir, method=method)
    g = (pw.stack([g1, g3], operator="/", margin=0) | pw.stack([g2, g4], operator="/", margin=0)) / (
        pw.stack([g5, g6], margin=0) | pw.stack([g7, g8], margin=0)
    )
    g.savefig(img_outdir / f"{fig}.png", dpi=700)


@cli.command("light-chain")
@click.pass_context
def light_chain_quality_control(ctx: click.Context) -> None:
    """Figure S1."""
    data = ctx.obj["data"]


@cli.command("methodology")
@click.pass_context
def methodology_comparison(ctx: click.Context) -> None:
    """Methodology comparison btw G001 (Sanger) and G002/G003 (10X)"""
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]

    outdir = Path(outdir)
    img_outdir = outdir / f"Sup/{fig}"
    metric_outdir = outdir / f"Sup-Metrics/{fig}"
    plot_methodology_comparision(data=data, img_outdir=img_outdir, metric_outdir=metric_outdir)


@cli.command("methodology2")
@click.pass_context
def methodology_comparison(ctx: click.Context) -> None:
    """Methodology comparison btw G001 (Sanger) and G002/G003 (10X)"""
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]

    outdir = Path(outdir)
    img_outdir = outdir / f"Sup/{fig}"
    metric_outdir = outdir / f"Sup-Metrics/{fig}"
    img_outdir.mkdir(parents=True, exist_ok=True)
    metric_outdir.mkdir(parents=True, exist_ok=True)
    plot_methodology_comparision2(data=data, img_outdir=img_outdir, metric_outdir=metric_outdir)


@cli.command("boost-clonality")
@click.pass_context
def boost_clonality(ctx: click.Context) -> None:
    """Plot boost clonality for supplementary figure 29."""
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    outdir = Path(outdir)
    img_outdir = outdir / "images"
    metric_outdir = outdir / "metrics"
    name = ctx.command.name
    img_outpath = img_outdir / f"{fig}-{name}.png"
    logging.info(f"Plotting {name}: {img_outpath}")
    g1 = plot_boost_clonality(
        data=data,
        img_outdir=outdir,
        metric_outdir=metric_outdir,
        fig="figA",
        y_axis="clonality",
    )
    g2 = plot_boost_clonality(
        data=data,
        img_outdir=outdir,
        metric_outdir=metric_outdir,
        fig="figB",
        y_axis="total_unique",
    )
    g3 = plot_boost_clonality(
        data=data,
        img_outdir=outdir,
        metric_outdir=metric_outdir,
        fig="figC",
        y_axis="median_cluster_size",
    )
    g = g1 / g2 / g3
    g.savefig(img_outpath, dpi=700)


@cli.command("number-of-BCR-clusters")
@click.pass_context
@click.option(
    "--use-last-week",
    "-u",
    is_flag=True,
)
def num_bcr_clusters(ctx: click.Context, use_last_week: bool) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    outdir = Path(outdir)
    img_outdir = outdir / "images"
    logging.info("Plotting: Number of BCR clusters")
    logging.info(f"Output directory: {outdir}")
    name = ctx.command.name
    img_outpath = img_outdir / f"{fig}-{name}.png"
    df = plot_multi_find_clonality(data, img_outpath, use_last_week=use_last_week)
    df.to_csv(outdir / "metrics" / f"{fig}-{name}.csv", index=False)


@cli.command("hierarchical-clustering-and-polyclonality")
@click.pass_context
@click.option(
    "--use-last-week",
    "-u",
    is_flag=True,
)
def hierarchical_clustering_and_polyclonality(ctx: click.Context, use_last_week: bool) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    outdir = Path(outdir)
    img_outdir = outdir / "images"
    logging.info("Plotting: Number of BCR clusters")
    logging.info(f"Output directory: {outdir}")
    name = ctx.command.name
    img_outpath = img_outdir / f"{fig}-{name}.png"
    metric_outdir = outdir / "metrics"
    plot_polyclonality(
        data,
        img_outpath,
        metric_outdir=metric_outdir,
        use_last_week=use_last_week,
    )


@cli.command("b-count-spr-g00x-eOD")
@click.pass_context
def s28(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    metric_outdir = Path(outdir) / "metrics"
    img_outpath = Path(outdir) / f"{fig}-spr-eOD.png"
    g1, vrc01_dfs = plot_spr_prime(
        data=data,
        analyte="eOD-GT8.1_His-Avi_mC",
        is_vrc01_class=True,
        igl_pad=9,
    )
    g2, nonvrc01_dfs = plot_spr_prime(
        data=data,
        analyte="eOD-GT8.1_His-Avi_mC",
        is_vrc01_class=False,
        letter="B",
        igl_pad=8,
    )
    g = pw.stack([g1, g2], operator="/", margin=0)
    save(g, vrc01_dfs + nonvrc01_dfs, img_outpath, metric_outdir)


@cli.command("b-count-spr-g00x-core")
@click.pass_context
def s33(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    metric_outdir = Path(outdir) / "metrics"
    img_outpath = Path(outdir) / f"{fig}-spr-core.png"
    g1, vrc01_dfs = plot_spr_prime(
        data=data,
        analyte="core-Hx_r4.0D_TH6_g28v2_pHLsecAvi",
        is_vrc01_class=True,
        igl_pad=8,
        remove_igl=True,
        mature_pad=7,
    )
    g2, nonvrc01_dfs = plot_spr_prime(
        data=data,
        analyte="core-Hx_r4.0D_TH6_g28v2_pHLsecAvi",
        is_vrc01_class=False,
        letter="B",
        remove_igl=True,
        mature_pad=7,
    )
    g = pw.stack([g1, g2], operator="/", margin=0)
    save(g, vrc01_dfs + nonvrc01_dfs, img_outpath, metric_outdir)


@cli.command("fig7")
@click.pass_context
def fig7(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outpath = Path(outdir) / f"{fig}.png"
    metric_outdir = Path(outdir) / "metrics"
    plot_spr_core_candidates(
        data=data,
        outpath=img_outpath,
        metric_outdir=metric_outdir,
        use_geomean=ctx.obj["use_geomean"],
        median_scale=ctx.obj["median_scale"],
    )


@cli.command("fig38")
@click.pass_context
def fig38(ctx: click.Context) -> None:
    """Figure 7."""
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outpath = Path(outdir) / f"{fig}.png"
    analyte_mapping = {
        "191084_SOSIP_MD39_mC": "191084",
        "1HD2_B4_S62_MD39v2_L14_mC2": "1HD2",
        "001428_MD39_L14_m2": "001428",
        "BG505_SOSIP_MD39_m": "BG505",
        "235_47_RnS_2G_L14_T278M_mC2": "235-T278M",
        "BG505_MD39.3_cd4bsHxB2_M278_mC2": "BG505-cd4bsHxB2-T278M",
    }
    metric_outdir = Path(outdir) / "metrics"
    plot_spr_core_candidates(
        data=data,
        outpath=img_outpath,
        is_cp=True,
        analyte_mapping=analyte_mapping,
        metric_outdir=metric_outdir,
        use_geomean=ctx.obj["use_geomean"],
        median_scale=ctx.obj["median_scale"],
    )


@cli.command("sup31")
@click.pass_context
def sup31(ctx: click.Context) -> None:
    """Figure 1."""
    logging.info("Figure S31")
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    metric_outdir = Path(outdir) / "metrics"
    img_outdir = Path(outdir) / "images"
    img_path = img_outdir / f"{fig}/{fig}"

    tdf, df = plot_light_chain_usage_boost(data, str(img_path) + "A")
    df = data.populate_psname(df.reset_index()).fillna(0)
    df.to_csv(metric_outdir / "figA.csv", index=False)
    tdf.to_csv(metric_outdir / "figA_seqs.csv", index=False)

    tdf, df = plot_qe_on_light_chain_boost(data, str(img_path) + "B")
    df = data.populate_psname(df)
    df.to_csv(metric_outdir / "figC.csv", index=False)
    tdf.to_csv(metric_outdir / "figC_seqs.csv", index=False)

    ktdf, kappa_df, ltdf, lambda_df = plot_light_chain_dist_boost(data, str(img_path) + "C")
    kappa_df = data.populate_psname(kappa_df)
    lambda_df = data.populate_psname(lambda_df)
    kappa_df.to_csv(metric_outdir / "figB_kappa.csv", index=False)
    lambda_df.to_csv(metric_outdir / "figB_lambda.csv", index=False)
    ktdf.to_csv(metric_outdir / "figB_kappa_seqs.csv", index=False)
    ltdf.to_csv(metric_outdir / "figB_lambda_seqs.csv", index=False)

    tdf, df = plot_has_100b_boost(data, str(img_path) + "D")
    df = data.populate_psname(df)
    df.to_csv(metric_outdir / "figD.csv", index=False)
    tdf.to_csv(metric_outdir / "figD_seqs.csv", index=False)


@cli.command("fig8")
@click.pass_context
def fig8(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outdir = Path(outdir)
    img_outpath = img_outdir / f"{fig}.png"
    metric_outdir = Path(outdir) / "metrics"
    plot_cp_frequency(data, img_outpath=img_outpath, metric_outdir=metric_outdir)


@cli.command("spr-boost")
@click.pass_context
def spr_boost(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outdir = Path(outdir)
    img_outpath = img_outdir / f"{fig}.png"
    metric_outdir = Path(outdir) / "metrics"
    plot_core_spr(
        data,
        outpath=img_outpath,
        metric_outdir=metric_outdir,
        use_geomean=ctx.obj["use_geomean"],
        median_scale=ctx.obj["median_scale"],
    )


def save(img, dfs, img_outpath, metric_outdir):
    img.savefig(img_outpath, dpi=700)
    key_list = ["trial", "pseudogroup", "weeks", "pubID"]
    for df in dfs:
        name = df.name
        kl = [k for k in key_list if k in df.columns]
        key = [k for k in df.key if k not in kl]
        df = df[kl + key]
        df = df.sort_values(kl)
        df.to_csv(metric_outdir / f"fig_{name}.csv", index=False)


@cli.command("b-count-flow-g002")
@click.pass_context
def s20(ctx: click.Context) -> None:
    flow_cytometry_plot = FlowCytometryPlot()
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outdir = Path(outdir)
    img_outpath = img_outdir / f"{fig}.png"
    metric_outdir = Path(outdir) / "metrics"
    img, dfs = flow_cytometry_plot.create_plot(
        x="group",
        df=data.get_g002_flow_and_seq_prime(),
        palette={
            -5: "#9567BD",
            4: "#17BFD0",
            8: "gold",
            16: "#E377C2",
            24: "#2078B4",
        },
    )
    save(img, dfs, img_outpath, metric_outdir)


@cli.command("b-count-boost-g002")
@click.pass_context
def s31(ctx: click.Context) -> None:
    flow_cytometry_plot = FlowCytometryPlot()
    flow_cytometry_plot.columns = {
        "B cells": "B cells",
        "IgD-IgM-/IgA-/IgG+ B cells": "IgG$^{+}$ B cells",
        "IgD-/IgM-/IgA-/IgG+/Antigen++ B cells": "core$^{++}$ IgG$^{+}$ B cells",
        "IgD-/IgM-/IgA-/IgG+/KO-/Antigen++ B cells": "core CD4bs-specific\nIgG$^{+}$ B cells",
    }
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outdir = Path(outdir)
    img_outpath = img_outdir / f"{fig}.png"
    metric_outdir = Path(outdir) / "metrics"
    img, dfs = flow_cytometry_plot.create_plot(
        x="group",
        df=data.get_g002_flow_and_seq_boost(),
        palette={
            -5: "#9567BD",
            4: "#17BFD0",
            8: "gold",
            16: "#E377C2",
            24: "#2078B4",
        },
    )
    save(img, dfs, img_outpath, metric_outdir)


@cli.command("b-count-flow-g003")
@click.pass_context
def s21(ctx: click.Context) -> None:
    flow_cytometry_plot = FlowCytometryPlot()
    flow_cytometry_plot.columns = {
        "IgD- B cells": "IgD$^{-}$ B cells",
        "IgD-IgM-IgG+ B cells": "IgG$^{+}$ B cells",
        "IgD-/IgM-/IgG+/Antigen++ B cells": "eOD$^{++}$ IgG$^{+}$ B cells",
        "IgM-/IgG+/KO-/Antigen++ B cells": "eOD CD4bs-specific\nIgG$^{+}$ B cells",
    }
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outdir = Path(outdir)
    img_outpath = img_outdir / f"{fig}.png"
    metric_outdir = Path(outdir) / "metrics"
    img, dfs = flow_cytometry_plot.create_plot(
        x="weeks",
        df=data.get_g003_flow_and_seq_prime(),
        palette={
            -5: "#9567BD",
            8: "#17BFD0",
            10: "gold",
            16: "#E377C2",
            21: "#2078B4",
        },
    )
    for ax in img.get_axes():
        ax.set_title("")
        ax.set_xlabel("")
        if not ax.get_xticklabels():
            continue
        print(ax.get_xticklabels())
        ax.set_xticklabels([f"Week {tl.get_text()}" for tl in ax.get_xticklabels()], rotation=0)
    save(img, dfs, img_outpath, metric_outdir)


@cli.command("s-freq-seq-g002")
@click.pass_context
def s22(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outdir = Path(outdir)
    img_outpath = img_outdir / f"{fig}.png"
    metric_outdir = Path(outdir) / "metrics"
    img, dfs = plot_isotype_data(seq=data.get_g002_sequences_prime())
    save(img, dfs, img_outpath, metric_outdir)


@cli.command("s-freq-seq-g003")
@click.pass_context
def s23(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outdir = Path(outdir)
    img_outpath = img_outdir / f"{fig}.png"
    metric_outdir = Path(outdir) / "metrics"
    img, dfs = plot_isotype_data(seq=data.get_g003_sequences_prime())
    save(img, dfs, img_outpath, metric_outdir)


@cli.command("b-count-spr-g002")
@click.pass_context
def s28(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outdir = Path(outdir)
    img_outpath = img_outdir / f"{fig}.png"
    metric_outdir = Path(outdir) / "metrics"
    img, dfs = plot_spr_core_candidates(
        data=data,
        outpath=img_outpath,
        is_cp=True,
        is_vrc01_class=True,
        metric_outdir=metric_outdir,
    )
    save(img, dfs, img_outpath, metric_outdir)


@cli.command("b-freq-prepost-g002")
@click.pass_context
def s32(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outdir = Path(outdir)
    img_outpath = img_outdir / f"{fig}.png"
    metric_outdir = Path(outdir) / "metrics"
    img, dfs = plot_pre_post_frequency(data=data)
    save(img, dfs, img_outpath, metric_outdir)


@cli.command("s-freq-boost-g002")
@click.pass_context
def s34(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outdir = Path(outdir)
    img_outpath = img_outdir / f"{fig}.png"
    metric_outdir = Path(outdir) / "metrics"
    img, dfs = plot_isotype_data_pseudogroups(seq_boost=data.get_g002_sequences_boost())
    save(img, dfs, img_outpath, metric_outdir)


@cli.command("b-shm-spr-g002")
@click.pass_context
def s40(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outdir = Path(outdir)
    img_outpath = img_outdir / f"{fig}.png"
    metric_outdir = Path(outdir) / "metrics"
    img, dfs = plot_shm_spr(data=data)
    save(img, dfs, img_outpath, metric_outdir)


@cli.command("s-properties-spr-g002")
@click.pass_context
def s41(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outdir = Path(outdir)
    img_outpath = img_outdir / f"{fig}.png"
    metric_outdir = Path(outdir) / "metrics"
    img, dfs = spr_properties(data=data)
    save(img, dfs, img_outpath, metric_outdir)


@cli.command("b-kon-koff-spr-g28v2-g002")
@click.pass_context
def s42(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    analyte_ln = "core-Hx_r4.0D_TH6_g28v2_pHLsecAvi"
    analyte_sn = "core-g28v2"
    img_outpath = Path(outdir) / f"{fig}.png"
    img, dfs = plot_core_spr_kon_off(
        data,
        analyte_long_name=analyte_ln,
        analyte_short_name=analyte_sn,
        lt=5e-5,
    )
    save(img, dfs, img_outpath, Path(outdir) / "metrics")


@cli.command("b-kon-koff-spr-N276-g002")
@click.pass_context
def s45(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    analyte_ln = "core-Hx_r4.0_TH6_g28v2_N276_pHLsecAvi"
    analyte_sn = "core-N276"
    img_outpath = Path(outdir) / f"{fig}.png"
    img, dfs = plot_core_spr_kon_off(
        data,
        analyte_long_name=analyte_ln,
        analyte_short_name=analyte_sn,
        lt=1e-4,
    )
    save(img, dfs, img_outpath, Path(outdir) / "metrics")


@cli.command("b-freq-spr-g002-nonvrc01")
@click.pass_context
def sup51(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outpath = Path(outdir) / f"{fig}.png"
    analyte_mapping = {
        "191084_SOSIP_MD39_N276D_mC": "191084-N276D",
        "1HD2_B4_S62_MD39v2_L14_N276Q_mC2": "1HD2-N276Q",
        "001428_MD39_L14_T278M_m2": "001428-T278M",
        "BG505_SOSIP_MD39_N276Q_mC": "BG505-N276Q",
        "235_47_RnS_2G_L14_T278M_mC2": "235-T278M",
        "BG505_MD39.3_cd4bsHxB2_M278_mC2": "BG505-cd4bsHxB2-T278M",
    }
    metric_outdir = Path(outdir) / "metrics"
    plot_spr_core_candidates(
        data=data,
        outpath=img_outpath,
        is_cp=False,
        is_vrc01_class=False,
        ylabel="$\mathregular{K_D}$\nnon VRC01-class",
        analyte_mapping=analyte_mapping,
        metric_outdir=metric_outdir,
        median_scale=ctx.obj["median_scale"],
    )


@cli.command("S24")
@click.pass_context
def s24(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outdir = Path(outdir)
    img_outpath = img_outdir / f"{fig}.png"
    metric_outdir = Path(outdir) / "metrics"
    img, dfs = run_percent_igg_responders(data)
    save(img, dfs, img_outpath, metric_outdir)


@cli.command("S27")
@click.pass_context
def s27(ctx: click.Context) -> None:
    data = ctx.obj["data"]
    fig = ctx.obj["fig"]
    outdir = ctx.obj["outdir"]
    img_outdir = Path(outdir)
    img_outpath = img_outdir / f"{fig}.png"
    metric_outdir = Path(outdir) / "metrics"
    img, dfs = plot_qc()
    save(img, dfs, img_outpath, metric_outdir)


@cli.command("test")
@click.pass_context
def plot_test(ctx: click.Context) -> None:
    """Plot all figures."""
    commands = [
        "echo 'Plotting main figures'",
    ]
    click.echo(ctx.obj)
    for cmd in commands:
        subprocess.run(cmd, shell=True, check=True)


@cli.command("main")
@click.pass_context
def plot_main(ctx: click.Context) -> None:
    """Plot all figures."""
    commands = [
        "g00x plot -m --fig fig2 flow-freq",
        "g00x plot -m --fig fig3-nearest prime-mut --aa --method nearest",
        "g00x plot -m --fig fig3-midpoint prime-mut --aa --method midpoint",
        "g00x plot -m --fig fig4 boost-freq",
        "g00x plot -m --fig fig5 boost-mut-aa --method nearest",
        "g00x plot -m --fig fig6 -g -s 1 spr-boost",
        "g00x plot -m --fig fig7 fig8",
    ]
    for cmd in commands:
        subprocess.run(cmd, shell=True, check=True)


@cli.command("supp")
@click.pass_context
def plot_sup(ctx: click.Context) -> None:
    """Plot all figures."""
    commands = [
        "g00x plot --fig S19 methodology2",
        "g00x plot --fig S20 b-count-flow-g002",
        "g00x plot --fig S21 b-count-flow-g003",
        "g00x plot --fig S22 s-freq-seq-g002",
        "g00x plot --fig S23 s-freq-seq-g003",
        "g00x plot --fig S24 S24",
        "g00x plot --fig S25 prime-mut --method nearest --no-panel-e",
        "g00x plot --fig S27 S27",
        "g00x plot --fig S28 b-count-spr-g00x-eOD",
        "g00x plot --fig S29 number-of-BCR-clusters -u",
        "g00x plot --fig S30 hierarchical-clustering-and-polyclonality -u",
        "g00x plot --fig S31 b-count-boost-g002",
        "g00x plot --fig S32 b-freq-prepost-g002",
        "g00x plot --fig S33 b-count-spr-g00x-core",
        "g00x plot --fig S34 s-freq-boost-g002",
        "g00x plot --fig S35 boost-clonality",
        "g00x plot --fig S40 b-shm-spr-g002",
        "g00x plot --fig S41 s-properties-spr-g002",
        "g00x plot --fig S42 b-kon-koff-spr-g28v2-g002",
        "g00x plot --fig S45 b-kon-koff-spr-N276-g002",
        "g00x plot --fig S46 -g -s 1 fig38",
        "g00x plot --fig S47 -g -s 1 fig7",
        "g00x plot --fig S51 -g -s 1 b-freq-spr-g002-nonvrc01",
    ]
    assert len(set([c.split(" ")[-1] for c in commands])) == len(commands)
    for cmd in commands:
        subprocess.run(cmd, shell=True, check=True)


@cli.command("quant")
@click.pass_context
def plot_sup(ctx: click.Context) -> None:
    """Plot all figures."""
    commands = [
        "g00x plot -m --fig fig3 prime-mut --aa --method nearest",
        "g00x plot -m --fig fig5 boost-mut-aa --method nearest",
    ]
    for cmd in commands:
        subprocess.run(cmd, shell=True, check=True)


@cli.command("spr")
@click.pass_context
def plot_sup(ctx: click.Context) -> None:
    """Plot all figures."""
    commands = [
        "g00x plot -m --fig fig6 -g -s 1 spr-boost",
        "g00x plot --fig S28 b-count-spr-g00x-eOD",
        "g00x plot --fig S33 b-count-spr-g00x-core",
        "g00x plot --fig S42 b-kon-koff-spr-g28v2-g002",
        "g00x plot --fig S45 b-kon-koff-spr-N276-g002",
        "g00x plot --fig S46 -g -s 1 fig38",
        "g00x plot --fig S47 -g -s 1 fig7",
        "g00x plot --fig S51 -g -s 1 b-freq-spr-g002-nonvrc01",
    ]
    assert len(set([c.split(" ")[-1] for c in commands])) == len(commands)
    for cmd in commands:
        subprocess.run(cmd, shell=True, check=True)
