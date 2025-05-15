"""This is our main entry point"""
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
from g00x_figures.box_and_scatter.flow_frequencies import (
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
from g00x_figures.features import (
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

# @click.group()
# @click.pass_context
# def cli(ctx):
#     """
#     Figures.
#     """
#     sns.set_context("paper", font_scale=1.2)  # type: ignore
#     sns.set_style("ticks")
#     sns.set_style({"font.family": "Arial"})
#     np.random.seed(1000)
#     ctx.obj = {}
#     ctx.obj = {
#         "data": Data(),
#         "transforms": Transforms(),
#         "output": Path(__file__).parent.parent.parent,
#     }
#     root_logger = logging.getLogger()
#     root_logger.setLevel(logging.INFO)


@click.group("plot", invoke_without_command=True)
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
    required=True,
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
def figures(
    ctx: click.Context,
    outdir: Path,
    fig: str,
    is_main: bool,
    use_geomean,
    median_scale,
) -> None:
    """Plot figures."""
    sns.set_context("paper", font_scale=1.2)  # type: ignore
    sns.set_style("ticks")
    sns.set_style({"font.family": "Arial"})
    np.random.seed(1000)
    ctx.obj = {}
    ctx.obj = {
        "data": Data(),
        "transforms": Transforms(),
        "output": Path(__file__).parent.parent.parent,
    }
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    data = ctx.obj["data"]
    ctx.obj["outdir"] = Path("../G00X-plots-test")  # outdir if outdir else data.paths.figure_outdir

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


@figures.command("test")
@click.pass_context
@click.option(
    "--use-alt",
    "-a",
    is_flag=True,
    help="Alternate version of figure",
)
def plot_test(ctx: click.Context, *args, **kwargs) -> None:
    """Plot all figures."""
    commands = [
        "echo 'Plotting main figures'",
    ]
    click.echo(ctx.obj)
    click.echo(args)
    click.echo(kwargs)
    for cmd in commands:
        subprocess.run(cmd, shell=True, check=True)


@figures.command("fig2")
@click.pass_context
def fig2(ctx: click.Context) -> None:
    """Flow frequencies."""
    fig = ctx.obj["fig"]
    data = ctx.obj["data"]
    img_outdir = ctx.obj["img_outdir"] / f"{fig}"
    metric_outdir = ctx.obj["metric_outdir"]
    img_outdir.mkdir(parents=True, exist_ok=True)
    metric_outdir.mkdir(parents=True, exist_ok=True)
    click.echo(f"Plotting figure {fig} to {img_outdir}")
    click.echo(f"Saving metrics to {metric_outdir}")
    # plot_flow_frequencies(
    #     data=data,
    #     img_outdir=img_outdir,
    #     metric_outdir=metric_outdir,
    #     use_alt=use_alt,
    # )


@figures.command("main")
@click.pass_context
def plot_main(ctx: click.Context) -> None:
    """Plot all figures."""
    commands = [
        "g00x plot -m --fig fig2 fig2",
        # "g00x plot -m --fig fig3-nearest prime-mut --aa --method nearest",
        # "g00x plot -m --fig fig3-midpoint prime-mut --aa --method midpoint",
        # "g00x plot -m --fig fig4 boost-freq",
        # "g00x plot -m --fig fig5 boost-mut-aa --method nearest",
        # "g00x plot -m --fig fig6 -g -s 1 spr-boost",
        # "g00x plot -m --fig fig7 fig8",
    ]
    for cmd in commands:
        subprocess.run(cmd, shell=True, check=True)
