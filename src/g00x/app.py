"""This is our main entry point"""
import logging
import subprocess
import warnings
from pathlib import Path

import click
import pandas as pd
from pandas.errors import PerformanceWarning

from g00x.analysis.flow import count_current_samples
from g00x.analysis.g003_report import g003_combine_seq_and_flow
from g00x.analysis.report import combine_seq_and_flow
from g00x.data import Data, PlotParameters
from g00x.flow import g003_flow
from g00x.flow.flow import parse_flow_data
from g00x.sequencing.airr import run_airr
from g00x.sequencing.g003_airr import g003_run_airr
from g00x.sequencing.g003_tenX import g003_run_cso, g003_run_demultiplex, g003_run_vdj
from g00x.sequencing.merge import merge_flow_and_sequencing
from g00x.sequencing.tenX import run_cso, run_demultiplex, run_vdj
from g00x.tools.path import cd, pathing, pd_expand_path, pd_replace_home_with_tilde
from g00x.validations.flow_validation import ValidateG00X
from g00x.validations.g003_flow_validation import validate_g003_sorting
from g00x.validations.g003_sequencing_validation import validate_g003_sequencing
from g00x_figures.app import cli as plot_cli


@click.group("g00x")
@click.option("--logging-level", default="INFO", help="Set logging level")
@click.pass_context
def g00x(ctx: click.Context, logging_level: str | int = 0) -> None:
    logging.basicConfig(level=logging_level)
    ctx.obj = {}
    ctx.obj = {"data": Data(), "params": PlotParameters()}


g00x.add_command(plot_cli, name="plot")


@g00x.group("g002")
@click.pass_context
def g002(ctx: click.Context) -> None:
    """Run the G002 commands of G00x"""
    pass


@g00x.group("g003")
@click.pass_context
def g003(ctx: click.Context) -> None:
    """Run the G003 commands of G00x"""
    pass


@g002.group("box")
@click.pass_context
def box(ctx: click.Context) -> None:
    """Commands that will interact with NIHBox"""
    pass


@g002.group("globus")
@click.pass_context
def globus(ctx: click.Context) -> None:
    """Commands that will interact with Globus"""
    pass


@g002.group("validate")
@click.pass_context
def validate(ctx: click.Context) -> None:
    """Commands that will validate the file have all appropriate data and structure"""
    pass


@g002.group("pipeline")
@click.pass_context
@click.option(
    "--cellranger-path",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="The path to the celranger binary",
)
def pipeline(ctx: click.Context, cellranger_path: str | None) -> None:
    """Run the 10x pipeline including the SADIE AIRR output

    Parameters
    ----------
    cellranger_path : str | None
        Optionally describe where the cellranger binary is located. If not provided, it will be searched for in the path.
        If using Jordan's AMI, it is in /usr/local/bin/cellranger
    """
    data: Data = ctx.obj["data"]
    if cellranger_path:
        data.set_cellranger_path(cellranger_path)


@g003.group("pipeline")
@click.pass_context
@click.option(
    "--cellranger-path",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="The path to the celranger binary",
)
def g003_pipeline(ctx: click.Context, cellranger_path: str | None) -> None:
    """Run the 10x pipeline including the SADIE AIRR output

    Parameters
    ----------
    cellranger_path : str | None
        Optionally describe where the cellranger binary is located. If not provided, it will be searched for in the path.
        If using Jordan's AMI, it is in /usr/local/bin/cellranger
    """
    data: Data = ctx.obj["data"]
    if cellranger_path:
        data.set_cellranger_path(cellranger_path)


@g002.group("analysis")
def analysis():
    """Commands that will analyze the complete pipeline"""
    pass


@g003.group("analysis")
def g003_analysis():
    """Commands that will analyze the complete pipeline"""
    pass


@g003.group("validate")
@click.pass_context
def g003_validate(ctx: click.Context) -> None:
    """Commands that will validate the file have all appropriate data and structure"""
    pass


#####################
# Globus Commands
#####################
@globus.command("setup")
@click.pass_context
@click.option(
    "--mnt",
    "-m",
    type=click.Path(file_okay=False, dir_okay=True),
    default="/mnt/box",
)
@click.option(
    "--user",
    "-u",
    type=click.STRING,
    default=Path.home().name,
)
@click.option(
    "--group",
    "-g",
    type=click.STRING,
    default="schief",
)
def setup_globus(ctx: click.Context, mnt: Path | str, user: str, group: str) -> None:
    """Setup the NIHBox for the pipeline"""
    data = ctx.obj["data"]

    conda_env_template = Path(__file__).parent / "systemd_templates/conda_env"
    globus_service_session_template = Path(__file__).parent / "systemd_templates/globusconnectpersonal.service"
    globus_service_transfer_template = Path(__file__).parent / "systemd_templates/globustransfer.service"
    globus_timer_transfer_template = Path(__file__).parent / "systemd_templates/globustransfer.timer"

    try:
        ls = subprocess.run("which g00x", shell=True, capture_output=True)
        ls.check_returncode()
        g00x = ls.stdout.decode().strip()
    except subprocess.CalledProcessError as e:
        print("Error!\nreturn code: ", e.returncode, "\nOutput: ", e.stderr)
        raise

    cmd_string = f"""
        cd ~
        wget https://downloads.globus.org/globus-connect-personal/linux/stable/globusconnectpersonal-latest.tgz
        tar xzf globusconnectpersonal-latest.tgz
        rm -r ~/.globusconnectpersonal || exit 0
        mv globusconnectpersonal-*.*.*/ ~/.globusconnectpersonal

        mkdir -p ~/.config/systemd/user/
        cp {globus_timer_transfer_template} ~/.config/systemd/user/globustransfer.timer
        cp {globus_service_transfer_template} ~/.config/systemd/user/globustransfer.service
        systemctl --user import-environment
        systemctl --user daemon-reload

        ~/.globusconnectpersonal/globusconnectpersonal -setup
        globus login
        globus whoami
        globus endpoint local-id

        systemctl --user enable globusconnectpersonal.service
        systemctl --user enable globustransfer.service
        systemctl --user enable globustransfer.timer

        systemctl --user stop globusconnectpersonal.service

        systemctl --user start globusconnectpersonal.service
        systemctl --user start globustransfer.timer
        systemctl --user start globustransfer.service

        systemctl --user status globusconnectpersonal.service --no-pager
        systemctl --user status globustransfer.timer --no-pager
        systemctl --user status globustransfer.service --no-pager
    """

    cmds = [cmd.strip() for cmd in cmd_string.split("\n") if cmd.strip()]
    for cmd in cmds:
        print("$", cmd)
        try:
            ls = subprocess.run(cmd, shell=True)
            ls.check_returncode()
        except subprocess.CalledProcessError as e:
            print("Error!\nreturn code: ", e.returncode, "\nOutput: ", e.stderr)
            raise


@globus.command("start")
@click.pass_context
@click.option(
    "--src-uuid",
    "-s",
    type=click.STRING,
    default="81a01b62-2a01-11ed-8dd0-9f359c660fbd",
)
@click.option(
    "--dest-uuid",
    "-d",
    type=click.STRING,
)
@click.option(
    "--env",
    "-e",
    type=click.Path(file_okay=False, dir_okay=True),
)
def start_globus(ctx: click.Context, src_uuid: str, dest_uuid: str, env: str) -> None:
    """Start the globus service"""
    data = ctx.obj["data"]

    try:
        ls = subprocess.run(f"globus endpoint local-id", shell=True, capture_output=True)
        ls.check_returncode()
        dest_uuid = ls.stdout.decode().strip()
    except subprocess.CalledProcessError as e:
        print("Error!\nreturn code: ", e.returncode, "\nOutput: ", e.stderr)
        raise

    cmd = ["mkdir -p ~/g002/G002/sequencing/G002/"]
    run_numbers = [str(n).rjust(4, "0") for n in [15, 16, 17, 18, 19, 20]]
    for run_number in run_numbers:
        cmd.append(
            f'globus transfer {src_uuid}:/run{run_number}/ {dest_uuid}:~/g002/G002/sequencing/G002/run{run_number} --recursive --label "run{run_number}" --sync-level checksum -v --verify-checksum'
        )
    cmd_string = "\n".join(cmd)
    # cmd_string = f"""
    #     mkdir -p ~/g002/G002/sequencing/G002/
    #     globus transfer {src_uuid}:/run0002/ {dest_uuid}:~/g002/G002/sequencing/G002/run0002 --recursive --label "run0002" --sync-level checksum -v --verify-checksum
    #     globus transfer {src_uuid}:/run0003/ {dest_uuid}:~/g002/G002/sequencing/G002/run0003 --recursive --label "run0003" --sync-level checksum -v --verify-checksum
    #     globus transfer {src_uuid}:/run0004/ {dest_uuid}:~/g002/G002/sequencing/G002/run0004 --recursive --label "run0004" --sync-level checksum -v --verify-checksum
    #     globus transfer {src_uuid}:/run0005/ {dest_uuid}:~/g002/G002/sequencing/G002/run0005 --recursive --label "run0005" --sync-level checksum -v --verify-checksum
    #     globus transfer {src_uuid}:/run0006/ {dest_uuid}:~/g002/G002/sequencing/G002/run0006 --recursive --label "run0006" --sync-level checksum -v --verify-checksum
    #     globus transfer {src_uuid}:/run0007/ {dest_uuid}:~/g002/G002/sequencing/G002/run0007 --recursive --label "run0007" --sync-level checksum -v --verify-checksum
    # """
    cmds = [cmd.strip() for cmd in cmd_string.split("\n") if cmd.strip()]
    for cmd in cmds:
        print("$", cmd)
        try:
            ls = subprocess.run(cmd, shell=True)
            ls.check_returncode()
        except subprocess.CalledProcessError as e:
            print("Error!\nreturn code: ", e.returncode, "\nOutput: ", e.stderr)
            raise


# TODO: g00x globus stop


@globus.command("status")
@click.pass_context
def globus_status(ctx: click.Context) -> None:
    """Check the status of the Globus Transfer Service"""
    data = ctx.obj["data"]

    cmd_string = """
        systemctl --user status globusconnectpersonal.service --no-pager
        systemctl --user status globustransfer.timer --no-pager
        systemctl --user status globustransfer.service --no-pager -l
    """

    cmds = [cmd.strip() for cmd in cmd_string.split("\n") if cmd.strip()]
    for cmd in cmds:
        print("$", cmd, "\n")
        try:
            ls = subprocess.run(cmd, shell=True)
            print()
            if ls.returncode != 3:
                ls.check_returncode()
        except subprocess.CalledProcessError as e:
            print("Error!\nreturn code: ", e.returncode, "\nOutput: ", e.stderr)
            raise


#####################
# Box Commands
#####################
@box.command("setup")
@click.pass_context
@click.option(
    "--mnt",
    "-m",
    type=click.Path(file_okay=False, dir_okay=True),
    default="/mnt/box",
)
@click.option(
    "--user",
    "-u",
    type=click.STRING,
    default=Path.home().name,
)
@click.option(
    "--group",
    "-g",
    type=click.STRING,
    default="schief",
)
def setup_box(ctx: click.Context, mnt: Path | str, user: str, group: str) -> None:
    """Setup the NIHBox for the pipeline"""
    data = ctx.obj["data"]

    box_service_template = Path(__file__).parent / "systemd_templates/box.service"
    box_timer_template = Path(__file__).parent / "systemd_templates/box.timer"

    cmd_string = f"""
        # rclone config create box box
        sudo mkdir -p {mnt}
        sudo chown -R {user}:{group} {mnt} || exit 0
        rclone mount --allow-non-empty --daemon box: {mnt}
        mkdir -p ~/.config/systemd/user/
        cp {box_service_template} ~/.config/systemd/user/box.service
        cp {box_timer_template} ~/.config/systemd/user/box.timer
        # systemctl --user import-environment
        # systemctl --user daemon-reload
        # systemctl --user enable box.timer
        # systemctl --user start box.timer
        # systemctl --user start box.service
    """

    cmds = [cmd.strip() for cmd in cmd_string.split("\n") if cmd.strip()]
    for cmd in cmds:
        print("$", cmd)
        try:
            ls = subprocess.run(cmd, shell=True)
            ls.check_returncode()
        except subprocess.CalledProcessError as e:
            print("Error!\nreturn code: ", e.returncode, "\nOutput: ", e.stderr)
            raise


@box.command("status")
@click.pass_context
@click.option(
    "--mnt",
    "-m",
    type=click.STRING,
    default="box",
)
def box_status(ctx: click.Context, mnt: str) -> None:
    """Check the status of the NIHBox"""
    data = ctx.obj["data"]

    cmd_string = """
        rclone about box:
        systemctl --user status box.timer --no-pager
        systemctl --user status box.service --no-pager
    """

    cmds = [cmd.strip() for cmd in cmd_string.split("\n") if cmd.strip()]
    for cmd in cmds:
        print("$", cmd, "\n")
        try:
            ls = subprocess.run(cmd, shell=True)
            print()
            if ls.returncode != 3:
                ls.check_returncode()
        except subprocess.CalledProcessError as e:
            print("Error!\nreturn code: ", e.returncode, "\nOutput: ", e.stderr)
            raise


#####################
# Validate Commands
#####################
@validate.command("flow")
@click.pass_context
@click.argument("folder", type=click.Path(exists=True), required=True, default=".")
@click.option("--print_scheme", "-p", is_flag=True, default=False, help="Print scheme to stdout")
@click.help_option("--help", "-h", is_flag=True, help="Show this message and exit.")
def validate_flow(ctx: click.Context, folder: Path, print_scheme: bool) -> None:
    """
    Validate the flow data from NIHBox

    example:
        g00x validate flow /path/to/box/G002

    """
    if print_scheme:
        validate_g00x = ValidateG00X()
        print(validate_g00x)
        return

    data = ctx.obj["data"]
    # in validation, we will just run the parse flow data but just dump to the ether of the space-time contiuum
    parse_flow_data(data, folder)


@g003_validate.command("flow")
@click.pass_context
@click.argument("folder", type=click.Path(exists=True), required=True, default=".")
@click.option("--print_scheme", "-p", is_flag=True, default=False, help="Print scheme to stdout")
@click.help_option("--help", "-h", is_flag=True, help="Show this message and exit.")
def g003_validate_flow(ctx: click.Context, folder: Path, print_scheme: bool) -> None:
    """
    Validate the G003 flow data

    example:
        g00x validate flow /path/to/flow/G003

    """
    data = ctx.obj["data"]
    # in validation, we will just run the parse flow data but just dump to the ether of the space-time contiuum
    validate_g003_sorting(folder)


@validate.command("merge")
@click.pass_context
@click.option(
    "-f",
    "--flow-path",
    type=click.Path(exists=True),
    required=True,
    default=".",
    help="The path to the flow data. If this is passed and no flow_file is passed, then the workflow will be repeated",
)
@click.option(
    "-s",
    "--sequencing-path",
    type=click.Path(exists=True),
    required=True,
    default=".",
)
@click.option(
    "--out",
    "-o",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="merged_output",
    help="The output the merged flow and sequencing data",
)
def merge(ctx: click.Context, flow_path: Path, sequencing_path: Path, out: Path) -> None:
    """Merge the sequencing and flow data into a single dataframe

    example:
        g00x sequencing --flow-file output/flow_output.feather -s /path/to/sequencing

    """
    data = ctx.obj["data"]

    # Merge but throw to space time
    df = merge_flow_and_sequencing(data, flow_path, sequencing_path)
    click.echo("Merged flow and sequencing data")
    click.echo(f"Writting to {out}.feather/.csv.gz")
    df.to_csv(str(out) + ".csv")
    df.to_feather(str(out) + ".feather")


@g003_validate.command("sequencing")
@click.pass_context
@click.option(
    "-s",
    "--sequencing-path",
    type=click.Path(exists=True),
    required=True,
    default=".",
    help="The path to G003 Sequencing directory",
)
@click.option(
    "--manifest-name",
    "-m",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="merged_sequencing_manifest",
    help="The output the merged sequencing manifests",
)
@click.option(
    "--out",
    "-o",
    type=click.Path(file_okay=False, dir_okay=True, writable=True),
    default=".",
    help="The output directory for the merged sequencing manifest",
)
def g003_merge(
    ctx: click.Context,
    sequencing_path: Path,
    manifest_name,
    out: Path,
) -> None:
    """Validate Sequencing paths and Merge the sequencing manifests into a single dataframe

    example:
        g00x g003 validate sequencing --s /path/to/sequencing

    """
    # data = ctx.obj["data"]
    manifest_name = Path(manifest_name)
    df = validate_g003_sequencing(sequencing_path)
    click.echo("Merged sequencing manifests")
    click.echo(f"Writting to {manifest_name.stem}.feather/.csv.gz")
    df.to_csv(f"{out}/{manifest_name.stem}.csv", index=False)
    df.to_feather(f"{out}/{manifest_name.stem}.feather")


#####################
# Pipeline Commands #
#####################
@pipeline.command("flow")
@click.pass_context
@click.option(
    "--out",
    "-o",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="flow_output",
    help="The output the flow data.",
)
@click.argument(
    "folder",
    type=click.Path(exists=True),
    required=True,
    default=".",
)
def parse_flow(ctx: click.Context, out: Path, folder: Path) -> None:
    """Parse the flow into a flow dataframe

    Parameters
    ----------
    out : Path
        The output path to the parsed flow datafraem
    folder : Path
        The flow folder path, .e.g /path/to/box/G002

    """
    # get the flow dataframe back
    data = ctx.obj["data"]
    flow_data = parse_flow_data(data, folder)
    out = Path(out)
    output_feather = Path(out.parent / (out.stem + ".feather"))
    output_csv = Path(out.parent / (out.stem + ".csv"))
    click.echo(f"Writing to {output_feather}")
    flow_data.to_feather(output_feather)
    click.echo(f"Writing to {output_csv}")
    flow_data.to_csv(output_csv)


@pipeline.command("demultiplex")
@click.pass_context
@click.option(
    "--out",
    "-o",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="demultiplex_output",
    help="The dataframe with demultiplexed paths",
)
@click.option(
    "-f",
    "--flow-path",
    type=click.Path(exists=True, dir_okay=True, resolve_path=True),
    required=True,
    default=None,
    help="The path to the flow data. If this is passed and no flow_file is passed, then the workflow will be repeated",
)
@click.option(
    "-s",
    "--sequencing-path",
    type=click.Path(exists=True, dir_okay=True, resolve_path=True),
    required=True,
    default=None,
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite the demultiplex and run again",
)
def demultiplex(
    ctx: click.Context,
    out: Path,
    flow_path: Path | None,
    sequencing_path: Path | None,
    overwrite: bool,
) -> None:
    """
    Run demultiplexing on the merged dataframe.
    Either the dataframe must be passed or the flow and sequencing paths must be passed.

    Parameters
    ----------
    out : Path
        The output path for the demultiplexed dataframe
    flow_path : Path
        The path to the flow data. Needs to be passed if no merged dataframe
    sequencing_path : Path
        The path to the sequencing data. Needs to be passed if no merged dataframe
    """
    data = ctx.obj["data"]
    click.echo(f"Merging data with flow path {flow_path} and sequencing path {sequencing_path}")
    merged_dataframe: pd.DataFrame = merge_flow_and_sequencing(data, flow_path, sequencing_path)  # type: ignore
    run_demultiplex(data, merged_dataframe, out, overwrite)


@g003_pipeline.command("flow")
@click.pass_context
@click.option(
    "--out",
    "-o",
    type=click.Path(file_okay=False, dir_okay=True, writable=True),
    default=".",
    help="The output folder the flow data csv/feather",
)
@click.option(
    "--flow-name",
    "-f",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="flow",
    help="The output prefix the flow data csv/feather",
)
@click.argument(
    "folder",
    type=click.Path(exists=True),
    required=True,
    default=".",
)
def g003_parse_flow(ctx: click.Context, out: Path, flow_name: Path, folder: Path) -> None:
    """Parse the flow sorts into a flow dataframe

    Parameters
    ----------
    out : Path
        The output path to the parsed flow csv & dataframe
    flow_output : str
        Stem name of flow files
    folder : Path
        The flow folder path, .e.g ~/g003_bucket/g003/g003/sorting/G003

    """
    # get the flow dataframe back
    data = ctx.obj["data"]
    # legacy code
    ptid2pubid = {}  # data.get_g003_pubids_lookup()
    ptid_prefix2group = data.get_g003_ptid_prefix_2_group()
    visit_id2week = data.get_g003_visit_id_2_week()

    # Input Paths
    out = pathing(out)
    folder = pathing(folder)
    # Output Paths
    output_feather = (out / flow_name).with_suffix(".feather")
    output_csv = (out / flow_name).with_suffix(".csv")
    output_unmerged_feather = (out / (flow_name + "-unmerged")).with_suffix(".feather")
    output_unmerged_csv = (out / (flow_name + "-unmerged")).with_suffix(".csv")

    validation = validate_g003_sorting(folder)

    flow_df = g003_flow.pull_flow_from_validation(
        validation=validation, ptid2pubid=ptid2pubid, ptid_prefix2group=ptid_prefix2group, visit_id2week=visit_id2week
    )
    click.echo(f"Writing to {output_feather}")
    flow_df.to_feather(output_feather)
    click.echo(f"Writing to {output_csv}")
    flow_df.to_csv(output_csv)

    # click.echo(f"Writing to {output_unmerged_feather}")
    # flow_unmerged_df.to_feather(output_unmerged_feather)
    # click.echo(f"Writing to {output_unmerged_csv}")
    # flow_unmerged_df.to_csv(output_unmerged_csv)


## demultiplex options for g003
@g003_pipeline.command("demultiplex")
@click.pass_context
@click.option(
    "--demultiplex-output",
    "-d",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="demultiplex_output_dataframe",
    help="The dataframe with demultiplexed paths",
)
@click.option(
    "-s",
    "--sequencing-path",
    type=click.Path(exists=True, dir_okay=True, resolve_path=True),
    required=True,
    default=None,
    help="The path to the G003 sequencing directory (Must be named G003)",
)
@click.option(
    "-o",
    "--out",
    type=click.Path(exists=False, dir_okay=True, resolve_path=False),
    required=False,
    default="working_directory",
    help="The path to an output directory (Default is a working_directory within the sequencing run",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite the demultiplex and run again",
)
@click.option(
    "--run",
    default=None,
    show_default=True,
    help="Demultiplex only this run. If not provided, all runs will be demultiplexed",
    multiple=False,
)
def g003_demultiplex(
    ctx: click.Context,
    demultiplex_output: str,
    out: Path,
    sequencing_path: Path | None,
    overwrite: bool,
    run: str | None,
) -> None:
    """
    Run demultiplexing on the merged dataframe.
    Either the dataframe must be passed or the flow and sequencing paths must be passed.

    Parameters
    ----------
    demulitplex_output: Path
        The file name for saving the demultiplex dataframe
    sequencing_path : Path
        The path to the sequencing data. Needs to be passed if no merged dataframe
    out : Path
        The output path for the working directory to use
    overwrite : bool
        Overwrite the demultiplex and run again
    run : list[str]
        Demultiplex only this run. If not provided, all runs will be demultiplexed


    >>> g00x g003 pipeline demultiplex -s ./g003_bucket/g003/g003/sequencing/G003/ --run run0001 --run --run0003
    """
    data = ctx.obj["data"]
    out = pathing(out)
    if not out.exists():
        out.mkdir(out)  # type: ignore

    merged_dataframe: pd.DataFrame = validate_g003_sequencing(sequencing_path)  # type: ignore
    if run:
        merged_dataframe = merged_dataframe[merged_dataframe["run_dir_path"].str.endswith(run)]
        demultiplex_output = run + "/" + demultiplex_output

    demultiplexed_dataframe = g003_run_demultiplex(data, merged_dataframe, out, overwrite)

    demultiplexed_dataframe = demultiplexed_dataframe.applymap(pd_replace_home_with_tilde)

    demultiplexed_dataframe.to_csv(out / f"{demultiplex_output}.csv")
    demultiplexed_dataframe.to_feather(out / f"{demultiplex_output}.feather")


@pipeline.command("vdj")
@click.pass_context
@click.option(
    "--out",
    "-o",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="vdj_demultiplex_output",
    help="The dataframe with demultiplexed paths",
)
@click.option(
    "-d",
    "--demultiplex-dataframe-path",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
    required=True,
    help="The path to the demultiplexed dataframe from the demultiplexed pipeline",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite the vdj files and run again",
)
def vdj(
    ctx: click.Context,
    out: Path,
    demultiplex_dataframe_path: Path,
    overwrite: bool,
) -> None:
    """
    Run vdj on the demultiplex dataframe.
    Either the dataframe must be passed or the flow and sequencing paths must be passed.
    If the flow and sequencing are passed, we will run the demultiplexing step first.

    Parameters
    ----------
    out : Path
        The output path for the demultiplexed dataframe
    demultiplex_dataframe_path : Path
        The demultiplexed dataframe from the demultiplexed pipeline. If not provided, the flow and sequencing paths must be provided
    """
    data = ctx.obj["data"]
    demultiplex_dataframe = pd.read_feather(Path(demultiplex_dataframe_path))
    click.echo("Running VDJ pipeline")
    run_vdj(data, demultiplex_dataframe, out, overwrite)


@g003_pipeline.command("vdj")
@click.pass_context
@click.option(
    "-o",
    "--out",
    type=click.Path(exists=False, dir_okay=True, resolve_path=False),
    required=False,
    default="working_directory",
    help="The path to an output directory (Default is a working_directory within the sequencing run",
)
@click.option(
    "-d",
    "--demultiplex-dataframe-path",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
    required=True,
    help="The path to the demultiplexed dataframe from the demultiplexed pipeline",
)
@click.option(
    "--vdj-frame-output",
    "-v",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="vdj_demultiplex_output",
    help="The merged dataframe with demultiplexed paths",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite the vdj files and run again",
)
def g003_vdj(
    ctx: click.Context,
    out: Path,
    vdj_frame_output: Path,
    demultiplex_dataframe_path: Path,
    overwrite: bool,
) -> None:
    """
    Run vdj on the demultiplex dataframe.
    Either the dataframe must be passed or the flow and sequencing paths must be passed.
    If the flow and sequencing are passed, we will run the demultiplexing step first.

    Parameters
    ----------
    out : Path
        The output path for the demultiplexed dataframe
    demultiplex_dataframe_path : Path
        The demultiplexed dataframe from the demultiplexed pipeline. If not provided, the flow and sequencing paths must be provided
    """
    data = ctx.obj["data"]
    out = pathing(out)
    if not out.exists():
        out.mkdir(out)

    demultiplex_dataframe = pd.read_feather(Path(demultiplex_dataframe_path))
    demultiplex_dataframe = demultiplex_dataframe.applymap(pd_expand_path)

    click.echo(f"Running VDJ pipeline in {out}")
    demultiplexed_dataframe = g003_run_vdj(data, demultiplex_dataframe, out, overwrite)
    demultiplexed_dataframe = demultiplexed_dataframe.applymap(pd_replace_home_with_tilde)

    demultiplexed_dataframe.to_csv(out / f"{vdj_frame_output}.csv")
    demultiplexed_dataframe.to_feather(out / f"{vdj_frame_output}.feather")


@pipeline.command("cso")
@click.pass_context
@click.option(
    "--out",
    "-o",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="cso_demultiplex_output",
    help="The dataframe with demultiplexed paths",
)
@click.option(
    "-d",
    "--demultiplex-dataframe-path",
    type=click.Path(exists=True),
    required=True,
    help="The path to the demultiplexed dataframe from the demultiplexed pipeline or the vdj pipeline",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite the cso files and run again",
)
def cso(
    ctx: click.Context,
    out: Path,
    demultiplex_dataframe_path: Path,
    overwrite: bool,
) -> None:
    """
    Run cso on the demultiplex dataframe.
    Either the dataframe must be passed or the flow and sequencing paths must be passed.
    If the flow and sequencing are passed, we will run the demultiplexing step first.

    Parameters
    ----------
    out : Path
        The output path for the demultiplexed dataframe
    demultiplex_dataframe_path : Path
        The demultiplexed dataframe from the demultiplexed pipeline. If not provided, the flow and sequencing paths must be provided
    overwrite : bool
        Overwrite the cso files and run again
    """
    data = ctx.obj["data"]
    click.echo("Reading Demultiplexed Dataframe")
    demultiplex_dataframe = pd.read_feather(demultiplex_dataframe_path)
    run_cso(data, demultiplex_dataframe, out, overwrite)


@g003_pipeline.command("cso")
@click.pass_context
@click.option(
    "-o",
    "--out",
    type=click.Path(exists=False, dir_okay=True, resolve_path=False),
    required=False,
    default="working_directory",
    help="The path to an output directory (Default is a working_directory within the sequencing run",
)
@click.option(
    "--cso-frame-output",
    "-c",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="cso_demultiplex_output",
    help="The dataframe with demultiplexed paths",
)
@click.option(
    "-d",
    "--demultiplex-dataframe-path",
    type=click.Path(exists=True),
    required=True,
    help="The path to the demultiplexed dataframe from the demultiplexed pipeline or the vdj pipeline",
)
@click.option(
    "-g",
    "--genome-reference",
    type=click.Path(exists=True),
    required=False,
    help="The path to the human genome",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite the cso files and run again",
)
def g003_cso(
    ctx: click.Context,
    demultiplex_dataframe_path: Path,
    cso_frame_output: Path,
    out: Path,
    genome_reference: Path,
    overwrite: bool,
) -> None:
    """
    Run cso on the demultiplex dataframe.
    Either the dataframe must be passed or the flow and sequencing paths must be passed.
    If the flow and sequencing are passed, we will run the demultiplexing step first.

    Parameters
    ----------
    out : Path
        The output path for the demultiplexed dataframe
    demultiplex_dataframe_path : Path
        The demultiplexed dataframe from the demultiplexed pipeline. If not provided, the flow and sequencing paths must be provided
    overwrite : bool
        Overwrite the cso files and run again
    """
    data = ctx.obj["data"]
    out = pathing(out)
    if not out.exists():
        out.mkdir(out)

    click.echo("Reading Demultiplexed Dataframe")
    demultiplex_dataframe = pd.read_feather(demultiplex_dataframe_path)
    demultiplex_dataframe = demultiplex_dataframe.applymap(pd_expand_path)

    demultiplexed_dataframe = g003_run_cso(data, demultiplex_dataframe, out, genome_reference, overwrite)
    demultiplexed_dataframe = demultiplexed_dataframe.applymap(pd_replace_home_with_tilde)

    demultiplexed_dataframe.to_csv(out / f"{cso_frame_output}.csv")
    demultiplexed_dataframe.to_feather(out / f"{cso_frame_output}.feather")


@pipeline.command("airr")
@click.pass_context
@click.option(
    "--vdj-out",
    "-v",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    required=True,
    help="The dataframe output from the vdj pipeline",
)
@click.option(
    "--cso-out",
    "-c",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    required=True,
    help="The dataframe output from cso pipeline",
)
@click.option(
    "--out",
    "-o",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="combined_airr",
    help="The dataframe output for the combined airr pipeline and meta",
)
@click.option(
    "--skip-mutation",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    "--cluster-heavy-only",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    "--cluster-n",
    "-k",
    type=click.INT,
    default=5,
    show_default=True,
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite the airr files and run again",
)
def airr(
    ctx: click.Context,
    vdj_out: Path,
    cso_out: Path,
    out: Path,
    skip_mutation: bool,
    cluster_heavy_only: bool,
    cluster_n: int,
    overwrite: bool,
) -> None:
    """
    Run the airr pipeline on the vdj and cso dataframes.

    Returning a combined dataframes with SADIE AIRR annotations and meta data

    Parameters
    ----------
    vdj_out : Path
        The vdj pipeline dataframe
    cso_out : Path
        The cso pipeline dataframe
    out : Path
        The output path for the combined airr dataframe
    """
    data = ctx.obj["data"]
    click.echo("Reading in vdj and cso dataframes")
    vdj_dataframe = pd.read_feather(vdj_out)
    cso_dataframe = pd.read_feather(cso_out)
    print(f"Running AIRR pipeline {cluster_n} {cluster_heavy_only}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)
        warnings.simplefilter("ignore", category=PerformanceWarning)
        _ = run_airr(
            data,
            vdj_dataframe,
            cso_dataframe,
            out,
            overwrite,
            skip_mutation,
            cluster_n=cluster_n,
            cluster_heavy_only=cluster_heavy_only,
        )


@g003_pipeline.command("airr")
@click.pass_context
@click.option(
    "--vdj-out",
    "-v",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    required=True,
    help="The dataframe output from the vdj pipeline",
)
@click.option(
    "--cso-out",
    "-c",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    required=True,
    help="The dataframe output from cso pipeline",
)
@click.option(
    "-o",
    "--out",
    type=click.Path(exists=True, dir_okay=True, resolve_path=False),
    required=False,
    default=".",
    help="The path to an output directory (Default is a working_directory within the sequencing run",
)
@click.option(
    "--airr-frame-output",
    "-a",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="combined_airr",
    help="The dataframe output for the combined airr pipeline and meta",
)
@click.option(
    "--skip-mutation",
    is_flag=True,
    default=False,
    show_default=True,
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite the airr files and run again",
)
def g003_airr(
    ctx: click.Context,
    vdj_out: Path,
    cso_out: Path,
    airr_frame_output: Path,
    out: Path,
    skip_mutation: bool,
    overwrite: bool,
) -> None:
    """
    Run the airr pipeline on the vdj and cso dataframes.

    Returning a combined dataframes with SADIE AIRR annotations and meta data

    Parameters
    ----------
    vdj_out : Path
        The vdj pipeline dataframe
    cso_out : Path
        The cso pipeline dataframe
    out : Path
        The output path for the combined airr dataframe
    """
    data = ctx.obj["data"]
    out = pathing(out)
    if not out.exists():
        out.mkdir(out)
    if (out / f"{airr_frame_output}.feather").exists() and not overwrite:
        click.echo(f"{out / f'{airr_frame_output}.feather'} exists. Skipping.")
        return
    click.echo("Reading in vdj and cso dataframes")
    vdj_dataframe = pd.read_feather(vdj_out).applymap(pd_expand_path)
    cso_dataframe = pd.read_feather(cso_out).applymap(pd_expand_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)
        warnings.simplefilter("ignore", category=PerformanceWarning)
        combined_airr = g003_run_airr(data, vdj_dataframe, cso_dataframe, out, overwrite, skip_mutation)
        combined_airr = combined_airr.applymap(pd_replace_home_with_tilde)
        combined_airr.to_csv(out / f"{airr_frame_output}.csv")
        combined_airr.to_feather(out / f"{airr_frame_output}.feather")


@g003_pipeline.command("merge")
@click.pass_context
@click.option(
    "--out",
    "-o",
    type=click.Path(file_okay=False, dir_okay=True, writable=True),
    default=".",
    help="The output folder the flow data csv/feather",
)
@click.option(
    "--name",
    "-n",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="flow_manifest",
    help="The output prefix the flow data csv/feather",
)
@click.option(
    "--flow-path",
    "-f",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    help="The output path to flow",
)
@click.option(
    "--seq-manifest-path",
    "-s",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    help="The output path to seq",
)
def g003_merge_flow_and_airr(
    ctx: click.Context, out: Path, name: str, flow_path: Path, seq_manifest_path: Path
) -> None:
    """Parse the flow sorts into a flow dataframe

    Parameters
    ----------
    out : Path
        The output path to the parsed flow csv & dataframe
    name: str
        Stem name of merged flow and airr files
    flow_path : Path
        Path of the flow output feather
    seq_manifest_path : Path
        Path of the sequence manifest feather
    """
    # get the flow dataframe back
    # data = ctx.obj["g003_data"]

    # Input Paths
    out = pathing(out)
    flow_path = pathing(flow_path)
    seq_manifest_path = pathing(seq_manifest_path)

    flow_df = pd.read_feather(flow_path) if flow_path.suffix == ".feather" else pd.read_csv(flow_path)
    seq_manifest_df = (
        pd.read_feather(seq_manifest_path) if seq_manifest_path.suffix == ".feather" else pd.read_csv(seq_manifest_path)
    )
    seq_manifest_df = seq_manifest_df
    seq_manifest_df["sorted_date"] = seq_manifest_df["sorted_date"].astype(str)
    flow_df["run_date"] = flow_df["run_date"].astype(str)
    unique_cols: list[str] = [
        "run_purpose",
        "ptid",
        "group",
        "weeks",
        "visit_id",
        "probe_set",
        "sample_type",
        "run_date",
        "sort_pool",
        "hashtag",
    ]
    flow_manifest = flow_df.drop_duplicates(subset=unique_cols, keep="first")[unique_cols].merge(
        seq_manifest_df,
        left_on=["ptid", "visit_id", "run_date", "sort_pool"],
        right_on=["ptid", "timepoint", "sorted_date", "pool_number"],
        how="outer",
    )
    # Output Paths
    output_feather = (out / name).with_suffix(".feather")
    output_csv = (out / name).with_suffix(".csv")

    click.echo(f"Writing to {output_feather}")
    flow_manifest.to_feather(output_feather)
    click.echo(f"Writing to {output_csv}")
    flow_manifest.to_csv(output_csv)


@pipeline.command("e2e")
@click.pass_context
@click.option(
    "-f",
    "--flow-path",
    type=click.Path(exists=True),
    required=True,
    default=None,
    help="The path to the flow data. If this is passed and no flow_file is passed, then the workflow will be repeated",
)
@click.option(
    "-s",
    "--sequencing-path",
    type=click.Path(exists=True),
    required=True,
    default=None,
)
@click.option(
    "--out",
    "-o",
    type=click.Path(file_okay=False, dir_okay=True, writable=True),
    default="g002/G002/output",
    help="The output the pipeline",
)
@click.option(
    "--overwrite",
    is_flag=True,
    default=False,
    show_default=True,
    help="Overwrite the demultiplex and run again",
)
def run_e2e(
    ctx: click.Context,
    flow_path: Path,
    sequencing_path: Path,
    out: Path,
    overwrite: bool,
) -> None:
    """Run the end to end pipeline: validation, flow, demultiplexing, vdj, cso and airr

    We expect the flow and sequences to be in their proper file scheme for this to work.

    Parameters
    ----------
    flow_path : Path
        The path to the flow files e.g, g002/G002/sorting/G002/
    sequencing_path : Path
        The path to the e.g., g002/G002/sequencing/G002/
    out: Path
        Path for the complete output for the g00x pipeline
    overwrite : bool
        even if output files exist, overwrite them anyway
    """
    data = ctx.obj["data"]
    flow_path = pathing(flow_path)
    sequencing_path = pathing(sequencing_path)

    # Pop into output directory
    with cd(out):
        # Flow
        ctx.invoke(parse_flow, folder=flow_path)
        # Merge
        ctx.invoke(merge, flow_path=flow_path, sequencing_path=sequencing_path)
        # Demultiplex
        ctx.invoke(demultiplex, flow_path=flow_path, sequencing_path=sequencing_path)
        # VDJ
        ctx.invoke(vdj, demultiplex_dataframe_path="demultiplex_output.feather")
        # CSO
        ctx.invoke(cso, demultiplex_dataframe_path="demultiplex_output.feather")
        # AIRR
        ctx.invoke(
            airr,
            vdj_out="vdj_demultiplex_output.feather",
            cso_out="cso_demultiplex_output.feather",
        )


#####################
# Analysis Commands #
#####################
@analysis.command("count")
@click.pass_context
@click.option(
    "-f",
    "--flow-path",
    type=click.Path(exists=True),
    required=True,
    default=None,
    help="The path to the flow data. From the flow pipeline",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="count",
    help="The output path for the count dataframe",
)
def count(ctx: click.Context, flow_path: Path, output: Path) -> None:
    """
    Count the samples that have been processed

    Parameters
    ----------
    flow_path : Path
        The path to the flow data. e.g. flow output from the pipeline flow command
    """
    data = ctx.obj["data"]
    click.echo(f"Counting samples in {flow_path}")
    flow_data = pd.read_feather(flow_path)
    count_current_samples(data, flow_data, output)
    click.echo(f"Counted samples written to {output}.png")


@analysis.command("report")
@click.pass_context
@click.option(
    "-o",
    "--output",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="flow_and_sequencing",
    help="The output path for the flow and sequencing dataframe",
)
@click.option(
    "-s",
    "--sequencing-dataframe-path",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
    default=None,
    help="The path to the sequencing dataframe. From the airr pipeline. Give feather",
)
@click.option(
    "-f",
    "--flow-dataframe-path",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
    default=None,
    help="The path to the flow dataframe. From the flow pipeline. Give feather",
)
def generate_report(
    ctx: click.Context,
    output: Path,
    sequencing_dataframe_path: Path,
    flow_dataframe_path: Path,
) -> None:
    """
    Generate a report of the flow and sequencing data. These will most likely be used to plot everything else
    """
    click.echo("Generating flow and sequencing report")
    data = ctx.obj["data"]
    sequencing_dataframe = pd.read_feather(sequencing_dataframe_path)
    click.echo(f"Filtering sequencing dataframe to only PBMC samples before {len(sequencing_dataframe)}")
    sequening_dataframe = sequencing_dataframe.query("sample_type=='PBMC'")
    click.echo(f"After {len(sequencing_dataframe)}")
    flow_dataframe = pd.read_feather(flow_dataframe_path)
    click.echo(f"Filtering sequencing dataframe to only PBMC samples before {len(flow_dataframe)}")
    flow_dataframe = flow_dataframe.query("sample_type=='PBMC'")
    click.echo(f"After {len(flow_dataframe)}")

    (
        seq_and_flow_df,
        seq_and_flow_df_long_name,
        seq_and_flow_df_long_form,
    ) = combine_seq_and_flow(data, sequencing_dataframe, flow_dataframe)

    # compact name pivot
    seq_and_flow_df.to_feather(str(output) + ".feather")
    seq_and_flow_df.to_csv(str(output) + ".csv")
    click.echo(f"Flow and sequencing report written to {output}.feather/.csv")

    # long name pivot
    seq_and_flow_df_long_name.to_feather(str(output) + "_long_names.feather")
    seq_and_flow_df_long_name.to_csv(str(output) + "_long_names.csv")
    click.echo(f"Flow and sequencing report written to {output}_long_names.feather/.csv")

    # Long form for the sane
    seq_and_flow_df_long_form.to_feather(str(output) + "_long_form.feather")
    seq_and_flow_df_long_form.to_csv(str(output) + "_long_form.csv")
    click.echo(f"Flow and sequencing report written to {output}_long_form.feather/.csv")


@g003_analysis.command("report")
@click.pass_context
@click.option(
    "--out",
    "-o",
    type=click.Path(file_okay=False, dir_okay=True, writable=True),
    default=".",
    help="The output the pipeline report",
)
@click.option(
    "-r",
    "--report-output",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
    default="flow_and_sequencing",
    help="The output path for the flow and sequencing dataframe",
)
@click.option(
    "-s",
    "--sequencing-dataframe-path",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
    default=None,
    help="The path to the sequencing dataframe. From the airr pipeline. Give feather",
)
@click.option(
    "-f",
    "--flow-dataframe-path",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
    default=None,
    help="The path to the flow dataframe. From the flow pipeline. Give feather",
)
def g003_generate_report(
    ctx: click.Context,
    out: Path,
    report_output: Path,
    sequencing_dataframe_path: Path,
    flow_dataframe_path: Path,
) -> None:
    """
    Generate a report of the flow and sequencing data. These will most likely be used to plot everything else
    """
    click.echo("Generating flow and sequencing reports")
    data = ctx.obj["data"]
    out = pathing(out)
    if not out.exists():
        out.mkdir(out)

    sequencing_dataframe_path = pathing(sequencing_dataframe_path)
    flow_dataframe_path = pathing(flow_dataframe_path)

    sequencing_dataframe = pd.read_feather(sequencing_dataframe_path)
    sequencing_dataframe = sequencing_dataframe.applymap(pd_expand_path)
    sequencing_dataframe["sorted_date"] = sequencing_dataframe["sorted_date"].astype(str)

    # sequencing_dataframe["run_date"] = sequencing_dataframe["run_date"].astype(str)

    click.echo(f"Filtering sequencing dataframe to only PBMC samples before {len(sequencing_dataframe)}")
    click.echo(f"After {len(sequencing_dataframe)}")

    flow_dataframe = pd.read_feather(flow_dataframe_path)
    # flow_dataframe["run_date"] = flow_dataframe["run_date"].astype(str)
    click.echo(f"Filtering sequencing dataframe to only PBMC samples before {len(flow_dataframe)}")

    flow_dataframe = flow_dataframe.query("sample_type=='PBMC'")
    click.echo(f"After {len(flow_dataframe)}")

    (
        seq_and_flow_df,
        seq_and_flow_df_long_name,
        seq_and_flow_df_long_form,
        # seq_and_flow_df_long_calc,
    ) = g003_combine_seq_and_flow(data, sequencing_dataframe, flow_dataframe)
    seq_and_flow_df = seq_and_flow_df.applymap(pd_replace_home_with_tilde)
    seq_and_flow_df_long_name = seq_and_flow_df_long_name.applymap(pd_replace_home_with_tilde)
    seq_and_flow_df_long_form = seq_and_flow_df_long_form.applymap(pd_replace_home_with_tilde)
    # seq_and_flow_df_long_calc = seq_and_flow_df_long_calc.applymap(pd_replace_home_with_tilde)

    # compact name pivot
    seq_and_flow_df.to_feather(out / f"{report_output}.feather")
    seq_and_flow_df.to_csv(out / f"{report_output}.csv")
    click.echo(f"Flow and sequencing report written to {report_output}.feather/.csv")

    # long name pivot
    seq_and_flow_df_long_name.to_feather(out / f"{report_output}_long_names.feather")
    seq_and_flow_df_long_name.to_csv(out / f"{report_output}_long_names.csv")
    click.echo(f"Flow and sequencing report written to {report_output}_long_names.feather/.csv")

    # Long form for the sane
    seq_and_flow_df_long_form.to_feather(out / f"{report_output}_long_form.feather")
    seq_and_flow_df_long_form.to_csv(out / f"{report_output}_long_form.csv")
    click.echo(f"Flow and sequencing report written to {report_output}_long_form.feather/.csv")

    # Long form for the sane
    # seq_and_flow_df_long_calc.to_feather(out / f"{report_output}_long_calc.feather")
    # seq_and_flow_df_long_calc.to_csv(out / f"{report_output}_long_calc.csv")
    # click.echo(f"Flow and sequencing report written to {report_output}_long_calc.feather/.csv")


if __name__ == "__main__":
    try:
        g00x()  # type: ignore
    except Exception as e:
        import traceback

        traceback.print_exc()
