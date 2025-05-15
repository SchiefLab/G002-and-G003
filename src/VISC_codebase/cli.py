import os
import subprocess
from pathlib import Path

import click


class cd:
    """Context manager for changing the current working directory"""

    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def run_r_script(rscript_path: str, sourcedir: str = "../../G00x-plots/") -> None:
    """Run an R script with the given source directory."""
    script_dir = Path(__file__).parent
    rscript_path = script_dir / Path(rscript_path)
    if not rscript_path.exists():
        click.echo(f"R script not found: {rscript_path}")
        raise FileNotFoundError(f"R script not found: {rscript_path}")
    if not rscript_path.is_file():
        click.echo(f"R script path is not a file: {rscript_path}")
        raise ValueError(f"R script path is not a file: {rscript_path}")
    if rscript_path.suffix.lower() != ".r":
        click.echo(f"R script does not have .R extension: {rscript_path}")
        raise ValueError(f"R script does not have .R extension: {rscript_path}")
    click.echo(f"Rscript '{rscript_path}' '{sourcedir}'")
    # Needs to be run from the directory of the script
    with cd(script_dir):
        click.echo(Path.cwd())
        rscript_cmd = f"Rscript '{rscript_path}' '{sourcedir}'"
        subprocess.run(rscript_cmd, shell=True, check=True)


def render_rmarkdown(rmarkdown_path: str, sourcedir: str = "../../G00x-plots/") -> None:
    """Render an R markdown file."""
    rmd = Path(rmarkdown_path)
    if not rmd.exists():
        raise FileNotFoundError(f"R markdown not found: {rmarkdown_path}")
    if not rmd.is_file():
        raise ValueError(f"R markdown path is not a file: {rmarkdown_path}")
    if rmd.suffix.lower() not in [".rmd"]:
        raise ValueError(f"R markdown does not have .Rmd extension: {rmarkdown_path}")
    rmarkdown_cmd = f"Rscript -e \"rmarkdown::render('{rmarkdown_path}',  params=list(data_location='{sourcedir}'))\""
    subprocess.run(rmarkdown_cmd, shell=True, check=True)


@click.command()
def cli():
    pass


@click.command("comparison-tables")
@click.option(
    "--sourcedir",
    default=str(Path(__file__).parent.parent.parent / "G00x-plots/"),
    show_default=True,
    help="Path to the source directory for the R script.",
)
def comparison_tables(sourcedir):
    """Run the comparison rcripts."""
    # render_rmarkdown(rmarkdown_path='FigS8_TabS31-34_TabS77.Rmd', sourcedir=sourcedir)
    # render_rmarkdown(rmarkdown_path='Figure_S10.Rmd', sourcedir=sourcedir)
    # render_rmarkdown(rmarkdown_path='TableS35_S36.Rmd', sourcedir=sourcedir)
    click.echo("Running mk_table_figure2.R")
    run_r_script(rscript_path="mk_table_figure2.R", sourcedir=sourcedir)
    run_r_script(rscript_path="FigureS9.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S40.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S54_S56.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S55.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S57.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S58.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S59.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S60.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S61.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S62.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S63.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S64.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S66.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S67.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S68_fig4.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S68_fig6.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S71.R", sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S78.R", sourcedir=sourcedir)


@click.command("make-tables")
@click.option(
    "--sourcedir", default="../../G00x-plots/", show_default=True, help="Path to the source directory for the R script."
)
def make_tables(sourcedir):
    """Run the Rmd scripts."""
    render_rmarkdown(rmarkdown_path="TableS40.Rmd", sourcedir=sourcedir)
    render_rmarkdown(rmarkdown_path="FigS8_TabS31-34_TabS77.Rmd", sourcedir=sourcedir)
    render_rmarkdown(rmarkdown_path="Figure_S10.Rmd", sourcedir=sourcedir)
    render_rmarkdown(rmarkdown_path="TableS35_S36.Rmd", sourcedir=sourcedir)


# Create a click group to include both commands
@click.group()
def main():
    """Run R scripts."""
    pass


main.add_command(cli)
main.add_command(comparison_tables)

if __name__ == "__main__":
    main()  # noqa: F403
