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
    script_dir = Path(__file__).parent
    rmarkdown_path = script_dir / Path(rmarkdown_path)
    if not rmarkdown_path.exists():
        click.echo(f"R markdown not found: {rmarkdown_path}")
        raise FileNotFoundError(f"R markdown not found: {rmarkdown_path}")
    if not rmarkdown_path.is_file():
        click.echo(f"R markdown path is not a file: {rmarkdown_path}")
        raise ValueError(f"R markdown path is not a file: {rmarkdown_path}")
    if rmarkdown_path.suffix.lower() not in [".rmd"]:
        click.echo(f"R markdown does not have .Rmd extension: {rmarkdown_path}")
        raise ValueError(f"R markdown does not have .Rmd extension: {rmarkdown_path}")
    with cd(script_dir):
        click.echo(Path.cwd())
        rmarkdown_cmd = f"Rscript -e \"rmarkdown::render('{rmarkdown_path}'))\""
        subprocess.run(rmarkdown_cmd, shell=True, check=True)


@click.command()
def cli():
    pass


@click.command("tables")
@click.option(
    "--sourcedir",
    default=str(Path(__file__).parent.parent.parent / "G00x-plots/"),
    show_default=True,
    help="Path to the source directory for the R script.",
)
def tables(sourcedir):
    render_rmarkdown(rmarkdown_path="Schief856_Bcell_tables_for_manuscript.Rmd", sourcedir=sourcedir)


# Create a click group to include both commands
@click.group()
def main():
    """Run R scripts."""
    pass


main.add_command(cli)
main.add_command(comparison_tables)

if __name__ == "__main__":
    main()  # noqa: F403
