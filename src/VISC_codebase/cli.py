import subprocess
from pathlib import Path

import click


def run_r_script(rscript_path: str, sourcedir: str = "../../G00x-plots/") -> None:
    """Run an R script with the given source directory."""
    rscript = Path(rscript_path)
    if not rscript.exists():
        raise FileNotFoundError(f"R script not found: {rscript_path}")
    if not rscript.is_file():
        raise ValueError(f"R script path is not a file: {rscript_path}")
    if rscript.suffix.lower() != ".r":
        raise ValueError(f"R script does not have .R extension: {rscript_path}")
    rscript_cmd = f"Rscript {rscript_path} '{sourcedir}'"
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


# will run the script and then the markdown
# def run_Rscript(rscript_path: str = None, rmarkdown_path: str = None, sourcedir: str = '../../G00x-plots') -> None:
#     # Fail if both rscript_path and rmarkdown_path are None
#     if rscript_path is None and rmarkdown_path is None:
#         raise ValueError("At least one of rscript_path or rmarkdown_path must be provided")

#     # Run R script if provided
#     if rscript_path:
#         run_r_script(rscript_path, sourcedir)

#     # Run R markdown if provided
#     if rmarkdown_path:
#         render_rmarkdown(rmarkdown_path, sourcedir)


@click.command()
# @click.argument('rscript_path')
# @click.argument('rmarkdown_path', required=False)
# @click.option('--sourcedir', default='../../G00x-plots', show_default=True, help='Path to the source directory for the R script.')
# def cli(rscript_path, rmarkdown_path, sourcedir):
#     run_Rscript(rscript_path, rmarkdown_path, sourcedir)
def cli():
    pass


@click.command("comparisons")
@click.option(
    "--sourcedir", default="../../G00x-plots/", show_default=True, help="Path to the source directory for the R script."
)
def comparisons(sourcedir):
    """Run the mk_table_figure2.R script."""
    # render_rmarkdown(rmarkdown_path='FigS8_TabS31-34_TabS77.Rmd', sourcedir=sourcedir)
    # render_rmarkdown(rmarkdown_path='Figure_S10.Rmd', sourcedir=sourcedir)
    # render_rmarkdown(rmarkdown_path='TableS35_S36.Rmd', sourcedir=sourcedir)
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


@click.command("test")
@click.option(
    "--sourcedir", default="../../G00x-plots/", show_default=True, help="Path to the source directory for the R script."
)
def test(sourcedir):
    """Run the mk_table_figure2.R script."""
    # render_rmarkdown(rmarkdown_path='TableS40.Rmd', sourcedir=sourcedir)
    run_r_script(rscript_path="comparisons_table_S62.R", sourcedir=sourcedir)


# Create a click group to include both commands
@click.group()
def main():
    """Run R scripts."""
    pass


main.add_command(cli)
main.add_command(comparisons)
main.add_command(test)


if __name__ == "__main__":
    main()  # noqa: F403
