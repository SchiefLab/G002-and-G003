import click

from g00x.data import Data, PlotParameters
from g00x.cli import g002, g003
from g00x_figures.cli import figures
from VISC_codebase.cli import comparison_tables


@click.group()
@click.pass_context
def main(ctx: click.Context):
    """Run All scripts."""
    ctx.ensure_object(dict)
    ctx.obj = {"data": Data(), "params": PlotParameters()} 
    #pass


# Supplementary Comparison Tables
figures.add_command(comparison_tables, name="comparison-tables")
# Supplementary Tables
# hardcode paths to internal VISC data; will have to rely on PDF already generated
# Main Figures
main.add_command(figures, name="plot")
# Pipeline
main.add_command(g002, name="g002")
main.add_command(g003, name="g003")

if __name__ == "__main__":
    # Run the CLI
    main()
    click.echo("This is a placeholder for the CLI. Please implement the necessary commands.")
