import click

from g00x_figures.cli import figures
from VISC_codebase.cli import comparison_tables


@click.group()
def main():
    """Run All scripts."""
    pass


# Supplementary Comparison Tables
figures.add_command(comparison_tables, name="comparison-tables")
# Supplementary Tables
# TODO:
# Main Figures
main.add_command(figures, name="plot")
# Pipeline
# main.add_command(pipeline, name="pipeline")

if __name__ == "__main__":
    # Run the CLI
    main()
    click.echo("This is a placeholder for the CLI. Please implement the necessary commands.")
