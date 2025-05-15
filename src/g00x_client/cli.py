import click

from g00x_figures.cli import figures, plot_test
from VISC_codebase.cli import comparisons, test


@click.group()
def main():
    """Run All scripts."""
    pass


main.add_command(figures, name="plot")
main.add_command(comparisons, name="comparisons")
main.add_command(test, name="test")

if __name__ == "__main__":
    # Run the CLI
    main()
    click.echo("This is a placeholder for the CLI. Please implement the necessary commands.")
