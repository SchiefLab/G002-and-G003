import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import ticker as mtick

from g00x_figures.box_and_scatter.flow_frequencies import (
    MinorSymLogLocator,
    adjust_boxplot,
)
from g00x_figures.data import Data, Transforms
from g00x_figures.plots import Plot


class FlowCytometryPlot:
    def __init__(self):
        self.data = Data()
        self.transforms = Transforms()
        self.plot = Plot()
        self.lt = 1e-3

        self.columns = {
            "B cells": "B cells",
            "IgD-IgM-/IgA-/IgG+ B cells": "IgG$^{+}$ B cells",
            "IgD-/IgM-/IgA-/IgG+/Antigen++ B cells": "eOD$^{++}$ IgG$^{+}$ B cells",
            "IgD-/IgM-/IgA-/IgG+/KO-/Antigen++ B cells": "eOD CD4bs-specific\nIgG$^{+}$ B cells",
        }
        self.ylims = {0: (1e5, 1e8), 1: (1e3, 1e7), 2: (1, 1e6), 3: (1, 1e5)}

    @staticmethod
    def format_y_axis(y: float, lt: float) -> str:
        """Format Y axis for log scale."""
        if y > lt:
            return "{:g}".format(y)
        if y == 0:
            return "0"
        if y == lt:
            return r"$\leq$" + r"$10^{" + str(round(math.log(y, 10))) + r"}$"
        return r"$10^{" + str(round(math.log(y, 10))) + r"}$"

    def remove_titles(self, fig):
        for ax in fig.get_axes():
            ax.set_title("")
        return fig

    def set_axis_properties(self, ax, i: int, j: int, week: int, col_label: str):
        """Set properties for a single axis"""
        # Set scales and labels
        ax.set_yscale("log")

        # Set titles
        ax.set_title(f"Week {week}" if i == 0 else "")
        ax.yaxis.set_major_locator(mtick.LogLocator(base=10, numticks=10))
        ax.yaxis.set_minor_locator(MinorSymLogLocator(linthresh=1))
        ax.yaxis.set_major_formatter(
            mtick.FuncFormatter(lambda y, _: r"$10^{" + str(int(np.log10(y))) + r"}$" if y > 0 else "0")
        )
        # Set y-axis properties
        if j == 0:
            ax.set_ylabel(col_label)

            # Set y-limits based on row
            if i in self.ylims:
                ax.set_ylim(*self.ylims[i])
        else:
            ax.set_ylabel("")

        # Set x-axis properties
        if i == len(self.columns) - 1:
            ax.set_xlabel("Group")
        else:
            ax.set_xlabel("")
            ax.set_xticklabels([])

    def create_plot(self, x: str, df: pd.DataFrame, palette: dict):
        """Create the main plot"""
        fig, ax = plt.subplots(
            len(self.columns),
            df.weeks.nunique(),
            figsize=(8, 10),
            sharex="col",
            sharey="row",
        )

        dfs = []
        for i, (col, col_label) in enumerate(self.columns.items()):
            for j, week in enumerate(sorted(df.weeks.unique())):
                _df = df.query(f"weeks=={week}")
                _df.name = f"i{i}_j{j}_wk{week}"
                _df.key = [x, col]
                dfs.append(_df)
                # Create strip plot
                self.plot.plot_stripbox(
                    df=_df,
                    x=x,
                    y=col,
                    ax=ax[i, j],
                    palette=palette[week],
                )

                # Set axis properties
                self.set_axis_properties(ax[i, j], i, j, week, col_label)
                adjust_boxplot(ax[i, j])

        sns.despine()

        letters = ["A", "B", "C", "D"]
        for i, (col, col_label) in enumerate(self.columns.items()):
            # ax[i, 0].set_title(letters[i], fontsize=14, fontweight="bold", x=-0.7, y=1)
            ax[i, 0].text(-0.7, 1, letters[i], transform=ax[i, 0].transAxes, size=14, weight="bold")
        return (
            fig,
            dfs,
        )


def main():
    plotter = FlowCytometryPlot()
    fig, dfs = plotter.create_plot()
    plt.show()


if __name__ == "__main__":
    main()
