import math

import matplotlib.pyplot as plt


def format_y_axis(y: float, lt: float):
    """Format Y axis for log scale."""
    if y > lt:
        return "{:g}".format(y)
    else:
        if y == 0:
            return "0"
        if y == lt:
            return r"$\leq$" + r"$10^{" + str(round(math.log(y, 10))) + r"}$"
        else:
            return r"$10^{" + str(round(math.log(y, 10))) + r"}$"


def adjust_boxplot(ax: plt.Axes) -> plt.Axes:
    """For box and whisker plots, adjust the appearance of the boxplot.

    Parameters
    ----------
    ax : plt.axis

    Returns
    -------
    plt.axis
        return back axis
    """
    # for boxplots we can make adjustments
    list_of_groups = zip(*(iter(ax.get_lines()),) * 6)
    for line_set in list_of_groups:
        for line_index, line in enumerate(line_set):
            # the median line of the box plot
            if line_index == 4:
                line.set_color("black")
                line.set_linewidth(2)
                line.set_solid_capstyle("butt")

    for mybox in ax.artists:
        # Change the appearance of the surrounding boxes, like the sides
        mybox.set_edgecolor("black")
        mybox.set_linewidth(0.4)
    return ax
