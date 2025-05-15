from math import sqrt

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def wilsonify(g, data, x, y, hue=None, z=1.96):
    """
    Adds Wilson confidence interval error bars to a Seaborn pointplot.

    Parameters:
    - g: the Axes object returned by sns.pointplot
    - data: original dataframe used in sns.pointplot
    - x: column name used for x-axis (categorical)
    - y: column name used for binary outcome (0/1)
    - hue: column name for hue groupings (if any)
    - z: z-score for confidence interval (default is 1.96 for 95%)
    """

    def wls_confidence_interval(values, z=1.96):
        n = len(values)
        if n == 0:
            return (np.nan, np.nan)
        p = np.mean(values)
        denominator = 1 + z**2 / n
        centre_adjusted_probability = p + z**2 / (2 * n)
        adjusted_standard_deviation = sqrt((p * (1 - p) + z**2 / (4 * n)) / n)
        lower_bound = (centre_adjusted_probability - z * adjusted_standard_deviation) / denominator
        upper_bound = (centre_adjusted_probability + z * adjusted_standard_deviation) / denominator
        return (lower_bound, upper_bound)

    group_cols = [x]
    if hue:
        group_cols.append(hue)

    # Group data and compute Wilson intervals
    summary = (
        data.groupby(group_cols)[y].agg(["mean", "count", "sum"]).rename(columns={"mean": "proportion"}).reset_index()
    )

    ci_bounds = summary.apply(
        lambda row: wls_confidence_interval(np.repeat([1, 0], [row["sum"], row["count"] - row["sum"]]), z=z), axis=1
    )

    summary["ci_low"] = ci_bounds.apply(lambda x: x[0])
    summary["ci_high"] = ci_bounds.apply(lambda x: x[1])
    summary["yerr_low"] = summary["proportion"] - summary["ci_low"]
    summary["yerr_high"] = summary["ci_high"] - summary["proportion"]

    # Map error bars to the correct x-axis positions
    if hue:
        offsets = g.get_children()[0].get_offsets()
        x_positions = offsets[:, 0]
    else:
        x_positions = np.arange(len(summary))

    for i, (_, row) in enumerate(summary.iterrows()):
        g.errorbar(
            x=x_positions[i],
            y=row["proportion"],
            yerr=[[row["yerr_low"]], [row["yerr_high"]]],
            fmt="none",
            capsize=5,
            color="black",
        )


if __name__ == "__main__":
    # Create example binary data (e.g., success/failure)
    np.random.seed(42)
    n_per_group = 100
    groups = ["A", "B", "C"]
    df = pd.DataFrame(
        {
            "group": np.repeat(groups, n_per_group),
            "success": np.concatenate(
                [
                    np.random.binomial(1, 0.6, n_per_group),
                    np.random.binomial(1, 0.4, n_per_group),
                    np.random.binomial(1, 0.5, n_per_group),
                ]
            ),
        }
    )
    import seaborn as sns

    print(sns)
    # Demonstration with simple binary data
    g = sns.pointplot(data=df, x="group", y="success", errorbar=None)
    wilsonify(g, df, x="group", y="success")
    plt.show()
