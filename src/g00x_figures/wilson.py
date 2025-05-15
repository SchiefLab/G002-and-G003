from math import sqrt

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def wilsonify(g, data, x, y, hue=None, z=1.96):
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

    # Compute Wilson intervals
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

    # Determine x-axis positions
    positions = {}
    for line in g.lines:
        x_data = line.get_xdata()
        y_data = line.get_ydata()
        for x_val, y_val in zip(x_data, y_data):
            positions[(x_val, y_val)] = x_val

    # Now match each row in summary to its x-position and draw error bars
    for line in g.lines:
        x_vals = line.get_xdata()
        y_vals = line.get_ydata()
        for x_val, y_val in zip(x_vals, y_vals):
            match = summary[np.isclose(summary["proportion"], y_val)]
            if hue:
                # If hue, filter further by category position
                hue_levels = g.legend_.texts
                for idx, row in match.iterrows():
                    g.errorbar(
                        x=x_val,
                        y=row["proportion"],
                        yerr=[[row["yerr_low"]], [row["yerr_high"]]],
                        fmt="none",
                        capsize=5,
                        color="black",
                    )
            else:
                for idx, row in match.iterrows():
                    g.errorbar(
                        x=x_val,
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
