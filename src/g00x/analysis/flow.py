from pathlib import Path

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

from g00x.data import Data


def make_legend(ax: plt.axis) -> None:
    custom_lines: list[Patch] = []
    for x, label in zip(["#91FCC0", "#E377C2", "#2078B4", "#FF7F0F"], ["1", "2", "3", "4"]):
        custom_lines.append(Patch(facecolor=x, edgecolor="black", linewidth=1, label=label))
        ax.legend(
            custom_lines,
            [i._label for i in custom_lines],
            loc="upper center",
            frameon=False,
            handlelength=0.8,
            ncol=1,
            bbox_to_anchor=(1.1, 1),
            labelspacing=0.1,
            title="Group",
        )


def get_frequency_dataframe(data: Data, count_dataframe: pd.DataFrame) -> pd.DataFrame:
    flow_df_indexes = ["ptid", "group", "weeks", "visit_id", "probe_set", "run_purpose"]
    min_flow_columns = [
        "ptid",
        "group",
        "weeks",
        "visit_id",
        "probe_set",
        "run_purpose",
        "gate",
        "phenotype",
        "value_type",
        "value",
        "branch",
        "easy_name",
    ]

    # get minimal columns and set an index
    flow_df_min = count_dataframe[min_flow_columns].set_index(flow_df_indexes)

    def get_flow_freq(entry_lookup: dict[str, str]) -> pd.DataFrame:
        # numerator gate
        num = entry_lookup["numerator"]

        # denominator gate
        dom = entry_lookup["denominator"]

        # a name for this frequency lookup
        axis_name = entry_lookup["axis_name"]

        return_value = flow_df_min.query(f"gate=='{num}'")["value"] / flow_df_min.query(f"gate=='{dom}'")["value"]
        return_value = return_value.to_frame().rename({"value": f"{axis_name}"}, axis=1)
        return return_value

    flow_frequency_df = pd.concat(
        [get_flow_freq(entry) for entry in data.get_frequency_measures()], axis=1
    ).reset_index()
    return flow_frequency_df


def count_current_samples(data: Data, flow_dataframe: pd.DataFrame, output: Path | str) -> None:
    """Count all the current samples we have so far

    Parameters
    ----------
    data : Data
        The data object housing G00x data
    flow_dataframe : pd.DataFrame
        The flow dataframe from the flow pipeline
    output : Path | str
        The output path to save the count dataframe
    """

    # get the unique samples, just a single gate per sample
    flow_df_indexes = ["run_purpose", "ptid", "group", "weeks", "visit_id", "probe_set"]
    flow_counts = flow_dataframe[flow_df_indexes]
    flow_unique = flow_counts[~flow_counts.duplicated(keep="first")].reset_index(drop=True)

    # change weeks to in case
    flow_unique["weeks"] = flow_unique["weeks"].astype(int)

    # groupby these to get counts
    count_coulumns = ["run_purpose", "group", "weeks", "visit_id", "probe_set"]
    count_df = flow_unique.groupby(count_coulumns).size().to_frame().rename({0: "size"}, axis=1).reset_index()

    fig, axes = plt.subplots(1, 4, figsize=(8, 3), sharey=True)
    for ax_index, ax in enumerate(axes):
        if ax_index < 2:
            purpose = "PreS"
        else:
            purpose = "Sort"
        if ax_index == 0 or ax_index == 2:
            probe_set = "eODGT8"
        else:
            probe_set = "Cg28v2"
        query_df = (
            count_df.query(f"run_purpose=='{purpose}'")
            .query(f"probe_set=='{probe_set}'")
            .pivot(["weeks"], "group", "size")  # type: ignore
        )

        query_df = query_df.reindex(index=count_df["weeks"].unique(), columns=count_df["group"].unique())

        query_df.plot(
            kind="bar",
            stacked=True,
            ax=ax,
            width=0.75,
            linewidth=1,
            edgecolor="black",
            color=["#91FCC0", "#E377C2", "#2078B4", "#FF7F0F"],
        )
        ax.yaxis.set_major_locator(plt.MultipleLocator(2))
        ax.legend_.set_visible(False)
        ax.set_xlabel("Weeks post vaccine")
        ax.set_title(f"{purpose}:{probe_set}")
        if ax_index == 0:
            ax.set_ylabel("Counts")
        if ax_index == 3:
            make_legend(ax)
    sns.despine()
    fig.tight_layout()
    fig.savefig(str(output) + ".png", dpi=300)
