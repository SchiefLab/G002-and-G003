import contextlib
import math
import warnings
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import ticker as mtick
from matplotlib.patches import Patch

from g00x_figures.data import Data

mpl.rcParams["font.sans-serif"] = [
    "Arial",
    "Liberation Sans",
    "Bitstream Vera Sans",
]
mpl.rcParams["font.family"] = "sans-serif"
# data = Data()
pal = {
    "G001": "#6C65FF",
    "G002": "#91FCC0",
    "G003": "#E377C2",
}
pallete = pal
week_ticks = ["-5", "4", "8", "10", "16", "24/21"]
week_order = [-5, 4, 8, 10, 16, 21]
yaxis_fontsize = 16
xaxis_fontsize = 16


def calculate_resonse(df: pd.DataFrame) -> pd.DataFrame:
    """calculate response rate given our critera from paper 1"""
    for _, group_df in df.groupby(["ptid"]):
        # get wk 5 first
        wk5 = group_df.query("weeks==-5")

        # fill seqs we didn't get with 0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            wk5["Percent of VRC01-class sequences among IgG"] = wk5[
                "Percent of VRC01-class sequences among IgG"
            ].fillna(0)

        # all wk -5 are 0
        df.loc[wk5.index, "Response x"] = 0

        # do we get percent VRC01-class or percent epitope specific
        try:
            pct_ig = wk5["Percent of VRC01-class sequences among IgG"].iloc[0]
        except:
            pct_ig = 0
        try:
            pct_es = wk5["Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}"].iloc[0]
        except:
            pct_es = 0

        # if we didn't find any vrc01 specific just use pct_es
        if pct_ig == 0:
            value = pct_es
        else:
            value = pct_ig

        # now for the rest of the weeks
        not_group_5 = group_df.query("weeks!=-5")

        # for every week lets calculate
        for index, week in not_group_5.iterrows():
            if week["Percent of VRC01-class sequences among IgG"] == 0:
                df.loc[index, "Response x"] = 0  # type: ignore
            elif week["Percent of VRC01-class sequences among IgG"] > value:
                df.loc[index, "Response x"] = 1  # type: ignore
            else:  # if it equals na i think
                df.loc[index, "Response x"] = 0  # type: ignore
    return df


def get_response_count_dataframe(combined) -> pd.DataFrame:
    # g002_flow = data.get_g002_flow_and_seq_prime()
    # combined = combined.rename(
    #     {"Response (missing seq to 0)": "Response x"}, axis=1
    # )
    # combined = pd.concat([g002_flow, flow_and_seq_merged]).reset_index(drop=True)
    combined["weeks"] = combined["weeks"].astype(int)
    # get rid of pre vaccine weeks
    combined = combined.drop(combined.query("weeks < 0").index)
    # combined = combined.query("weeks > 0")
    combined.loc[combined.query("weeks < 16").index, "shot"] = "first"
    combined.loc[combined.query("weeks > 8").index, "shot"] = "second"

    # create a new dataframe for the response
    new_df = []

    # do it once grouping by shot
    for (shot, trial), group_df in combined.groupby(["shot", "trial"]):
        # group_df = group_df.query("weeks > 0")

        responders = len(group_df.query("`Response x`==1")["pubID"].unique())
        print(responders)
        total = len(group_df["pubID"].unique())
        if trial == "G001":
            total -= 1
        new_df.append(
            {
                "trial": trial,
                "shot": shot,
                "responders": responders,
                "total": total,
                "percentage_responders": responders / total,
            }
        )

    # do it again without grouping by shot
    for trial, group_df in combined.groupby("trial"):
        # group_df = group_df.query("weeks > 0")

        responders = len(group_df.query("`Response x`==1")["pubID"].unique())
        if trial == "G001":
            responders += 1
        total = len(group_df["pubID"].unique())
        new_df.append(
            {
                "trial": trial,
                "shot": "all",
                "responders": responders,
                "total": total,
                "percentage_responders": responders / total,
            }
        )

    return pd.DataFrame(new_df)


def run_percent_igg_responders(data: Data):
    flow_and_seq_merged = data.get_g00x_flow_and_seq()
    set_letters = False
    dfs = []

    @contextlib.contextmanager
    def figure_context(*args, **kwargs):
        fig, axes = plt.subplots(*args, **kwargs)
        yield fig, axes
        plt.close("all")

    with figure_context(
        1,
        2,
        figsize=(6, 3),
        sharey=True,
        gridspec_kw={
            "wspace": 0.1,
            "left": 0.15,
            "right": 0.99,
            "top": 0.83,
        },
    ) as (fig, axes):
        df = get_response_count_dataframe(flow_and_seq_merged)
        # first shot
        ax = axes[0]
        ax.letter = "A"
        pal = {k: v for k, v in pallete.items() if k in df["trial"].unique()}

        _df = df.query("shot=='first'").copy(deep=True)
        sns.barplot(
            x="trial",
            y="percentage_responders",
            data=_df,
            linewidth=1,
            edgecolor="black",
            ax=ax,
            hue="trial",
            order=pal.keys(),
            dodge=False,
            palette=pal,
        )
        _df.name = "first_shot"
        _df.key = ["percentage_responders"]
        dfs.append(_df)
        ax.yaxis.set_major_locator(mtick.MultipleLocator(base=0.1))
        ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: str(int(float(y) * 100))))
        ax.set_ylim(0, 1)
        sns.despine()
        ax.legend_.remove()

        ax = axes[1]
        ax.letter = "A"
        _df = df.query("shot=='all'").copy(deep=True)
        sns.barplot(
            x="trial",
            y="percentage_responders",
            data=_df,
            linewidth=1,
            edgecolor="black",
            ax=ax,
            hue="trial",
            # order=["G002", "G003", "CFHR", "Aurum", "NAC"],
            order=pal.keys(),
            dodge=False,
            palette=pal,
        )
        _df.name = "all_shots"
        _df.key = ["percentage_responders"]
        dfs.append(_df)

        sns.despine()
        ax.set_ylabel("")
        ax.legend_.remove()

        for ax_index, ax in enumerate(axes):
            for x in ax.get_xticklabels():
                (x1, _) = x.get_position()
                lookup = x.get_text()
                if ax_index == 0:
                    value = "first"
                else:
                    value = "all"
                sub_df = df.query("shot==@value").query("trial==@lookup")
                responders = sub_df["responders"].iloc[0]
                total = sub_df["total"].iloc[0]
                percentage = sub_df["percentage_responders"].iloc[0]
                ax.text(
                    x1,
                    percentage + 0.01,
                    f"{responders}/{total}",
                    ha="center",
                    fontsize=xaxis_fontsize,
                )

        # ax.yaxis.set_major_locator(mtick.MultipleLocator(base=0.1))
        # ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: str(int(float(y) * 100))))
        # ax.set_ylim(0, 1)

        axes[0].set_title("After 1 vaccination\n(wk 4,8)", pad=20, fontsize=yaxis_fontsize)  #
        axes[1].set_title(
            "After 1 or 2 vaccinations\n(wk 4,8,10,16,21,24)",
            pad=20,
            fontsize=yaxis_fontsize,
        )
        axes[0].set_ylabel(
            "% VRC01-class\nIgG$^{+}$ responders",
            labelpad=3,
            fontsize=yaxis_fontsize,
        )
        axes[0].set_xlabel("")
        axes[1].set_xlabel("")
        axes[0].tick_params(axis="x", labelsize=xaxis_fontsize)
        axes[0].tick_params(axis="y", labelsize=xaxis_fontsize)
        axes[1].tick_params(axis="x", labelsize=xaxis_fontsize)
        axes[1].tick_params(axis="y", labelsize=xaxis_fontsize)
        if set_letters:
            axes[0].set_title("I", fontsize=24, fontweight="bold", x=-0.2, y=1.1)
        return fig, dfs
        # fig.savefig(
        #     local_out / f"percent_igg_responders_{suffix_name}.png",
        #     bbox_inches="tight",
        #     dpi=1000,
        # )
        # plt.show()
