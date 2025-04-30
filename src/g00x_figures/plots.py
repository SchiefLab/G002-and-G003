# %%
import argparse
import glob
import os
import re
import subprocess
import sys
from collections import defaultdict
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import ticker as mtick
from matplotlib.patches import Patch
from scipy.stats import gmean

from g00x_figures.data import Data, Transforms


class Plot:
    colors = [
        "#91FCC0",
        "#2078B4",
        "#E377C2",
        "#9567BD",
        "#17BFD0",
        "#FF7F0F",
        "#BDBD23",
        "white",
    ]
    explicit_mapping = {
        "IGKV1-33": "#91FCC0",
        "IGKV1-5": "#E377C2",
        "IGKV3-20": "#2078B4",
        "IGLV2-14": "#17BFD0",
        "IGKV3-15": "#9567BD",
        "IGLV2-11": "#BDBD23",
        "IGLV2-23": "#FF7F0F",
        "Other": "white",
    }
    week_palette = {
        -5: "#9567BD",
        4: "#17BFD0",
        8: "gold",
        10: "gold",
        16: "#E377C2",
        21: "#2078B4",
        24: "#2078B4",
    }
    psuedogroups_palette = {
        2.0: "gold",
        3.0: "#E377C2",
        4.0: "#E377C2",
        5.0: "#2078B4",
        6.0: "#2078B4",
        7.0: "#2078B4",
    }
    psuedogroups_legend_palette = {False: "#91FCC0", True: "#6C65FF"}
    fixed_locator_lookup = {
        1e-13: "100 fM",
        1e-12: "1 pM",
        1e-11: "10 pM",
        1e-10: "100 pM",
        1e-09: "1 nM",
        1e-08: "10 nM",
        1e-07: "100 nM",
        1e-06: "1 $\mathregular{\mu}$M",  # noqa: W605
        1e-05: "10 $\mathregular{\mu}$M",  # noqa: W605
        0.0001: "100$\mathregular{\mu}$M",  # noqa: W605
        # 0.0001: "$\mathregular{\geq}$100$\mathregular{\mu}$M",  # noqa: W605
        0.001: "1 mM",
        0.01: "10 mM",
        0.1: "100 mM",
        0: "1 M",
        1: "10 M",
        2: "100 M",
        10: "1 G",
    }

    def get_1x2_grid(self) -> tuple[Any, Any]:
        figure, axes = plt.subplots(
            1,
            2,
            figsize=(7, 4),
            gridspec_kw={
                "width_ratios": [1, 0.3],
                "hspace": 0.12,
                "wspace": 0.1,
                "bottom": 0.25,
                "top": 0.85,
                "left": 0.18,
                "right": 0.9,
            },
        )
        return figure, axes

    def get_1x3_grid(self, size: tuple[float, float] = (9, 6)) -> tuple[Any, Any]:
        figure, axes = plt.subplots(
            1,
            3,
            figsize=size,
            gridspec_kw={
                "width_ratios": [1, 1, 0.1],
                "hspace": 0.12,
                "wspace": 0.1,
                "bottom": 0.25,
                "top": 0.85,
                "left": 0.18,
                "right": 0.9,
            },
        )
        return figure, axes

    def get_1x4_grid(self, size: tuple[float, float] = (9, 6)) -> tuple[Any, Any]:
        figure, axes = plt.subplots(
            1,
            4,
            figsize=size,
            gridspec_kw={
                "width_ratios": [1, 1, 0.2, 0.1],
                "hspace": 0.12,
                "wspace": 0.1,
                "bottom": 0.25,
                "top": 0.85,
                "left": 0.18,
                "right": 0.9,
            },
        )
        return figure, axes

    def get_color_mapping(self, list_of_labels: list[str], colors: list[str] | None = None) -> dict[str, str]:
        colors = colors or self.colors
        assert len(list_of_labels) <= len(
            colors
        ), f"Too many labels {len(list_of_labels)} {list_of_labels} for the number of colors {len(self.colors)} {self.colors}"
        return {label: color for label, color in zip(list_of_labels, colors)}

    def get_week_color_mapping(self, list_of_labels: list[str], expand: bool = False):
        pal = {week: self.week_palette[week] for week in list_of_labels}
        if expand:
            pal = {f"wk {week}": color for week, color in pal.items()}
        return pal

    def side_legend(
        self,
        axes: list[Any],
        palette: dict[str, str],
        axe_loc: int = -1,
        fontsize: int = 8,
        bbox_to_anchor: tuple[float, float] = (0.5, 1),
    ) -> None:
        custom_lines = []
        for label, color in palette.items():
            custom_lines.append(Patch(facecolor=color, edgecolor="black", linewidth=1, label=label))
        axes[axe_loc].legend(
            custom_lines,
            [i._label for i in custom_lines],
            loc="upper center",  # occupies a whole subgrid so top center is just top
            frameon=False,
            handlelength=0.8,
            ncol=1,
            fontsize=fontsize,
            bbox_to_anchor=bbox_to_anchor,
            labelspacing=0.1,
        )
        # axes[axe_loc].set_axis_off()

    def bottom_legend(
        self,
        ax: plt.Axes,
        palette: dict[str, str],
        title: str | None = None,
        distance: float = -0.4,
    ) -> None:
        custom_lines = []
        for label, color in palette.items():
            custom_lines.append(Patch(facecolor=color, edgecolor="black", linewidth=1, label=label))
        ax.legend(
            custom_lines,
            [i._label for i in custom_lines],
            loc="center",
            frameon=False,
            handlelength=0.8,
            ncol=3,
            bbox_to_anchor=(0.5, distance),
            labelspacing=0.1,
            title=title,
        )

    def bottom_legend_random_selected(self, ax: plt.Axes, *args, **kwargs) -> None:
        self.bottom_legend(
            ax,
            palette={"Random": "#91FCC0", "Selected": "#6C65FF"},
            *args,
            **kwargs,
        )

    def stacked_bar_plot_pivot_df(
        self,
        pivot_df: pd.DataFrame,
        title: str,
        ylabel: str | None = None,
        xlabel: str | None = None,
        yaxis_mtick_major: float | None = None,
        yaxis_mtick_minor: float | None = None,
        remove_yticklabels: bool = False,
        xticklabels: list[str] | None = None,
        xticklabel_rotation: int = 0,
        remove_legend: bool = True,
        yaxis_mtick_normalize: bool = True,
        colors: list[str] | None = None,
        fig: plt.Figure | None = None,
        ax: plt.Axes | None = None,  # type: ignore
        *args: Any,  # type: ignore
        **kwargs: Any,  # type: ignore
    ) -> plt.Axes:  # -> Any | Axes | None:
        """stacked bar plot from a pivot dataframe.
        should be limited indexes and columns to keep size down.

        Parameters
        ----------
        pivot_df : pd.DataFrame
            dataframe with multiindex and columns
        title : str
            title of the plot
        ylabel : str | None, optional
            yaxis label, by default None
        xlabel : str | None, optional
            xaxis label, by default None
        yaxis_mtick_major : float | None, optional
            physical yaxis tick - combines with minor, by default None
        yaxis_mtick_minor : float | None, optional
            physical yaxis tick - combines with major, by default None
        remove_yticklabels : bool, optional
            remove values for each ytick, by default False
        xticklabels : list[str] | None, optional
            xtick values, by default None
        xticklabel_rotation : int, optional
            xtick value rotation - 0 is horizonal; 90 is vertical, by default 0
        remove_legend : bool, optional
            removes legen incase you have a shared legend you need to custom make later, by default True
        yaxis_mtick_normalize : bool, optional
            i.e. make ytick value .5 to 50 if False, by default True
        colors : list[str] | None, optional
            hex or color name list with exact number of elements as xaxis, by default None
        ax : plt.Axes | None, optional
            matplotlib graph index location to edit specific subgraph, by default None

        Returns
        -------
        plt.Axes
            matplotlib graph index location to edit specific subgraph
        >>> stacked_bar_plot_pivot_df(
                pivot_df=g00x_seq_prime_igg_vrc01_class_pivot_trial_df,
                title="VRC01-class",
                ylabel="% BCRs with VRC01-class \n" + r"bnAb $\mathregular{V_{K/L}}$",
                yaxis_mtick_major=0.2,
                yaxis_mtick_minor=0.1,
                yaxis_mtick_normalize=False,
                colors=v_call_top_gene_color_mapping.values(),
                ax=axes[0],)
        >>> stacked_bar_plot_pivot_df(
                pivot_df=deskosky_ready,
                title="Control",
                xticklabels=["Dekosky\nVH1-2"],
                remove_yticklabels=True,
                yaxis_mtick_major=0.2,
                yaxis_mtick_minor=0.1,
                colors=v_call_top_gene_color_mapping.values(),
                ax=axes[1],)
        """
        ax: plt.Axes = pivot_df.plot(
            kind="bar",
            stacked=True,
            color=colors,
            width=0.85,
            linewidth=1,
            edgecolor="black",
            ax=ax,
            *args,
            **kwargs,
        )
        ysize = None
        xsize = None
        # if fig:
        #     size = fig.get_size_inches()
        #     ysize = size[1] * 2
        #     xsize = size[0] * 0.75

        ax.set_title(title, fontsize=ysize or None)
        ax.set_ylabel(ylabel or "", fontsize=ysize or None)
        ax.set_xlabel(xlabel or "", fontsize=xsize or None)
        ax.set_xticklabels(
            xticklabels or ax.get_xticklabels(),
            rotation=xticklabel_rotation,
            fontsize=xsize or None,
            ha="right" if xticklabel_rotation else "center",
            rotation_mode="anchor",
        )
        if remove_legend:
            ax.legend_.remove()  # type: ignore
        if yaxis_mtick_major:
            ax.yaxis.set_major_locator(mtick.MultipleLocator(yaxis_mtick_major))
        if yaxis_mtick_minor:
            ax.yaxis.set_minor_locator(mtick.MultipleLocator(yaxis_mtick_minor))
        if remove_yticklabels:
            ax.set_yticklabels([])
        if not yaxis_mtick_normalize:
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
        sns.despine()
        return ax

    def plot_stripbox(
        self,
        df: pd.DataFrame,
        x: str,
        y: str,
        hue: str | None = None,
        title: str | None = None,
        xlabel: str | None = None,
        ylabel: str | None = None,
        xticklabels: list[str] | None = None,
        yaxis_mtick_normalize: bool = False,
        yaxis_mtick_major: float | None = None,
        yaxis_mtick_minor: float | None = None,
        remove_legend: bool = False,
        remove_yticklabels: bool = False,
        remove_xticklabels: bool = False,
        palette: dict[str, str] | None = None,
        dodge: bool = True,
        tilt: bool = True,
        ax: plt.Axes | None = None,  # type: ignore
        *args: Any,
        **kwargs: Any,
    ) -> tuple[plt.Axes, pd.DataFrame]:
        """Dynamic stripplot with boxplot overlay.

        Parameters
        ----------
        df : pd.DataFrame
            Full dataframe
        x : str
            xaxis column name
        y : str
            yaxis column name
        hue : str | None, optional
            hue column name
        xlabel : str | None, optional
            xaxis label, by default None
        ylabel : str | None, optional
            yaxis label, by default None
        xticklabels : list[str] | None, optional
            xaxis tick labels, by default None
        yaxis_mtick_normalize : bool, optional
            yaxis tick % if False decimal if True, by default False
        palette : dict[str, str] | None, optional
            label:color hash, by default None

        Returns
        -------
        plt.Axes
            matplotlib graph index location to edit specific subgraph
        """
        palette = palette or None
        color = None
        if isinstance(palette, str):
            color = palette
            palette = None
        ax = sns.stripplot(
            x=x,
            y=y,
            hue=hue,
            edgecolor="black",
            linewidth=1,
            # s=7,
            jitter=0.15,
            dodge=dodge,
            palette=palette,
            color=color,
            data=df,
            ax=ax,
            *args,
            **kwargs,
        )
        ax = sns.boxplot(
            x=x,
            y=y,
            hue=hue,
            dodge=dodge,
            fliersize=0,
            palette=palette,
            color=color,
            whis=[10, 90],
            data=df,
            ax=ax,
            *args,
            **kwargs,
        )
        if title:
            ax.set_title(title)
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        if remove_yticklabels:
            ax.set_yticklabels([])
        if remove_xticklabels:
            ax.set_xticklabels([])
        if xticklabels:
            if tilt:
                _ = ax.set_xticklabels(
                    xticklabels,
                    rotation=45,
                    ha="right",
                    rotation_mode="anchor",
                    fontsize=8,
                )
            else:
                _ = ax.set_xticklabels(
                    xticklabels,
                    verticalalignment="top",
                    horizontalalignment="center",
                    fontsize=8,
                )
        if remove_legend:
            ax.legend().set_visible(False)  # type: ignore
        if not yaxis_mtick_normalize:
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals=0, symbol=None))
        if yaxis_mtick_major:
            ax.yaxis.set_major_locator(mtick.MultipleLocator(yaxis_mtick_major))
        if yaxis_mtick_minor:
            ax.yaxis.set_minor_locator(mtick.MultipleLocator(yaxis_mtick_minor))
        sns.despine(ax=ax)
        return ax, df

    def plot_stacked_bar_g00x_seq_prime_igg_vrc01_class_pivot_trial_light_value_counts_df(
        self,
        v_call_top_light_genes: list[str],
        transforms: Transforms | None = None,
        ax: Any | None = None,
    ) -> plt.Axes:
        plot = Plot()
        transforms = transforms or Transforms()
        g00x_seq_prime_igg_vrc01_class_pivot_trial_light_value_counts_df = (
            transforms.get_g00x_seq_prime_igg_vrc01_class_pivot_trial_light_value_counts_df()
        )
        v_call_top_gene_color_mapping = plot.get_color_mapping(list_of_labels=v_call_top_light_genes)
        ax = plot.stacked_bar_plot_pivot_df(
            pivot_df=g00x_seq_prime_igg_vrc01_class_pivot_trial_light_value_counts_df,
            title="VRC01-class",
            ylabel="% BCRs with VRC01-class \n" + r"bnAb $\mathregular{V_{K/L}}$",
            yaxis_mtick_major=0.2,
            yaxis_mtick_minor=0.1,
            yaxis_mtick_normalize=False,
            colors=list(v_call_top_gene_color_mapping.values()),
            ax=ax,
        )
        sns.despine()
        return ax

    def plot_stacked_bar_dekosky_pivot_replicate_donors_by_light_value_counts_df(
        self,
        v_call_top_light_genes: list[str],
        transforms: Transforms | None = None,
        ax: Any | None = None,
        seperate_legend: bool = False,
    ) -> plt.Axes:
        plot = Plot()
        transforms = transforms or Transforms()
        dekosky_pivot_replicate_donors_by_light_value_counts_df = (
            transforms.get_dekosky_pivot_replicate_donors_by_light_value_counts_df()
        )
        ylabel = ("% BCRs with VRC01-class \n" + r"bnAb $\mathregular{V_{K/L}}$") if ax is None else None
        v_call_top_gene_color_mapping = plot.get_color_mapping(list_of_labels=v_call_top_light_genes)
        ax = plot.stacked_bar_plot_pivot_df(
            pivot_df=dekosky_pivot_replicate_donors_by_light_value_counts_df,
            title="Control",
            ylabel=ylabel,
            xticklabels=["Dekosky\nVH1-2"],
            remove_yticklabels=True if ax is not None else False,
            yaxis_mtick_major=0.2,
            yaxis_mtick_minor=0.1,
            colors=v_call_top_gene_color_mapping.values(),
            remove_legend=True if ax is not None else False,
            ax=ax,
        )
        if seperate_legend:
            plot.legend(axes=[ax], color_mapping=v_call_top_gene_color_mapping)
        sns.despine()
        return ax

    def get_table(
        self,
        analyte: str = "core-Hx_r4.0D_TH6_g28v2_pHLsecAvi",
        is_vrc01_class: bool = True,
        lim: float = 1e-4,
    ):
        round = 1
        spr_df = Data().get_g002_spr_df_boost()
        spr_df = spr_df.query(f"is_vrc01_class=={is_vrc01_class}").query("top_c_call=='IGHG'")
        spr_df = spr_df.sort_values("Chi2").groupby(["Ligand", "Analyte"]).head(1)
        spr_df = spr_df.query("Analyte==@analyte")
        spr_df.loc[spr_df[spr_df["estimated"] == " TRUE"].index, "KD_fix"] = lim
        spr_df.loc[spr_df[spr_df["KD_fix"] > lim].index, "KD_fix"] = lim

        n_abs_tested = spr_df.groupby(["pseudogroup", "is_cp"]).size().to_frame("n Abs tested").astype(int).astype(str)

        kd_name = "% $\mathregular{K_D}$s < 100 $\mu$M"
        if lim == 5e-5:
            kd_name = "% $\mathregular{K_D}$s < 50 $\mu$M"
        spr_df["kd_lt_100m"] = spr_df["KD_fix"] < lim
        kd_lt_100um = (
            (spr_df.groupby(["pseudogroup", "is_cp"]).kd_lt_100m.mean().to_frame(kd_name) * 100).astype(int).astype(str)
        )

        median_kd_um = (spr_df.groupby(["pseudogroup", "is_cp"]).KD_fix.agg(["median"]) * 1e9).round(round)
        median_kd_um["median"] = median_kd_um["median"].apply(lambda x: str(x) if x < 1 else str(int(x)))
        median_kd_um["median"] = median_kd_um["median"].apply(lambda x: "-" if (float(x) * 1e-9) >= lim else x)
        median_kd_um.rename(columns={"median": "Median $\mathregular{K_D}$ (nM)"}, inplace=True)

        n_participants = (
            spr_df.groupby(["pseudogroup", "is_cp"]).ptid.nunique().to_frame("N participants").astype(int).astype(str)
        )

        df = pd.concat([n_abs_tested, kd_lt_100um, median_kd_um, n_participants], axis=1).T
        return df

    def get_dyn_table(
        self,
        spr_df: pd.DataFrame,
        groupby: list[str],
        lim: float = 1e-4,
        median_scale: float = 1e9,
        geomean: bool = False,
        KD_fix: str = "KD_fix",
        add_kd_lt: bool = True,
        add_binder_median: bool = False,
        median_name: str | None = None,
        sci_notation: bool = False,
    ):
        def get_stats(df, kd_lookup="KD_fix"):
            """Get table stats, counts, median and donors represented"""
            _count = int(len(df))
            less_than_100uM = (len(df[df[kd_lookup] < 1e-4]) / _count) * 100
            median = df[kd_lookup].median()
            subjects = int(len(df["ptid"].unique()))
            return pd.Series(
                {
                    "count": _count,
                    "greater_than_100": less_than_100uM,
                    "median": median,
                    "num_subjects": subjects,
                }
            )

        gt_median = "-"  # r"$\mathregular{-}$"

        def scale_median(
            x, median_scale: float = 1e9, fr: bool = False, skip: bool = False, catch_zero: bool = False
        ) -> str:
            if skip:
                return gt_median
            if np.isnan(x):
                return gt_median
            # if catch_zero and x <= lim:
            #     return gt_median
            if x != lim:
                # print(median_scale)
                i = median_scale * x
                if median_scale == 1:
                    return "{:.1e}".format(i).replace("e-0", "e-").replace("e+0", "e+").replace("e0", "e")
                if i < 1 or fr:
                    return "{:.1e}".format(i).replace("e-0", "e-").replace("e+0", "e+").replace("e0", "e")
                    # return str(round(i, 3))
                else:
                    return str(int(i))
            else:
                return gt_median

        n_abs_tested = spr_df.groupby(groupby).size().to_frame("n Abs tested").astype(int).astype(str)

        kd_name = "% $\mathregular{K_D}$s < 100 $\mu$M"
        if lim == 5e-5:
            kd_name = "% $\mathregular{K_D}$s < 50 $\mu$M"

        # TODO: if uM set to 1e6
        if median_scale == 1e6:
            scale_name = "$\mathregular{\mu}$M"
        elif median_scale == 1e9:
            scale_name = "nM"
        elif median_scale == 1e3:
            scale_name = "mM"
        elif median_scale == 1:
            scale_name = "M"
        else:
            ValueError("median_scale must be 1e6, 1e9 or 1")

        if add_kd_lt:
            spr_df["kd_lt_100m"] = spr_df[KD_fix] < lim
            kd_lt_100um = (spr_df.groupby(groupby).kd_lt_100m.mean().to_frame(kd_name) * 100).astype(int).astype(str)

        median_kd_um = spr_df.groupby(groupby)[KD_fix].agg(["median"])  # .round(1)
        # anything small as one decimal place is a zero
        median_kd_um["median"] = median_kd_um["median"].apply(scale_median, median_scale=median_scale)
        median_name = median_name or "Median $\mathregular{K_D}$" + f" ({scale_name})"
        median_kd_um.rename(columns={"median": median_name}, inplace=True)

        # for each cell run format_cell
        def format_cell(cell):
            def convert_numeric_string(value):
                try:
                    # Try converting to int first
                    return int(value)
                except (ValueError, TypeError):
                    try:
                        # If int fails, try float
                        return float(value)
                    except (ValueError, TypeError):
                        # If both conversions fail, return original value
                        return value

            cell = convert_numeric_string(cell)
            if isinstance(cell, (int, float)):
                cell = "{:.1e}".format(cell).replace("e+", "e")
                p, s = cell.split("e")
                if "-" in s:
                    sep = "e-"
                    s = s[1:]
                else:
                    sep = "e"
                v = p + sep + str(int(s))
                return v.replace("e0", "")
                # remove any zeros
            return cell

        # breakpoint()
        # numeric_columns = df.select_dtypes(include=["int64", "float64"]).columns
        # print(median_kd_um.columns)
        if sci_notation:
            for col in median_kd_um.columns:
                median_kd_um[col] = median_kd_um[col].apply(format_cell)

        binder_median_kd_um = spr_df[spr_df[KD_fix] < lim].groupby(groupby)[KD_fix].agg(["median"])
        binder_median_kd_um["median"] = binder_median_kd_um["median"].apply(scale_median, median_scale=1e6, fr=True)
        binder_median_kd_um.rename(
            columns={"median": "Binder Median $\mathregular{K_D}$" + f" ({scale_name})"},
            inplace=True,
        )

        n_participants = spr_df.groupby(groupby).pubID.nunique().to_frame("N participants").astype(int).astype(str)

        skip = False
        spr_df["KD_fix_lim"] = spr_df[KD_fix]  # .apply(lambda x: x if x < lim else np.nan)

        if geomean:
            # spr_df.loc[spr_df[spr_df["KD_fix_lim"] >= lim].index, "KD_fix_lim"] = np.nan
            df = spr_df[spr_df["KD_fix_lim"] < lim].copy(deep=True)
            # df = spr_df.copy(deep=True)
            if df.empty:
                df["KD_fix_lim"] = lim
                skip = True
                # breakpoint()
            _geomean = (
                df.groupby(groupby)
                .KD_fix_lim.agg(gmean)
                .apply(scale_median, median_scale=median_scale, fr=True, skip=skip, catch_zero=True)
                .to_frame("Binder geomean $\mathregular{K_D}$" + f" ({scale_name})")
            ).fillna(gt_median)

            metrics = [
                n_abs_tested,
                kd_lt_100um,
                n_participants,
                median_kd_um,
                _geomean,
            ]
        elif add_kd_lt:
            metrics = [n_abs_tested, kd_lt_100um, n_participants, median_kd_um]
        else:
            metrics = [n_abs_tested, n_participants, median_kd_um]

        if add_binder_median:
            metrics = metrics + [binder_median_kd_um]

        df = pd.concat(
            metrics,
            axis=1,
            join="outer",
        ).T.fillna(gt_median)

        # Ensure all columns are present in the final dataframe
        # all_columns = set()
        # for metric in metrics:
        #     all_columns.update(metric.index)

        # for metric in metrics:
        #     for row in all_columns:
        #         if row not in metric.index:
        #             metric.loc[row] = 0

        # df = pd.concat(
        #     metrics,
        #     axis=1,
        #     join="inner",
        # ).T
        # breakpoint()

        return df

    def plot_dyn_table(
        self,
        spr_df: pd.DataFrame,
        groupby: list[str],
        # analyte: str,
        ax: plt.Axes | None = None,
        remove_label: bool = False,
        lim: float = 1e-4,
        fontsize: int | None = 8,
        median_scale: float = 1e9,
        geomean: bool = False,
        KD_fix: str = "KD_fix",
        add_kd_lt: bool = True,
        add_binder_median: bool = False,
        bbox: list[float] | None = None,
        median_name: str | None = None,
        colWidths: list[float] | None = None,
        ingect_spaces: list[int] | None = None,
        fill_na: bool = False,
        loc: str = "bottom",
        pad: int = 5,
        sci_notation: bool = False,
    ) -> tuple[plt.Axes, pd.DataFrame]:
        df = spr_df.copy()
        # df.fillna(lim, inplace=True)
        df = self.get_dyn_table(
            spr_df=df[~df[KD_fix].isna()],
            groupby=groupby,
            lim=lim,
            median_scale=median_scale,
            geomean=geomean,
            KD_fix=KD_fix,
            add_kd_lt=add_kd_lt,
            add_binder_median=add_binder_median,
            median_name=median_name,
            sci_notation=sci_notation,
        )
        # from IPython import embed

        # embed()

        if fill_na:
            weeks_by_group = {}

            try:
                for week, group in df.columns:
                    if week not in weeks_by_group:
                        weeks_by_group[week] = []
                    weeks_by_group[week].append(group)

                # Make sure all groups exist for each week
                groups = set(["G001", "G002", "G003"])
                result = []
                for week in sorted(weeks_by_group.keys()):
                    existing_groups = weeks_by_group[week]
                    for group in groups:
                        if group not in existing_groups:
                            result.append((week, group))
            except:
                # same but for single index
                groups = set(["G001", "G002", "G003"])
                result = []
                for group in groups:
                    if group not in df.columns:
                        result.append(group)
            # from IPython import embed

            # embed()
            try:
                for r in result:
                    df[r] = [" ", " ", " ", " "]
            except:
                for r in result:
                    df[r] = [" ", " ", " ", " ", " "]
            df = df.sort_index(axis=1, level=[0, 1])
        # breakpoint()
        # df = df.applymap(lambda x: f"{x:>{pad}}")

        ax = ax or plt.subplot(111, frame_on=False)
        ax.xaxis.set_visible(False)  # hide the x axis
        ax.yaxis.set_visible(False)  # hide the y axis

        ax.set_axis_off()
        ax.grid(False)

        table = ax.table(
            cellText=df.applymap(lambda x: f"{x:>{pad}}").values,
            cellLoc="center",
            rowLoc="right",
            rowLabels=df.index if not remove_label else None,
            edges="open",
            loc=loc,
            colWidths=colWidths,
            bbox=bbox,
        )
        # table.scale(1, 0.5)  # Adjust the scale to remove padding/margins
        # table.auto_set_font_size(False)
        # table.set_fontsize(fontsize or 10)
        # table.auto_set_font_size(True)
        if fontsize:
            table.set_fontsize(fontsize)

        # fix place holder
        df = df.replace("-", lim)
        return table, df

    def plot_table(
        self,
        analyte: str = "core-Hx_r4.0D_TH6_g28v2_pHLsecAvi",
        is_vrc01_class: bool = True,
        ax: plt.Axes | None = None,
        lim: float = 1e-4,
    ) -> plt.Axes:
        df = self.get_table(analyte, is_vrc01_class=is_vrc01_class, lim=lim)
        ax = ax or plt.subplot(111, frame_on=False)
        ax.xaxis.set_visible(False)  # hide the x axis
        ax.yaxis.set_visible(False)  # hide the y axis
        ax.set_axis_off()
        ax.grid(False)
        table = ax.table(
            cellText=df.values,
            cellLoc="center",
            rowLoc="right",
            rowLabels=df.index,
            edges="open",
            loc="bottom",
            # bbox=[0, 1, 1, 0.25],
        )
        # table.auto_set_font_size(False)
        table.set_fontsize(10)
        return table
