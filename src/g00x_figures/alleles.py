from pathlib import Path

import pandas as pd

from g00x_figures.data import Data


def plot_allele_table(data: Data, outpath: str | Path):
    personalized_alleles = data.get_personalized_alleles()
    g002_seq_and_flow = data.get_g002_flow_and_seq_prime()
    merged_seq_and_allele = g002_seq_and_flow.merge(personalized_alleles, on="ptid", how="left")
    series = []
    for g, g_df in merged_seq_and_allele.groupby("ptid"):
        allele = (g_df["allele_1"] + "/" + g_df["allele_2"]).iloc[0]
        if (g_df["Response x"] == 1).any():
            s = pd.Series({"ptid": g, "global_response": 1, "allele": allele})
        else:
            s = pd.Series({"ptid": g, "global_response": 0, "allele": allele})
        series.append(s)

    allele_table = pd.DataFrame(series)
    allele_table = (
        allele_table.groupby(["allele", "global_response"])
        .size()
        .reset_index()
        .rename({0: "count"}, axis=1)
        .sort_values("count")
        .reset_index(drop=True)
    )
    allele_table.replace({"global_response": {0: "No", 1: "Yes"}}, inplace=True)
    allele_table = allele_table.replace("IGHV1-2\*", "", regex=True)
    allele_table.columns = ["Allele", "Response", "Count"]
    allele_table.to_excel(str(outpath) + ".xlsx", index=False)


def plot_allele_table_group4(data: Data, outpath: str | Path) -> None:
    personalized_alleles = data.get_personalized_alleles()
    g002_seq_and_flow = data.get_pre_post_core_df().query("group==4")
    merged_seq_and_allele = g002_seq_and_flow.merge(personalized_alleles, on="ptid", how="left")
    series = []
    for g, g_df in merged_seq_and_allele.groupby("ptid"):
        allele = (g_df["allele_1"] + "/" + g_df["allele_2"]).iloc[0]
        if (g_df["Response y"] == 1).any():
            s = pd.Series({"ptid": g, "global_response": 1, "allele": allele})
        else:
            s = pd.Series({"ptid": g, "global_response": 0, "allele": allele})
        series.append(s)

    allele_table = pd.DataFrame(series)
    allele_table = (
        allele_table.groupby(["allele", "global_response"])
        .size()
        .reset_index()
        .rename({0: "count"}, axis=1)
        .sort_values("count")
        .reset_index(drop=True)
    )
    allele_table.replace({"global_response": {0: "No", 1: "Yes"}}, inplace=True)
    allele_table = allele_table.replace("IGHV1-2\*", "", regex=True)
    allele_table.columns = ["Allele", "Response", "Count"]
    allele_table.to_excel(str(outpath) + ".xlsx", index=False)
