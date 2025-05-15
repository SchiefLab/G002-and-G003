import re
import warnings
from ast import literal_eval
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pandera.typing import Series
from pydantic import BaseModel, validator


def get_fraction_eq(s: Series[Any], col: str, eq: Any, name: str) -> Series[float]:
    l = len(s[s[col] == eq])  # type: ignore
    return pd.Series({name or "fraction_eq_to": l / len(s)})  # type: ignore


def get_fraction_gt(s: Series[Any], col: str, gt: Any, name: str) -> Series[float]:
    l = len(s[s[col] > gt])  # type: ignore
    return pd.Series({name or "fraction_gt_to": l / len(s)})  # type: ignore


def get_fraction_lt(s: Series[Any], col: str, lt: Any, name: str) -> Series[float]:
    l = len(s[s[col] < lt])  # type: ignore
    return pd.Series({name or "fraction_lt_to": l / len(s)})  # type: ignore


def get_better_than(s, gt: int):
    """Get the fraction of sequences that are better than the Cottrell et al. threshold"""
    s = s.query("top_c_call=='IGHG'")
    if s.empty:
        return
    l = len(s[s["cottrell_focused_v_common_score"] > gt])
    return pd.Series({"fraction_of_cottrell_gt_n": l / len(s)})


def get_better_than_hcdr2(s, gt: int):
    """Get the fraction of sequences that are better than the Cottrell et al. threshold"""
    s = s.query("top_c_call=='IGHG'")
    if s.empty:
        return
    l = len(s[s["num_hcdr2_mutations"] > gt])
    return pd.Series({"fraction_of_cottrell_gt_hcdr2_n": l / len(s)})


def calculate_resonse_boost(df: pd.DataFrame) -> pd.DataFrame:
    """calculate response rate given our critera from paper 1"""
    for index, data in df.iterrows():
        if data["Percent of VRC01-class sequences among IgG"] > 0:  # add sub population to filter
            df.loc[index, "Response y"] = 1  # type: ignore
        elif data["Percent of VRC01-class sequences among IgG"] == 0:
            df.loc[index, "Response y"] = 0  # type: ignore
        else:
            df.loc[index, "Response y"] = 0  # type: ignore
    return df


def calculate_resonse(df: pd.DataFrame) -> pd.DataFrame:
    """calculate response rate given our critera from paper 1"""
    for name, group_df in df.groupby("ptid"):
        # get wk 5 first
        wk5 = group_df.query("weeks==-5")
        resp_using_es = False

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
        # CD4bs
        try:
            pct_es = wk5["Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}"].iloc[0]
        except:
            pct_es = 0

        # if we didn't find any vrc01 specific just use pct_es
        if pct_ig == 0:
            value = pct_es
            resp_using_es = True
        else:
            value = pct_ig
        # print(name, resp_using_es, pct_ig, pct_es, value)

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


def determine_if_vrc01(df: pd.DataFrame):
    # has vh12
    df["has_vh12"] = False
    df.loc[df[df["v_call_heavy"].str.contains(r"IGHV1-2\*")].index, "has_vh12"] = True

    df["has_vh12_02_04"] = False
    df.loc[
        df[df["v_call_heavy"].str.contains(r"IGHV1-2\*0[24]")].index,
        "has_vh12_02_04",
    ] = True

    df["has_vh12_05_06"] = False
    df.loc[
        df[df["v_call_heavy"].str.contains(r"IGHV1-2\*0[56]")].index,
        "has_vh12_05_06",
    ] = True

    # has five leng
    df["has_5_len"] = False
    df.loc[df[df["cdr3_aa_light"].str.len() == 5].index, "has_5_len"] = True

    df["is_vrc01_class"] = False
    df.loc[df.query("has_vh12").query("has_5_len").index, "is_vrc01_class"] = True

    df["is_vrc01_class_02_04"] = False
    df.loc[
        df.query("has_vh12_02_04").query("has_5_len").index,
        "is_vrc01_class_02_04",
    ] = True

    df["is_vrc01_class_05_06"] = False
    df.loc[
        df.query("has_vh12_05_06").query("has_5_len").index,
        "is_vrc01_class_05_06",
    ] = True

    # remove 05/06 from vrc01_class
    df.loc[df.query("has_vh12_05_06").query("has_5_len").index, "is_vrc01_class"] = False

    return df


class DataPaths(BaseModel):
    data_base_path: Path = Path(__file__).parent
    g001_sequences: Path = data_base_path / Path("g001/mutational_analysis.feather")
    g001_clustered_sequences: Path = data_base_path / Path("g001/cluster_df_w_len.feather")
    g001_flow_and_seq: Path = data_base_path / Path("g001/flow_and_seq.csv")
    # g002_flow_and_seq: Path = data_base_path / Path("g002/flow_and_sequencing_long_names_complete.feather")
    g002_flow_and_seq: Path = data_base_path / Path("g002/flow_and_sequencing_long_names_complete.csv")
    g002_flow_and_seq_gates: Path = data_base_path / Path("g002/flow_and_sequencing_long_form.feather")
    g002_sequences: Path = data_base_path / Path("g002/final_df_complete.feather")
    # g002_sequences: Path = data_base_path / Path("g002/final_df_jan24.feather")
    g002_cluster: Path = data_base_path / Path(
        "g002/g002-clusters/final_df_BOOST_pad_somaticFalse_CDR3_link_average_group_v_call-pubID_cluster_n5.feather"
    )
    g002_cluster_prime_nowk24: Path = data_base_path / Path(
        "g002/g002-clusters/final_df_prime_pad_somaticFalse_CDR3_link_average_group_v_call_cluster_n3_no_week24.feather"
    )
    g003_cluster_prime_nowk21: Path = data_base_path / Path(
        "g003/g003-clusters/final_df_prime_pad_somaticFalse_CDR3_link_average_group_v_call_cluster_n3_no_week21.feather"
    )
    g002_cluster_prime: Path = data_base_path / Path(
        "g002/g002-clusters/final_df_prime_pad_somaticFalse_CDR3_link_average_group_v_call_cluster_n3.feather"
    )
    g003_cluster_prime: Path = data_base_path / Path(
        "g003/g003-clusters/final_df_prime_pad_somaticFalse_CDR3_link_average_group_v_call_cluster_n3.feather"
    )
    g003_merged = data_base_path / Path("g003/flow_manifest.feather")
    g003_flow_and_seq: Path = data_base_path / Path("g003/flow_and_sequencing_long_names.feather")
    g003_sequences: Path = data_base_path / Path("g003/final_df.feather")
    g003_methodology_flow_and_seq = data_base_path / Path("g001/flow_and_sequencing_long_names_g001_mg003.csv")
    g003_methodology_seq = data_base_path / Path("g001/final_df_g001_mg003.feather")
    dekosky_vh12: Path = data_base_path / Path("controls/dekosky_vh12.feather")
    oas_5_len: Path = data_base_path / Path("controls/five_length_oas_sample.feather")
    oas_vh12: Path = data_base_path / Path("controls/oas_vh12_sampled.feather")
    vrc01_class_bnabs_path: Path = data_base_path / Path("controls/vrc1-2_mabs_extended_airr_new.feather")
    human_naive_path: Path = data_base_path / Path("controls/human_naive_airr.feather")
    select_kappa_vh12_path: Path = data_base_path / Path("controls/Selected_VRC01c_kappa_chains_segment.feather")
    select_lambda_vh12_path: Path = data_base_path / Path("controls/Selected_VRC01c_lambda_chains_segment.feather")
    personalize_path: Path = data_base_path / Path("g002/personalize.csv")
    g003_personalize_path: Path = data_base_path / Path("g003/genotype.csv")
    # g003_spr: Path = data_base_path / Path("g003/spr_df_g003_w_airr.feather")
    g003_spr: Path = data_base_path / Path("g003/g003_spr_v2.feather")
    g002_spr: Path = data_base_path / Path("g002/spr_df_g002_w_airr.feather")
    # g001_spr_w_id: Path = data_base_path / Path(
    #     "g001/spr_df_g001_w_airr.feather"
    # )
    g001_spr: Path = data_base_path / Path("g001/spr_for_g001_paper2.csv")
    # g001_spr_core_vrc01: Path = data_base_path / Path("g001/Data_S3_20231214.csv")
    g001_spr_core_vrc01: Path = data_base_path / Path("g001/G001-coreg28v2_SPR.csv")
    g001_spr_core_nonvrc01: Path = data_base_path / Path("g001/G001_nonVRC01c_coreg28v2.csv")
    g001_alleles: Path = data_base_path / Path("g001/IGHVI-2_allele_study_output.csv")
    # g001_cluster: Path = data_base_path / Path("g001/cluster_df_w_len.feather")
    g001_cluster: Path = data_base_path / Path("g001/G001_clusters.feather")

    g001_methodology: Path = data_base_path / Path("g001/G001_Re-Analysis_231121_MP_JRW.csv")
    g002_methodology_flow_and_seq: Path = data_base_path / Path("g001/flow_and_seq_g001_mg002.csv")
    g001_10X_g003_flow_and_seq: Path = data_base_path / Path("g001/flow_and_seq_g001_mg002.csv")
    g001_methodology_sequences: Path = data_base_path / Path("g001/g001_sequences_from_10x.feather")
    g002_spr_predictions = data_base_path / "g002/spr/spr_predictions_G002_19April2024.csv"
    g002_spr_tests = data_base_path / "g002/spr/spr_tests_G002_19April2024.csv"
    figure_outdir = Path(__file__).parent.parent.parent.parent / Path("G00X-plots")
    g001_igl_df = data_base_path / Path("g001/igl_remake.feather")

    @validator(
        "*",
        pre=True,
        always=True,
    )
    def validate_path(cls, v: str) -> str:
        if not Path(v).exists():
            raise ValueError(f"{v} does not exist")
        return v


class Data:
    METHODOLOGY_pubIDs = [
        "PubID_023",
        "PubID_005",
        "PubID_051",
        "PubID_187",
        "PubID_153",
        "PubID_079",
        "PubID_047",
        "PubID_154",
        "PubID_110",
        "PubID_046",
    ]

    def __init__(self):
        self.paths = DataPaths()

    # Warning: will be deprecated and will eventually map pubID to pubID
    # with col PTID being named PTID with values of pubID
    def create_PTID2pubID(self) -> str:
        df = pd.read_csv(self.paths.g001_flow_and_seq)
        return {row.PTID: row.PubID for row in df.itertuples()}

    def populate_psname(self, df: pd.DataFrame) -> pd.DataFrame:
        df["pseudogroup"] = df.pseudogroup.astype(int)
        ps2name = {
            k + 1: v.replace(r"$\rightarrow$", "->")
            for k, v in enumerate(
                [
                    "eOD",
                    "core",
                    r"eOD$\rightarrow$eOD",
                    r"eOD$\rightarrow$core",
                    r"eOD$\rightarrow$eOD",
                    r"eOD$\rightarrow$core",
                    r"eOD$\rightarrow$eOD$\rightarrow$core",
                ]
            )
        }
        df.insert(1, "psname", "")
        df["psname"] = df["pseudogroup"].map(ps2name)
        return df

    def populate_psname_spr(self, df: pd.DataFrame) -> pd.DataFrame:
        df["pseudogroup"] = df.pseudogroup.astype(int)
        ps2name = {
            k + 1: v.replace(r"$\rightarrow$", "->")
            for k, v in enumerate(
                [
                    "eOD",
                    r"eOD$\rightarrow$eOD",
                    r"eOD$\rightarrow$core",
                    r"eOD$\rightarrow$eOD",
                    r"eOD$\rightarrow$core",
                    r"eOD$\rightarrow$eOD$\rightarrow$core",
                ]
            )
        }
        df.insert(1, "psname", "")
        df["psname"] = df["pseudogroup"].map(ps2name)
        return df

    def get_g001_seqs(self) -> pd.DataFrame:
        df = pd.read_feather(self.paths.g001_sequences)
        df["pubID"] = df["Subject"].map(self.create_PTID2pubID())
        df["Subject"] = df["pubID"]
        df["PTID"] = df["pubID"]
        df["top_c_call"] = "IGHG"
        df = df
        return df

    def get_g001_sequences_prime(self, allow_week_10: bool = False) -> pd.DataFrame:
        """Get G001 Sequences for prime"""
        df = pd.read_feather(self.paths.g001_sequences)
        df["trial"] = "G001"
        df.rename({"PubID": "pubID", "weeks_post": "weeks"}, axis=1, inplace=True)
        df["weeks"] = df["weeks"].replace({-4: -5})
        if allow_week_10:
            weeks = [-5, 4, 8, 10, 16]
        else:
            weeks = [-5, 4, 8, 16]
        df = df[df["weeks"].isin(weeks)].reset_index(drop=True)
        df = df.query("vaccine_group=='vaccine'").reset_index(drop=True)
        # has a five length
        df["has_5_len"] = False
        df.loc[
            df[df["cdr3_aa_light"].str.len() == 5].index,
            "has_5_len",
        ] = True

        # has_vh12
        df["has_vh1-2"] = False
        df.loc[
            df[df["v_call_top_heavy"].str.split("*").str.get(0) == "IGHV1-2"].index,
            "has_vh1-2",
        ] = True
        df = df[~df["cdr3_aa_light"].isna()].reset_index(drop=True)
        return df

    def get_g001_10x_sequences(self) -> pd.DataFrame:
        df = pd.read_feather(self.paths.g001_methodology_sequences)
        df["pubID"] = df["PTID"].map(self.create_PTID2pubID())
        df["PTID"] = df["pubID"]
        return df

    def get_g001_flow_and_seq_prime(self, allow_week_10: bool = False) -> pd.DataFrame:
        """Get G001 Flow and Seq Summary for prime"""
        flow_and_seq = pd.read_csv(self.paths.g001_flow_and_seq)
        # only did IGG for G001 for this dataset
        flow_and_seq["top_c_call"] = "IGHG"
        visits = ["V02", "V05", "V06", "V07", "V07A", "V08", "V09", "V10"]
        weeks = [-5, 3, 4, 8, 9, 10, 11, 16]
        flow_and_seq["weeks"] = flow_and_seq["Visit"].map(dict(zip(visits, weeks)))
        flow_and_seq.rename({"PubID": "pubID"}, axis=1, inplace=True)
        if allow_week_10:
            weeks = [-5, 4, 8, 10, 16]
        else:
            weeks = [-5, 4, 8, 16]
        flow_and_seq = flow_and_seq[flow_and_seq["weeks"].isin(weeks)].reset_index(drop=True)
        flow_and_seq = flow_and_seq.query("Treatment!='DPBS sucrose'").reset_index(drop=True)
        flow_and_seq["pseudogroup"] = "G001"
        flow_and_seq["trial"] = "G001"

        # Add fields used by G002 and G003 to G001
        g001_mapping = {
            "Percent of GT8++IgG+ B cells that are KO-": "Percent IgG^{+}KO^- among Ag^{++}",
            "Number of epitope-specific (KO-GT8++) sequenced IgG BCRs that are VRC01-class": "Number of IGHG sequences that are VRC01-class",
            "Percent of IgG+ B cells detected as VRC01-class (missing seq to 0)": "Percent of VRC01-class sequences among IgG",
            "Response (missing seq to 0)": "Response x",
            "Percent of GT8++ IgG+ B cells detected as VRC01-class (missing seq to 0)": "Percent VRC01-class among eOD-specific IgG+ memory BCR sequences",
            # "Percent of B cells detected as VRC01-class": "Percent of IGD- sequences that are VRC01-class",
            "Percent of B cells detected as VRC01-class": "Percent of VRC01-class sequences among IgD-",
            # abc IGHG
            "Percent of epitope-specific (KO-GT8++) sequenced IgG BCRs that are VRC01-class": "Percent of IGHG sequences that are VRC01-class",
            "Percent of IgG+ B cells that are epitope-specific (KO-GT8++)": "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}",
            "Percent of IgG+ B cells that are GT8++ (without regard to KO binding status)": "Percent antigen-specific among IgG^{+}",
            # abc IGHD
            # "Percent of epitope-specific (KO-GT8++) sequenced IgD BCRs that are VRC01-class": "Percent of IGHD sequences that are VRC01-class",
        }
        for k, v in g001_mapping.items():
            flow_and_seq[v] = flow_and_seq[k]

        # flow_and_seq["Percent of VRC01-class sequences among IgD-"] = (
        #     flow_and_seq["Percent of VRC01-class sequences among IgD-"] * 100
        # )
        return flow_and_seq

    def get_filtered_g002_flow_and_seq(self) -> pd.DataFrame:
        # df = pd.read_feather(self.paths.g002_flow_and_seq)
        df = pd.read_csv(self.paths.g002_flow_and_seq)

        df["weeks"] = df["weeks"].astype(int)
        pubids_filter = [
            ("G002-479", [16, 24]),
            ("G002-254", [16, 24]),
            ("G002-733", [16, 24]),
            ("G002-632", [24]),
            ("G002-462", [24]),
            ("G002-477", [24]),
            ("G002-689", [16, 24]),
            ("G002-810", [24]),
        ]

        indexes = []
        for pubid, weeks in pubids_filter:
            indexes.extend(df.query(f"weeks in {weeks} and pubID == '{pubid}'").index)

        df = df.drop(indexes)
        return df

    def get_filtered_g002_sequences(self) -> pd.DataFrame:
        df = pd.read_feather(self.paths.g002_sequences)

        df["weeks"] = df["weeks"].astype(int)
        pubids_filter = [
            ("G002-479", [16, 24]),
            ("G002-254", [16, 24]),
            ("G002-733", [16, 24]),
            ("G002-632", [24]),
            ("G002-462", [24]),
            ("G002-477", [24]),
            ("G002-689", [16, 24]),
            ("G002-810", [24]),
        ]

        indexes = []
        for pubid, weeks in pubids_filter:
            indexes.extend(df.query(f"weeks in {weeks} and pubID == '{pubid}'").index)

        df = df.drop(indexes)

        # remove *05/*06 from the v_call_top_heavy
        # df = determine_if_vrc01(df)

        return df

    def get_g002_final_df_raw(self) -> pd.DataFrame:
        return self.get_filtered_g002_sequences()

    def get_g002_sequences(self) -> pd.DataFrame:
        """Get G002 Sequences"""
        df = self.get_filtered_g002_sequences()
        df["trial"] = "G002"
        df = df[~df["cdr3_aa_heavy"].isna()].reset_index(drop=True)
        df = df[~df["cdr3_aa_light"].isna()].reset_index(drop=True)
        df["has_vh1-2"] = False
        df.loc[
            df[df["v_call_top_heavy"].str.split("*").str.get(0) == "IGHV1-2"].index,
            "has_vh1-2",
        ] = True

        df["has_5_len"] = False
        df.loc[
            df[df["cdr3_aa_light"].str.len() == 5].index,
            "has_5_len",
        ] = True
        df["weeks"] = df["weeks"].astype(int)

        hcdr2_range = range(52, 58)

        def find_hcdr2_sets(mutations):
            # if isinstance(mutations, float):
            #     return 0
            hcdr2_mutations = 0
            for x in mutations:
                digit = int(re.findall(r"\d+", x)[0])
                if digit in hcdr2_range:
                    hcdr2_mutations += 1
            return hcdr2_mutations

        df["num_hcdr2_mutations"] = df["cottrell_focused_v_common_heavy_positive"].apply(find_hcdr2_sets)

        return df

    def get_g002_sequences_prime(
        self,
        use_cluster_file: bool = False,
        use_cluster_with_wk24: bool = False,
    ) -> pd.DataFrame:
        """Get G002 Sequences for just the prime"""
        if use_cluster_file:
            if not use_cluster_with_wk24:
                df = pd.read_feather(self.paths.g002_cluster_prime_nowk24)
            else:
                df = pd.read_feather(self.paths.g002_cluster_prime)
        else:
            df = self.get_filtered_g002_sequences()
        # df = pd.read_feather(self.get_filtered_g002_sequences())  # TODO: finish this
        df["trial"] = "G002"
        df = df[~df["cdr3_aa_heavy"].isna()].reset_index(drop=True)
        df = df[~df["cdr3_aa_light"].isna()].reset_index(drop=True)
        df["has_vh1-2"] = False
        df.loc[
            df[df["v_call_top_heavy"].str.split("*").str.get(0) == "IGHV1-2"].index,
            "has_vh1-2",
        ] = True
        df["has_5_len"] = False
        df.loc[
            df[df["cdr3_aa_light"].str.len() == 5].index,
            "has_5_len",
        ] = True
        df["weeks"] = df["weeks"].astype(int)
        df = df.query("group!=4").reset_index(drop=True)
        df_1 = df.query("weeks in [-5,4,8]").reset_index(drop=True)
        df_2 = df.query("weeks==16 and group in [1,3]").reset_index(drop=True)
        df_3 = df.query("weeks==24 and group in [1]").reset_index(drop=True)
        df = pd.concat([df_1, df_2, df_3]).reset_index(drop=True)
        df = df.query("probe_set=='eODGT8'").reset_index(drop=True)
        return df

    @cache
    def get_g00x_sequences_prime(self) -> pd.DataFrame:
        g001_seq = self.get_g001_sequences_prime(allow_week_10=True)
        g001_seq["weeks"] = g001_seq["weeks"].astype(int)
        g001_seq["pseudogroup"] = "G001"
        g001_seq["trial"] = "G001"
        g001_seq["top_c_call"] = "IGHG"
        g001_seq["hcdr3_len"] = g001_seq["cdr3_aa_heavy"].str.len()

        g002_seq = self.get_g002_sequences_prime()
        g002_seq["pseudogroup"] = "G002"
        # Warning: this is a hack to make the weeks match up for plotinging only
        # if sync_weeks:
        #     g002_seq["weeks"] = g002_seq["weeks"].replace({24: 21})

        g003_seq = self.get_g003_sequences_prime()

        g00x_seq = pd.concat(
            [
                g001_seq,
                g002_seq,
                g003_seq,
            ]
        ).reset_index(drop=True)
        return g00x_seq

    def get_g00x_flow_and_seq(self, allow_wk10: bool = True) -> pd.DataFrame:
        g001_flow_and_seq = self.get_g001_flow_and_seq_prime(allow_week_10=allow_wk10)
        g002_flow_and_seq = self.get_g002_flow_and_seq_prime()
        g003_flow_and_seq = self.get_g003_flow_and_seq_prime()

        g00x_flow_and_seq = pd.concat([g001_flow_and_seq, g002_flow_and_seq, g003_flow_and_seq]).reset_index(drop=True)

        g00x_flow_and_seq.loc[g00x_flow_and_seq.query("weeks < 16").index, "shot"] = "first"
        g00x_flow_and_seq.loc[g00x_flow_and_seq.query("weeks > 8").index, "shot"] = "second"

        return g00x_flow_and_seq

    def get_g00x_seq_prime_igg_df(self) -> pd.DataFrame:
        g001_seq_prime_df = self.get_g001_sequences_prime()
        g002_seq_prime_df = self.get_g002_sequences_prime()
        g002_seq_prime_igg_df = g002_seq_prime_df.query("top_c_call=='IGHG'")
        # can only use igg for g001 comparison
        g00x_seq_prime_igg_df = pd.concat([g001_seq_prime_df, g002_seq_prime_igg_df]).reset_index(drop=True)
        g00x_seq_prime_igg_df["v_call_top_light_gene"] = (
            g00x_seq_prime_igg_df["v_call_top_light"].str.split("*").str.get(0)
        )

        g00x_seq_prime_igg_df["v_call_top_heavy_gene"] = (
            g00x_seq_prime_igg_df["v_call_top_heavy"].str.split("*").str.get(0)
        )
        g00x_seq_prime_igg_df.loc[
            g00x_seq_prime_igg_df[g00x_seq_prime_igg_df["v_call_top_heavy_gene"] == "IGHV1-2"].index,
            "has_vh1-2",
        ] = True
        return g00x_seq_prime_igg_df

    def get_g003_merged_prime(self) -> pd.DataFrame:
        df = pd.read_feather(self.paths.g003_merged)
        df["trial"] = "G003"
        df = df.query("probe_set=='eODGT8'").reset_index(drop=True)
        df = df[~df["cdr3_aa_heavy"].isna()].reset_index(drop=True)
        df = df[~df["cdr3_aa_light"].isna()].reset_index(drop=True)
        df["has_vh1-2"] = False
        df.loc[
            df[df["v_call_top_heavy"].str.split("*").str.get(0) == "IGHV1-2"].index,
            "has_vh1-2",
        ] = True

        df["has_5_len"] = False
        df.loc[
            df[df["cdr3_aa_light"].str.len() == 5].index,
            "has_5_len",
        ] = True

        tp2weeks = {
            "V91": -5,
            "V101": 0,
            "V108": 1,
            "V115": 2,
            "V122": 3,
            "V201": 8,
            "V208": 9,
            "V215": 10,
            "V222": 11,
            "V229": 12,
            "V257": 16,
            "V292": 21,
        }

        df["weeks"] = [tp2weeks[tp] for tp in df["timepoint"]]
        df["weeks"] = df["weeks"].astype(int)

        df["group"] = df["pubID"]
        df["group"] = df["group"].replace(
            {pubID: group for pubID, group in self.get_g003_flow_and_seq_prime()[["pubID", "group"]].values}
        )
        # df = df.query("group!=4").reset_index(drop=True)
        # df_1 = df.query("weeks in [-5,4,8]").reset_index(drop=True)
        # df_2 = df.query("weeks==16 and group in [1,3]").reset_index(drop=True)
        # df_3 = df.query("weeks==24 and group in [1]").reset_index(drop=True)
        # df = pd.concat([df_1, df_2, df_3]).reset_index(drop=True)
        # df = df.query("probe_set=='eODGT8'").reset_index(drop=True)
        return df

    def get_g003_sequences_prime(
        self,
        use_cluster_file: bool = False,
        use_cluster_with_wk21: bool = False,
    ) -> pd.DataFrame:
        """Get G002 Sequences for just the prime"""
        if use_cluster_file:
            if not use_cluster_with_wk21:
                df = pd.read_feather(self.paths.g003_cluster_prime_nowk21)
            else:
                df = pd.read_feather(self.paths.g003_cluster_prime)
        else:
            df = pd.read_feather(self.paths.g003_sequences)
        df["trial"] = "G003"
        df["pseudogroup"] = "G003"
        df = df[~df["cdr3_aa_heavy"].isna()].reset_index(drop=True)
        df = df[~df["cdr3_aa_light"].isna()].reset_index(drop=True)
        df["has_vh1-2"] = False
        df.loc[
            df[df["v_call_top_heavy"].str.split("*").str.get(0) == "IGHV1-2"].index,
            "has_vh1-2",
        ] = True

        df["has_5_len"] = False
        df.loc[
            df[df["cdr3_aa_light"].str.len() == 5].index,
            "has_5_len",
        ] = True

        tp2weeks = {
            "V91": -5,
            "V101": 0,
            "V108": 1,
            "V115": 2,
            "V122": 3,
            "V201": 8,
            "V208": 9,
            "V215": 10,
            "V222": 11,
            "V229": 12,
            "V257": 16,
            "V292": 21,
        }

        df["weeks"] = [tp2weeks[tp] for tp in df["timepoint"]]
        df["weeks"] = df["weeks"].astype(int)

        df["group"] = df["pubID"]
        df["group"] = df["group"].replace(
            {pubID: group for pubID, group in self.get_g003_flow_and_seq_prime()[["pubID", "group"]].values}
        )
        # TODO: find out why this is here
        # df = calculate_resonse(df)

        # df = df.query("group!=4").reset_index(drop=True)
        # df_1 = df.query("weeks in [-5,4,8]").reset_index(drop=True)
        # df_2 = df.query("weeks==16 and group in [1,3]").reset_index(drop=True)
        # df_3 = df.query("weeks==24 and group in [1]").reset_index(drop=True)
        # df = pd.concat([df_1, df_2, df_3]).reset_index(drop=True)
        # df = df.query("probe_set=='eODGT8'").reset_index(drop=True)
        return df

    def get_g002_flow_and_seq_prime(self) -> pd.DataFrame:
        """Get G002 Flow and Seq Summary for the prime"""
        flow_and_seq = self.get_filtered_g002_flow_and_seq()
        flow_and_seq = flow_and_seq.query("run_purpose !='PreS'").reset_index(drop=True)
        flow_and_seq["weeks"] = flow_and_seq["weeks"].astype(int)
        groups = [1, 2, 3]
        flow_and_seq_to_8 = flow_and_seq.query("weeks in [-5,4,8] and group in @groups").reset_index(drop=True)
        groups = [1, 3]
        flow_and_seq_to_16 = flow_and_seq.query("weeks==16 and group in @groups").reset_index(drop=True)
        flow_and_seq_to_24 = flow_and_seq.query("weeks==24 and group in [1]").reset_index(drop=True)
        flow_and_seq = pd.concat([flow_and_seq_to_8, flow_and_seq_to_16, flow_and_seq_to_24]).reset_index(drop=True)
        flow_and_seq["pseudogroup"] = "G002"
        flow_and_seq["trial"] = "G002"
        flow_and_seq = flow_and_seq.query("probe_set=='eODGT8'").reset_index(drop=True)

        # for the ones that are redid, take the one with more sequences
        # TODO: tail
        flow_and_seq = (
            flow_and_seq.sort_values("Number of sequences", ascending=False)
            .groupby(["ptid", "weeks", "group", "probe_set"])
            .head(1)
        )

        flow_and_seq = calculate_resonse(flow_and_seq)
        return flow_and_seq

    def get_g003_flow_and_seq_prime(self) -> pd.DataFrame:
        """Get G003 Flow and Seq Summary for the prime"""
        flow_and_seq = pd.read_feather(self.paths.g003_flow_and_seq)
        flow_and_seq = flow_and_seq.query("run_purpose !='PreS'").reset_index(drop=True)
        flow_and_seq["weeks"] = flow_and_seq["weeks"].astype(int)
        flow_and_seq["group"] = flow_and_seq["group"].astype(int)
        groups = [6, 10]
        flow_and_seq_to_8 = flow_and_seq.query("weeks in [-5,8,10] and group in [6, 10]").reset_index(drop=True)
        flow_and_seq_to_16 = flow_and_seq.query("weeks==16 and group in @groups").reset_index(drop=True)
        flow_and_seq_to_21 = flow_and_seq.query("weeks==21 and group in @groups").reset_index(drop=True)
        flow_and_seq = pd.concat([flow_and_seq_to_8, flow_and_seq_to_16, flow_and_seq_to_21]).reset_index(drop=True)
        flow_and_seq["trial"] = "G003"
        flow_and_seq["pseudogroup"] = "G003"
        flow_and_seq = flow_and_seq.query("probe_set=='eODGT8'").reset_index(drop=True)

        # for the ones that are redid, take the one with more sequences
        flow_and_seq = (
            flow_and_seq.sort_values("Number of sequences", ascending=False)
            .groupby(["ptid", "weeks", "group", "probe_set"])
            .head(1)
        )

        flow_and_seq = calculate_resonse(flow_and_seq)
        return flow_and_seq

    def get_g002_flow_and_seq_boost(self) -> pd.DataFrame:
        """Fig 4"""
        flow_and_seq = self.get_filtered_g002_flow_and_seq()
        flow_and_seq["trial"] = "G002"
        flow_and_seq = flow_and_seq.query("run_purpose !='PreS'").reset_index(drop=True)
        flow_and_seq = flow_and_seq.query("probe_set!='eODGT8'").reset_index(drop=True)
        flow_and_seq["weeks"] = flow_and_seq["weeks"].astype(int)

        # Remove reruns filters!
        flow_and_seq = (
            flow_and_seq.sort_values("Number of sequences", ascending=False)
            .groupby(["ptid", "weeks", "group", "probe_set"])
            .head(1)
        )

        # eOD
        flow_and_seq.loc[
            flow_and_seq.query("group in [1,2,3]").query("weeks==8").index,
            "pseudogroup",
        ] = 1
        # core
        flow_and_seq.loc[
            flow_and_seq.query("group ==4").query("weeks==8").index,
            "pseudogroup",
        ] = 2
        # eOD -> eOD
        flow_and_seq.loc[
            flow_and_seq.query("group in [1,3]").query("weeks==16").index,
            "pseudogroup",
        ] = 3
        # eOD -> core
        flow_and_seq.loc[
            flow_and_seq.query("group ==2").query("weeks==16").index,
            "pseudogroup",
        ] = 4
        # eOD -> eOD
        flow_and_seq.loc[
            flow_and_seq.query("group ==1").query("weeks==24").index,
            "pseudogroup",
        ] = 5
        # eOD -> core
        flow_and_seq.loc[
            flow_and_seq.query("group ==2").query("weeks==24").index,
            "pseudogroup",
        ] = 6
        # eOD -> eOD -> core
        flow_and_seq.loc[
            flow_and_seq.query("group ==3").query("weeks==24").index,
            "pseudogroup",
        ] = 7

        flow_and_seq = calculate_resonse_boost(flow_and_seq)

        # flow_and_seq = flow_and_seq[flow_and_seq.pseudogroup.notna()]

        return flow_and_seq

    def get_g002_flow_and_seq_boost_plus(self) -> pd.DataFrame:
        flow_and_seq = self.get_g002_flow_and_seq_boost()
        seq = self.get_g002_sequences_boost()
        flow_and_seq = flow_and_seq.merge(
            seq.groupby(
                [
                    # "run_date",
                    "pubID",
                    "ptid",
                    "group",
                    "weeks",
                    "visit_id",
                    "probe_set",
                    "sample_type",
                ]
            )
            .apply(lambda x: get_better_than(x, 4))  # type:ignore
            .reset_index(),
            on=[
                # "run_date",
                "pubID",
                "ptid",
                "group",
                "weeks",
                "visit_id",
                "probe_set",
                "sample_type",
            ],
            how="left",
        )
        flow_and_seq = flow_and_seq.merge(
            seq.groupby(
                [
                    # "run_date",
                    "pubID",
                    "ptid",
                    "group",
                    "weeks",
                    "visit_id",
                    "probe_set",
                    "sample_type",
                ]
            )
            .apply(lambda x: get_better_than_hcdr2(x, 1))  # type:ignore
            .reset_index(),
            on=[
                # "run_date",
                "pubID",
                "ptid",
                "group",
                "weeks",
                "visit_id",
                "probe_set",
                "sample_type",
            ],
            how="left",
        )

        flow_and_seq["Percent VRC01 gt n among IgG"] = (
            flow_and_seq["Percent of VRC01-class sequences among IgG"] / 100 * flow_and_seq["fraction_of_cottrell_gt_n"]
        ) * 100

        flow_and_seq["Percent VRC01 gt n among IgG hcdr2"] = (
            flow_and_seq["Percent of VRC01-class sequences among IgG"]
            / 100
            * flow_and_seq["fraction_of_cottrell_gt_hcdr2_n"]
        ) * 100

        # lt marker is our threshold
        # lt = 1e-5
        # flow_and_seq.loc[
        #     flow_and_seq[flow_and_seq["Percent VRC01 gt n among IgG"] < lt].index, "Percent VRC01 gt n among IgG"
        # ] = 1e-5
        # flow_and_seq.loc[
        #     flow_and_seq[flow_and_seq["Percent VRC01 gt n among IgG hcdr2"] < lt].index,
        #     "Percent VRC01 gt n among IgG hcdr2",
        # ] = 1e-5
        return flow_and_seq
        # account for duplicate week 4 & 8; group 4 runs
        # flow_and_seq = flow_and_seq.drop(
        #     index=flow_and_seq.query("weeks in [4,8]")
        #     .query("group==4")
        #     .query('pubID in ["G002-611", "G002-710"]')
        #     .sort_values("run_date")
        #     .drop_duplicates(["ptid", "weeks", "group"])
        #     .index
        # )
        # sequences_boost = self.get_g002_sequences_boost()
        # pseudogroup_cols = [
        #     "run_date",
        #     "pubID",
        #     "ptid",
        #     "group",
        #     "weeks",
        #     "visit_id",
        #     "probe_set",
        #     "sample_type",
        # ]
        # flow_and_seq = flow_and_seq.merge(
        #     sequences_boost.query("top_c_call=='IGHG'")
        #     .groupby(pseudogroup_cols)
        #     .apply(
        #         lambda x: get_fraction_gt(
        #             x,
        #             "cottrell_focused_v_common_score",
        #             4,
        #             "fraction_of_cottrell_gt_4",
        #         )
        #     )
        #     .reset_index(),
        #     on=pseudogroup_cols,
        #     how="left",
        # )
        # flow_and_seq = flow_and_seq.merge(
        #     sequences_boost.query("top_c_call=='IGHG'")
        #     .groupby(pseudogroup_cols)
        #     .apply(
        #         lambda x: get_fraction_gt(
        #             x,
        #             "num_hcdr2_mutations",
        #             1,
        #             "fraction_of_cottrell_gt_hcdr2_1",
        #         )
        #     )
        #     .reset_index(),
        #     on=pseudogroup_cols,
        #     how="left",
        # )
        # flow_and_seq = flow_and_seq.merge(
        #     flow_and_seq.groupby("pseudogroup")
        #     .apply(lambda x: get_fraction_eq(x, "Response y", 1.0, "fraction_of_responders"))
        #     .reset_index(),
        #     on=["pseudogroup"],
        #     how="left",
        # )

    def get_g002_sequences_boost(self, use_cluster_file: bool = False, use_filtered: bool = True) -> pd.DataFrame:
        if use_cluster_file:
            sequences = pd.read_feather(self.paths.g002_cluster)
        elif use_filtered:
            sequences = self.get_filtered_g002_sequences()
        else:
            print("Using unfiltered sequences")
            sequences = pd.read_feather(self.paths.g002_sequences)

        sequences["weeks"] = sequences["weeks"].astype(int)
        sequences["trial"] = "G002"
        sequences = sequences.query("probe_set!='eODGT8'").reset_index(drop=True)
        # eOD
        sequences.loc[
            sequences.query("group in [1,2,3]").query("weeks==8").index,
            "pseudogroup",
        ] = 1
        # core
        sequences.loc[sequences.query("group ==4").query("weeks==8").index, "pseudogroup"] = 2
        # eOD -> eOD
        sequences.loc[
            sequences.query("group in [1,3]").query("weeks==16").index,
            "pseudogroup",
        ] = 3
        # eOD -> core
        sequences.loc[sequences.query("group ==2").query("weeks==16").index, "pseudogroup"] = 4
        # eOD -> eOD
        sequences.loc[sequences.query("group ==1").query("weeks==24").index, "pseudogroup"] = 5
        # eOD -> core
        sequences.loc[sequences.query("group ==2").query("weeks==24").index, "pseudogroup"] = 6
        # eOD -> eOD -> core
        sequences.loc[sequences.query("group ==3").query("weeks==24").index, "pseudogroup"] = 7

        hcdr2_range = range(52, 58)

        def find_hcdr2_sets(mutations):
            # if isinstance(mutations, float):
            #     return 0
            hcdr2_mutations = 0
            for x in mutations:
                digit = int(re.findall(r"\d+", x)[0])
                if digit in hcdr2_range:
                    hcdr2_mutations += 1
            return hcdr2_mutations

        sequences["num_hcdr2_mutations"] = sequences["cottrell_focused_v_common_heavy_positive"].apply(find_hcdr2_sets)
        return sequences

    def get_g002_sequences_boost_plus(self) -> pd.DataFrame:
        g002_sequences = self.get_g002_sequences_boost()
        g002_sequences = g002_sequences.merge(
            g002_sequences.groupby(
                [
                    # "run_date",
                    "pubID",
                    "ptid",
                    "group",
                    "weeks",
                    "visit_id",
                    "probe_set",
                    "sample_type",
                ]
            )
            .apply(lambda x: get_better_than(x, 4))  # type:ignore
            .reset_index(),
            on=[
                # "run_date",
                "pubID",
                "ptid",
                "group",
                "weeks",
                "visit_id",
                "probe_set",
                "sample_type",
            ],
            how="left",
        )
        g002_sequences = g002_sequences.merge(
            g002_sequences.groupby(
                [
                    # "run_date",
                    "pubID",
                    "ptid",
                    "group",
                    "weeks",
                    "visit_id",
                    "probe_set",
                    "sample_type",
                ]
            )
            .apply(lambda x: get_better_than_hcdr2(x, 1))  # type:ignore
            .reset_index(),
            on=[
                # "run_date",
                "pubID",
                "ptid",
                "group",
                "weeks",
                "visit_id",
                "probe_set",
                "sample_type",
            ],
            how="left",
        )

        return g002_sequences

    def get_trial_palette(self) -> dict[str, str]:
        """Get color pallette for G001 and G002 Trials"""
        return {"G001": "#6C65FF", "G002": "#91FCC0"}

    def get_trial_g001_g002_g003_palette(self) -> dict[str, str]:
        """Get color pallette for G001, G002, and G003 Trials"""
        return {
            "G001": "#6C65FF",
            "G002": "#91FCC0",
            "G003": "#E377C2",
        }

    def get_g003_site_palette(self) -> dict[str, str]:
        return {
            "G003": "#bf9b30",
            "CFHR": "#17BFD0",
            "Aurum": "#91FCC0",
            "NAC": "#FF6F61",
        }

    def get_g003_trial_palette(self) -> dict[str, str]:
        """Get color pallette for G002 and G003 Trials"""
        return {"G002": "#91FCC0", "G003": "#6C65FF"}

    def get_pseudogroup_palette(self) -> dict[int, str]:
        pseudogroup_colors = {
            1: "#91FCC0",
            2: "#E377C2",
            3: "#2078B4",
            5: "#9567BD",
            6: "#BDBD23",
            7: "#FF7F0F",
            8: "white",
        }
        return pseudogroup_colors

    def get_week_palette(self) -> dict[int, str]:
        week = {
            -5: "#9567BD",
            4: "#17BFD0",
            8: "gold",
            10: "gold",
            16: "#E377C2",
            21: "#2078B4",
            24: "#2078B4",
        }
        return week

    def get_dekosky_vh12(self) -> pd.DataFrame:
        """Get Dekosky VH12"""
        df = pd.read_feather(self.paths.dekosky_vh12)
        df["v_call_top_light_gene"] = df["v_call_top_light"]
        df["v_call_top_heavy_gene"] = df["v_call_top_heavy"]
        return df

    def get_oas_5_len(self) -> pd.DataFrame:
        """Get OAS 5 Length"""
        df = pd.read_feather(self.paths.oas_5_len)
        df = df.query("oas_subject != 'None'").reset_index(drop=True)
        return df

    def get_vrc01_class_bnabs(self) -> pd.DataFrame:
        """Get VRC01 Class Bnabs"""
        df = pd.read_feather(self.paths.vrc01_class_bnabs_path)
        hcdr2_range = range(52, 58)

        def find_hcdr2_sets(mutations):
            # if isinstance(mutations, float):
            #     return 0
            hcdr2_mutations = 0
            for x in mutations:
                digit = int(re.findall(r"\d+", x)[0])
                if digit in hcdr2_range:
                    hcdr2_mutations += 1
            return hcdr2_mutations

        df["num_hcdr2_mutations"] = df["cottrell_focused_v_common_heavy_positive"].apply(find_hcdr2_sets)

        return df

    def get_oas_vh12(self) -> pd.DataFrame:
        """Get OAS VH12"""
        df = pd.read_feather(self.paths.oas_vh12)
        return df

    def get_human_naive_5_len_df(self) -> pd.DataFrame:
        """Get Human Naive 5 Length"""
        df = pd.read_feather(self.paths.human_naive_path)
        return df

    def get_kappa_vrc01_select_df(self) -> pd.DataFrame:
        """Get Kappa VRC01 Select"""
        df = pd.read_feather(self.paths.select_kappa_vh12_path)
        return df

    def get_lambda_vrc01_select_df(self) -> pd.DataFrame:
        """Get Lambda VRC01 Select"""
        df = pd.read_feather(self.paths.select_lambda_vh12_path)
        return df

    def get_clustered_g001_seqs(self) -> pd.DataFrame:
        """Get Clustered G001 Sequences"""
        df = pd.read_feather(self.paths.g001_clustered_sequences)
        df["ptid"] = df["PubID"]
        df["pubID"] = df["PubID"]
        df["pubid"] = df["PubID"]
        df["timepoint"] = df["Timepoint"]
        df["dose_group"] = df["Dose_Group"]
        return df

    def get_response_count_dataframe(self, include_g003: bool = False, allow_g001_wk_10: bool = False) -> pd.DataFrame:
        g001_flow = self.get_g001_flow_and_seq_prime(allow_week_10=allow_g001_wk_10)
        g002_flow = self.get_g002_flow_and_seq_prime()
        g003_flow = self.get_g003_flow_and_seq_prime()

        if not include_g003:
            combined = pd.concat([g001_flow, g002_flow]).reset_index(drop=True)
        else:
            combined = pd.concat([g001_flow, g002_flow, g003_flow]).reset_index(drop=True)

        # get rid of pre vaccine weeks
        combined = combined.drop(combined.query("weeks < 0").index)
        # combined = combined[combined.weeks > 0]
        print(combined["weeks"].value_counts(dropna=False))

        combined.loc[combined.query("weeks < 16").index, "shot"] = "first"
        combined.loc[combined.query("weeks > 8").index, "shot"] = "second"

        # create a new dataframe for the response
        new_df = []

        # do it once grouping by shot
        for (shot, trial), group_df in combined.groupby(["shot", "trial"]):
            responders = len(group_df.query("`Response x`==1")["pubID"].unique())
            total = len(group_df["pubID"].unique())
            # TODO: manually fixed this to reflect participant change; should be more robust
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
            responders = len(group_df.query("`Response x`==1")["pubID"].unique())
            # TODO: manually fixed this to reflect participant change; should be more robust
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

    def get_g003_response_count_dataframe(self) -> pd.DataFrame:
        g003_flow = self.get_g003_flow_and_seq_prime()
        g002_flow = self.get_g002_flow_and_seq_prime()
        g003_flow = g003_flow.rename({"Response (missing seq to 0)": "Response x"}, axis=1)
        combined = pd.concat([g003_flow, g002_flow]).reset_index(drop=True)

        # get rid of pre vaccine weeks
        combined = combined.drop(combined.query("weeks < 0").index)

        combined.loc[combined.query("weeks < 16").index, "shot"] = "first"
        combined.loc[combined.query("weeks > 8").index, "shot"] = "second"

        # create a new dataframe for the response
        new_df = []

        # do it once grouping by shot
        for (shot, trial), group_df in combined.groupby(["shot", "trial"]):
            responders = len(group_df.query("`Response x`==1")["pubID"].unique())
            total = len(group_df["pubID"].unique())
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
            responders = len(group_df.query("`Response x`==1")["pubID"].unique())
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

    def get_personalized_alleles(self) -> pd.DataFrame:
        """Get Personalized Alleles"""
        df = pd.read_csv(self.paths.personalize_path, index_col=0)
        return df

    def get_g003_personalized_alleles(self) -> pd.DataFrame:
        """Get Personalized Alleles"""
        df = pd.read_csv(self.paths.g003_personalize_path, index_col=0)
        return df

    def get_g003_spr_df(self) -> pd.DataFrame:
        df = pd.read_feather(self.paths.g003_spr)
        # df = df[~df.weeks.isna()]
        return df

    def get_g002_spr_df(self) -> pd.DataFrame:
        df = pd.read_feather(self.paths.g002_spr)
        # df = df[~df.weeks.isna()]
        return df

    def get_g002_spr_predictions(self) -> pd.DataFrame:
        df = pd.read_csv(self.paths.g002_spr_predictions)
        return df

    def get_g002_spr_tests(self) -> pd.DataFrame:
        df = pd.read_csv(self.paths.g002_spr_tests)
        return df

    def get_g002_spr_df_eod_to_core(self) -> pd.DataFrame:
        g002_df = pd.read_feather(self.paths.g002_spr)
        g002_df["is_vrc01_class"] = g002_df.is_vrc01_class.astype(bool)
        g002_df["num_hcdr2_mutations"] = g002_df["num_hcdr2_mutations"].fillna(0).astype(int)
        g002_df["v_call_top_light_gene"] = g002_df["v_call_top_light"].str.split("*").str.get(0)
        g002_df["v_call_top_heavy_gene"] = g002_df["v_call_top_heavy"].str.split("*").str.get(0)
        g002_df["has_5_len"] = False
        g002_df.loc[
            g002_df[g002_df["cdr3_aa_light"].str.len() == 5].index,
            "has_5_len",
        ] = True
        # g002_df = g002_df.sort_values("Chi2").groupby(["Ligand", "Analyte"]).head(1)

        return g002_df

    def get_g002_spr_df_boost(self) -> pd.DataFrame:
        """Get G002 SPR Dataframe"""
        g002_df = pd.read_feather(self.paths.g002_spr)
        g002_df = g002_df[g002_df.Ligand.str.startswith("G002")]
        # g002_df = g002_df.sort_values("Chi2").groupby(["Ligand", "Analyte"]).head(1)
        g002_df["is_vrc01_class"] = g002_df.is_vrc01_class.astype(bool)
        g002_df["num_hcdr2_mutations"] = g002_df["num_hcdr2_mutations"].fillna(0).astype(int)

        g002_df["v_call_top_light_gene"] = g002_df["v_call_top_light"].str.split("*").str.get(0)
        g002_df["v_call_top_heavy_gene"] = g002_df["v_call_top_heavy"].str.split("*").str.get(0)
        g002_df.loc[
            g002_df.query("group in [1,2,3]").query("weeks==8").query("probe_set!='eODGT8'").index,
            "pseudogroup",
        ] = 1
        g002_df.loc[
            g002_df.query("group in [1,3]").query("weeks==16").query("probe_set!='eODGT8'").index,
            "pseudogroup",
        ] = 2
        g002_df.loc[
            g002_df.query("group ==2").query("weeks==16").query("probe_set!='eODGT8'").index,
            "pseudogroup",
        ] = 3
        g002_df.loc[
            g002_df.query("group ==1").query("weeks==24").query("probe_set!='eODGT8'").index,
            "pseudogroup",
        ] = 4
        g002_df.loc[
            g002_df.query("group ==2").query("weeks==24").query("probe_set!='eODGT8'").index,
            "pseudogroup",
        ] = 5
        g002_df.loc[
            g002_df.query("group ==3").query("weeks==24").query("probe_set!='eODGT8'").index,
            "pseudogroup",
        ] = 6
        g002_df = g002_df.query("probe_set!='eODGT8'").reset_index(drop=True)
        g002_df["has_5_len"] = False
        g002_df.loc[
            g002_df[g002_df["cdr3_aa_light"].str.len() == 5].index,
            "has_5_len",
        ] = True

        # g002_df = g002_df.sort_values("Chi2").groupby(["Ligand", "Analyte"]).head(1)
        g002_df["mutations_heavy_count"] = (
            g002_df["mutations_heavy"]
            .str.replace("\n", "")
            .str.replace("' '", "', '")
            .fillna("[]")
            .apply(literal_eval)
            .str.len()
        )
        g002_df["mutations_light_count"] = (
            g002_df["mutations_light"]
            .str.replace("\n", "")
            .str.replace("' '", "', '")
            .fillna("[]")
            .apply(literal_eval)
            .str.len()
        )
        g002_df["mutations_heavy_and_light_count"] = g002_df["mutations_heavy_count"] + g002_df["mutations_light_count"]

        hcdr2_range = range(52, 58)

        def find_hcdr2_sets(mutations):
            # if isinstance(mutations, float):
            #     return 0
            hcdr2_mutations = 0
            for x in mutations:
                digit = int(re.findall(r"\d+", x)[0])
                if digit in hcdr2_range:
                    hcdr2_mutations += 1
            return hcdr2_mutations

        g002_df["cottrell_focused_v_common_heavy_positive"] = (
            g002_df["cottrell_focused_v_common_heavy_positive"]
            .str.replace("\n", "")
            .str.replace("' '", "', '")
            .fillna("[]")
            .apply(literal_eval)
            # .str.len()
        )

        g002_df["num_hcdr2_mutations"] = g002_df["cottrell_focused_v_common_heavy_positive"].apply(find_hcdr2_sets)
        return g002_df

    def get_g003_spr_df_prime(self) -> pd.DataFrame:
        """Get G002 SPR Dataframe"""
        g003_timepoint_to_week_mapping = {
            "V91": -5,
            "V101": 0,
            "V108": 1,
            "V115": 2,
            "V122": 3,
            "V201": 8,
            "V208": 9,
            "V215": 10,
            "V222": 11,
            "V229": 12,
            "V257": 16,
            "V292": 21,
        }
        df = pd.read_feather(self.paths.g003_spr)
        df = df[df.Ligand.str.startswith("G003")]
        df["weeks"] = df.Ligand.str.split("_").str[1].map(g003_timepoint_to_week_mapping)
        df["is_igl"] = df["Ligand"].apply(lambda x: True if "_iGL_" in x else False)
        df["is_vrc01_class"] = df["Ligand"].apply(lambda x: True if "_VRC01" in x else False)
        df["trial"] = "G003"
        df["pseudogroup"] = "G003"
        df["is_cp"] = False
        df.loc[df["is_igl"], "KD_fix_iGL"] = df["KD_fix"]
        # df["KD_fix_iGL"] = df["KD_fix"]
        # df = df[~df.weeks.isna()]
        df["weeks"] = df["weeks"].astype(int)
        # df = df.query("group!=4").reset_index(drop=True)
        # df_1 = df.query("weeks in [-5,4,8]").reset_index(drop=True)
        # df_2 = df.query("weeks==16 and group in [1,3]").reset_index(drop=True)
        # df_3 = df.query("weeks==21 and group in [1]").reset_index(drop=True)
        # df = pd.concat([df_1, df_2, df_3]).reset_index(drop=True)
        # df = df.query("probe_set=='eODGT8'").reset_index(drop=True)

        # df = df.sort_values("Chi2").groupby(["Ligand", "Analyte"]).head(1)
        return df

    def get_g002_spr_df_prime(self) -> pd.DataFrame:
        """Get G002 SPR Dataframe"""
        df = pd.read_feather(self.paths.g002_spr)
        df = df[~df.weeks.isna()]
        df["weeks"] = df["weeks"].astype(int)
        df = df.query("group!=4").reset_index(drop=True)
        df_1 = df.query("weeks in [-5,4,8]").reset_index(drop=True)
        df_2 = df.query("weeks==16 and group in [1,3]").reset_index(drop=True)
        df_3 = df.query("weeks==24 and group in [1]").reset_index(drop=True)
        df = pd.concat([df_1, df_2, df_3]).reset_index(drop=True)
        df = df.query("probe_set=='eODGT8'").reset_index(drop=True)
        df["trial"] = "G002"
        df["pseudogroup"] = "G002"
        # df = df.sort_values("Chi2").groupby(["Ligand", "Analyte"]).head(1)
        return df

    def get_g001_spr_df(self) -> pd.DataFrame:
        timepoint2weeks = {
            "CLK": 0,
            "V05": 3,
            "V07": 8,
            "V07A": 9,
            "V09": 11,
            "V02": -4,
            "V06": 4,
            "V10": 16,
            "V08": 10,
        }
        """Get G001 SPR Dataframe"""
        # df = pd.read_feather(self.paths.g001_spr)
        df = pd.read_csv(self.paths.g001_spr)
        df = df[~df["cellid"].isna()]
        df = df[df.cellid.str.startswith("G001")]

        Ligand2cellid = {r["Ligand"]: r["cellid"] for i, r in df.iterrows()}

        def get_core(path: Path) -> pd.DataFrame:
            core = pd.read_csv(path).query('Analyte=="core-g28v2"').reset_index(drop=True)  # .iloc[:, :30]
            core = core[core.Ligand.str.startswith("G001")]
            core["timepoint"] = core["Ligand"].str.split("_").str[2]
            core["weeks"] = core["timepoint"].map(timepoint2weeks)
            core["weeks_post"] = core["weeks"]
            core["Analyte"] = "core-Hx_r4.0D_TH6_g28v2_pHLsecAvi"
            core["KD_fix"] = core["KD_fix"].apply(
                lambda x: (float(x.strip().lower().replace("na", "nan")) if type(x) == str else float(x))
            )
            core["Chi2"] = core["Chi2"].apply(
                lambda x: (float(x.strip().lower().replace("na", "nan")) if type(x) == str else float(x))
            )
            # core["KD_fix"] = core["KD_fix"].astype(float)
            # core["is_vrc01_class"] = True
            core["is_igl"] = False
            # for i, row in core.iterrows():
            #     if row["Type"].strip() == "VRC01c":
            #         core.loc[i, "is_vrc01_class"] = True
            #     elif row["Type"].strip() == "VRC01c-iGL":
            #         core.loc[i, "is_vrc01_class"] = True
            #         core.loc[i, "is_igl"] = True
            # core["Chi2"] = core["Chi2"].astype(float)
            core["cellid"] = core["Ligand"].map(Ligand2cellid)
            core = core.reset_index(drop=True)
            return core

        def repair_spr(df: pd.DataFrame) -> pd.DataFrame:
            df["cellid"] = df["Ligand"].map(Ligand2cellid)
            df["weeks"] = df["weeks_post"].astype(int)
            return df

        core_vrc01 = get_core(self.paths.g001_spr_core_vrc01)
        core_vrc01["is_vrc01_class"] = True
        core_nonvrc01 = get_core(self.paths.g001_spr_core_nonvrc01)
        core_nonvrc01["KD_fix"] = core_nonvrc01["KD_fix"].fillna(1e-4)
        core_nonvrc01["is_vrc01_class"] = False

        df = pd.concat([df, core_vrc01, core_nonvrc01]).reset_index(drop=True)

        key = pd.read_csv(self.paths.data_base_path / "g001/G001_SPR_ligands-key-mut.csv")
        key = repair_spr(key)
        # breakpoint()
        df = pd.merge(df, key[["Ligand", "key_mutations"]], on="Ligand", how="left")
        df["cottrell_focused_v_common_score"] = df["key_mutations"]

        df["KD_fix"] = df["KD_fix"].apply(
            lambda x: (float(x.strip().lower().replace("na", "nan")) if type(x) == str else float(x))
        )
        df["Chi2"] = df["Chi2"].apply(
            lambda x: (float(x.strip().lower().replace("na", "nan")) if type(x) == str else float(x))
        )
        df["KD_fix_iGL"] = np.nan
        df["KD_fix_iGL"] = df["KD_fix"].where(df["is_igl"] == True, df["KD_fix_iGL"])
        df["hcdr3_len"] = df["cdr3_aa_heavy"].str.len()
        df["is_cp"] = False
        df["top_c_call"] = "IGHG"
        df["pubID"] = df.cellid.str.split("_").str[0]  # .str.join("_")
        df["pubID"] = df["pubID"].map(self.create_PTID2pubID())
        df["trial"] = "G001"
        df["weeks"] = df["weeks_post"].astype(int)
        df["pseudogroup"] = "G001"

        # remove missing nonvrc01 (no cellid mappings) rows that should not be counted
        df = df[~df.cellid.isna()]

        # df = df.query("weeks in [4,8,10,16]")
        return df

    def get_g001_cluster_ses(self) -> pd.DataFrame:
        """Get G001 Clustered Sequences"""
        df = pd.read_feather(self.paths.g001_cluster)
        df["PTID"] = df["pubid"]
        df["pubID"] = df["pubid"]
        return df

    def get_pre_post_core_df(self) -> pd.DataFrame:
        """supp fig 24"""
        df = self.get_filtered_g002_flow_and_seq()
        df = df.query("run_purpose=='Sort'").query("probe_set=='Cg28v2'").reset_index(drop=True)
        df["weeks"] = df["weeks"].astype(int)

        # Remove reruns filter! : It was missing this filter!
        df = (
            df.sort_values("Number of sequences", ascending=False)
            .groupby(["ptid", "weeks", "group", "probe_set"])
            .head(1)
        )

        # Gp 4 week -5
        df.loc[df.query("group==4").query("weeks==-5").index, "pseudogroup"] = 1
        # Gp 4 week 4
        df.loc[df.query("group==4").query("weeks==4").index, "pseudogroup"] = 2
        # Gp 4 week 8
        df.loc[df.query("group==4").query("weeks==8").index, "pseudogroup"] = 3
        # Gp 2 week 8
        df.loc[df.query("group==2").query("weeks==8").index, "pseudogroup"] = 4

        df = calculate_resonse_boost(df)

        return df

    def get_g002_methodology_flow_and_seq(self) -> pd.DataFrame:
        df = pd.read_csv(self.paths.g002_methodology_flow_and_seq)
        # g001_mapping = {
        #     # "Percent of IgG+ B cells that are GT8++ (without regard to KO binding status)": "Percent antigen-specific among IgG^{+}",
        #     # "Percent of IgG+ B cells that are epitope-specific (KO-GT8++)": "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}",
        #     "Percent of GT8++IgG+ B cells that are KO-": "Percent IgG^{+}KO^- among Ag^{++}",
        #     # "Number of epitope-specific (KO-GT8++) sequenced IgG BCRs that are VRC01-class": "Number of IGHG sequences that are VRC01-class",
        #     # "Percent of IgG+ B cells detected as VRC01-class (missing seq to 0)": "Percent of VRC01-class sequences among IgG",
        #     # "Response (missing seq to 0)": "Response x",
        #     # "Percent of epitope-specific (KO-GT8++) sequenced IgG BCRs that are VRC01-class": "Percent of IGHG sequences that are VRC01-class",
        #     # "Percent of GT8++ IgG+ B cells detected as VRC01-class (missing seq to 0)": "Percent VRC01-class among eOD-specific IgG+ memory BCR sequences",
        # }
        # for k, v in g001_mapping.items():
        #     df[v] = df[k]
        df["pubID"] = df["PTID"].map(self.create_PTID2pubID())
        df["PTID"] = df["pubID"]
        return df

    def get_g003_methodology_flow_and_seq(self) -> pd.DataFrame:
        df = pd.read_csv(self.paths.g003_methodology_flow_and_seq)
        df["pubID"] = df["pubID"].map(self.create_PTID2pubID())
        df["PTID"] = df["pubID"]
        column_g003_2_g001_mapping = {
            "Percent antigen-specific among IgG^{+}": "percentage_of_antigen++",
            # "Percent epitope-specific (KO^-Ag^{++}) among IgD^-": "percentage_of_epitope_specific",
            "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}": "percentage_of_epitope_specific",
            "pubID": "PTID",
            "Percent of VRC01-class sequences among IgG": "Percent of IgG+ B cells detected as VRC01-class (missing seq to 0)",
            # "v_mutation_aa_heavy": "v_mutation_aa_heavy",
        }
        for column_g003, column_g001 in column_g003_2_g001_mapping.items():
            df[column_g001] = df[column_g003]

        df["method"] = df["run_purpose"].apply(lambda x: "KWTRP" if x == "G001Sort" else "NAC")
        df = df.query('method!="NAC"')
        return df

    def get_g003_methodology_seq(self) -> pd.DataFrame:
        df = pd.read_feather(self.paths.g003_methodology_seq)
        df["pubID"] = df["ptid"].map(self.create_PTID2pubID())
        df.pop("ptid")
        df["PTID"] = df["pubID"]
        return df
