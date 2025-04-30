import pandas as pd

from .data import Data


class Transforms:
    data = Data()
    colors = ["#91FCC0", "#2078B4", "#E377C2", "#9567BD", "#17BFD0", "#FF7F0F", "#BDBD23", "white"]
    minimum_set = [
        "12A12",
        "12A21",
        "N6",
        "VRC27",
        "N17",
        "N60P1.1",
        "N60P25.1",
        "N60P23",
        "PCIN63_71I",
        "PCIN63_66B",
        "PCIN63_71G",
        "NIH45-46",
        "VRC07b",
        "VRC23",
        "VRC01",
        "VRC02",
        "VRC18",
        "VRC08",
        "VRC-PG19",
    ]

    def __init__(self) -> None:
        self.g00x_seq_prime_igg_df = self.data.get_g00x_seq_prime_igg_df()
        self.dekosky_vh12_df = self.data.get_dekosky_vh12()

    def get_g00x_seq_prime_igg_vrc01_class_pivot_trial_light_value_counts_df(
        self, light_genes: list[str] | None = None, normalize: bool = True
    ) -> pd.DataFrame:
        light_genes = [
            "IGKV1-33",
            "IGKV3-20",
            "IGKV1-5",
            "IGKV3-15",
            "IGLV2-14",
            "IGLV2-23",
            "IGLV2-11",
            "Other",  # anything not in this list is dumped into "Other"
        ]
        g00x_seq_prime_igg_df = self.g00x_seq_prime_igg_df.copy()
        g00x_seq_prime_igg_df.loc[
            g00x_seq_prime_igg_df[g00x_seq_prime_igg_df["v_call_top_light_gene"] == "IGKV3D-15"].index,
            "v_call_top_light_gene",
        ] = "IGKV3-15"  # TODO: might want to ingest this to data
        g00x_seq_prime_igg_df.loc[
            g00x_seq_prime_igg_df[~g00x_seq_prime_igg_df["v_call_top_light_gene"].isin(light_genes)].index,
            "v_call_top_light_gene",
        ] = "Other"
        g00x_seq_prime_igg_vrc01_class_df = g00x_seq_prime_igg_df.query("is_vrc01_class==True")
        g00x_seq_prime_igg_vrc01_class_pivot_trial_df = (
            g00x_seq_prime_igg_vrc01_class_df.groupby("trial")
            .apply(lambda x: x["v_call_top_light_gene"].value_counts(normalize=normalize))
            .reset_index()
            .set_index("trial")
            .loc[:, light_genes]
        )
        return g00x_seq_prime_igg_vrc01_class_pivot_trial_df

    def g00x_seq_prime_igg_non_vrc01_class_has_vh1_2_pivot_trial_light_value_counts_df(
        self, light_genes: list[str] | None = None, normalize: bool = True
    ) -> pd.DataFrame:
        """_summary_

        Parameters
        ----------
        light_genes : list[str] | None, optional
            _description_, by default None
        normalize : bool, optional
            _description_, by default True

        Returns
        -------
        pd.DataFrame
            _description_
        """
        light_genes = [
            "IGKV1-33",
            "IGKV3-20",
            "IGKV1-5",
            "IGKV3-15",
            "IGLV2-14",
            "IGLV2-23",
            "IGLV2-11",
            "Other",  # anything not in this list is dumped into "Other"
        ]
        g00x_seq_prime_igg_df = self.g00x_seq_prime_igg_df.copy()
        g00x_seq_prime_igg_df.loc[
            g00x_seq_prime_igg_df[g00x_seq_prime_igg_df["v_call_top_light_gene"] == "IGKV3D-15"].index,
            "v_call_top_light_gene",
        ] = "IGKV3-15"  # TODO: might want to ingest this to data
        g00x_seq_prime_igg_df.loc[
            g00x_seq_prime_igg_df[~g00x_seq_prime_igg_df["v_call_top_light_gene"].isin(light_genes)].index,
            "v_call_top_light_gene",
        ] = "Other"
        g00x_seq_prime_igg_non_vrc01_class_has_vh1_2_df = g00x_seq_prime_igg_df.query("is_vrc01_class==False").query(
            "`has_vh1-2`==True"
        )
        g00x_seq_prime_igg_non_vrc01_class_has_vh1_2_pivot_trial_light_value_counts_df = (
            g00x_seq_prime_igg_non_vrc01_class_has_vh1_2_df.groupby("trial")
            .apply(lambda x: x["v_call_top_light_gene"].value_counts(normalize=normalize))
            .reset_index()
            .pivot(index="trial", columns="level_1", values="v_call_top_light_gene")
            .loc[:, light_genes]
        )
        return g00x_seq_prime_igg_non_vrc01_class_has_vh1_2_pivot_trial_light_value_counts_df

    def get_dekosky_pivot_replicate_donors_by_light_value_counts_df(
        self, light_genes: list[str] | None = None, normalize: bool = True
    ) -> pd.DataFrame:
        light_genes = [
            "IGKV1-33",
            "IGKV3-20",
            "IGKV1-5",
            "IGKV3-15",
            "IGLV2-14",
            "IGLV2-23",
            "IGLV2-11",
            "Other",  # anything not in this list is dumped into "Other"
        ]
        dekosky_df = self.dekosky_vh12_df
        dekosky_df.loc[
            dekosky_df[dekosky_df["v_call_top_light_gene"] == "IGKV3D-15"].index,
            "v_call_top_light_gene",
        ] = "IGKV3-15"

        # anything not in explicit mapping is other
        dekosky_df.loc[
            dekosky_df[~dekosky_df["v_call_top_light_gene"].isin(light_genes)].index,
            "v_call_top_light_gene",
        ] = "Other"
        dekosky_pivot_replicate_donors_by_light_value_counts_df = (
            dekosky_df.groupby(["Replicate", "donor"])
            .apply(lambda x: x["v_call_top_light_gene"].value_counts(normalize=normalize))
            .reset_index()
            .groupby(["level_2"])
            .mean()
            .drop("Replicate", axis=1)
            .reset_index()
        )
        # Throw group into x-axis for easy pivoting
        dekosky_pivot_replicate_donors_by_light_value_counts_df["x-axis"] = "Dekosky"
        dekosky_pivot_replicate_donors_by_light_value_counts_df = (
            dekosky_pivot_replicate_donors_by_light_value_counts_df.pivot(
                index="x-axis", columns="level_2", values="v_call_top_light_gene"
            ).loc[:, light_genes]
        )
        return dekosky_pivot_replicate_donors_by_light_value_counts_df
