import json
import logging
from dataclasses import dataclass
from pathlib import Path
from shutil import which

import pandas as pd
from pydantic import BaseModel, validator

logger = logging.getLogger()


def find_cellranger() -> str | None:
    which_cellranger = which("cellranger")
    if which_cellranger:
        logger.info("Cellranger found at: %s", which_cellranger)
        return which_cellranger
    else:
        raise FileNotFoundError("Cellranger not found in path. Please have Cellranger in path")


def find_bcl2fastq() -> str | None:
    which_bcl2fastq = which("bcl2fastq")
    if which_bcl2fastq:
        logger.info("bcl2fastq found at: %s", which_bcl2fastq)
        return which_bcl2fastq
    else:
        raise FileNotFoundError("bcl2fastq not found in path. Please have bcl2fastq in path")


class DataPaths(BaseModel):
    data_base_path: Path = Path(__file__).parent / Path("data")
    g001_sequences: Path = data_base_path / Path("g001/mutational_analysis.feather")
    g001_flow_and_seq: Path = data_base_path / Path("g001/flow_and_seq.csv")
    g002_pbmc_gates: Path = data_base_path / Path("pbmc_gates.json")
    g002_lfna_gates: Path = data_base_path / Path("lfna_gates.json")
    g003_pbmc_gates: Path = data_base_path / Path("pbmc_gates_g003.json")
    g002_frequency_measures: Path = data_base_path / Path("frequency_measures.json")
    g003_frequency_measures: Path = data_base_path / Path("frequency_measures_g003.json")
    g002_pub_ids_path: Path = data_base_path / Path("g002_pubids.xlsx")
    # g003_pub_ids_path: Path = data_base_path / Path("g003/g003_pubids.xlsx")
    hto_gates: Path = data_base_path / Path("hto_gates.csv")
    g003_hto_gates: Path = data_base_path / Path("g003_hto_gates.csv")
    vdj_reference: Path = data_base_path / Path("vdj_reference")
    vh12_reference_mabs_path: Path = data_base_path / Path("vrc1-2_mabs_extended_airr.feather")
    cotrell_focus_path: Path = data_base_path / Path("cottrell.json")
    personalized_vh12: Path = data_base_path / Path("personalize.csv")
    g003_personalized_vh12: Path = data_base_path / Path("g003/genotype.csv")
    g003_visit_id_2_week: Path = data_base_path / Path("g003/visit_id2week.json")
    g003_ptid_prefix_2_group: Path = data_base_path / Path("g003/ptid2group.json")

    @validator(
        "*",
        pre=True,
        always=True,
    )
    def validate_path(cls, v: str) -> str:
        if not Path(v).exists():
            raise ValueError(f"{v} does not exist")
        return v


@dataclass
class PlotParameters:
    pass


class Data:
    def __init__(self) -> None:
        self.data_paths = DataPaths()
        self.plot_parameters = PlotParameters()
        self.cellranger_path = ""
        self.genome_reference: Path | None = None

    def get_g001_sequences(self) -> pd.DataFrame:
        """The  g001 sequnces that can be used compariatiely in this package"""
        return (
            pd.read_feather(self.data_paths.g001_sequences)
            .query("Plate_Type_heavy=='Probe Specific'")
            .reset_index(drop=True)
        )

    def get_g001_flow_and_seq(self) -> pd.DataFrame:
        "The g001 flow and sequence analysis that can be use compariteiely"
        return pd.read_csv(self.data_paths.g001_flow_and_seq, index_col=0)

    def get_pbmc_gates(self) -> list[dict[str, str]]:
        """Get the gates for the PBMC sort for G002"""
        return json.load(open(self.data_paths.g002_pbmc_gates))

    def get_lfna_gates(self) -> list[dict[str, str]]:
        """Get the gates for the lfna sort for G002"""
        return json.load(open(self.data_paths.g002_lfna_gates))

    def get_pbmc_gates_g003(self) -> list[dict[str, str]]:
        """Get the gates for the PBMC sort for G003"""
        return json.load(open(self.data_paths.g003_pbmc_gates))

    def get_frequency_measures(self) -> list[dict[str, str]]:
        """The frequency measures that can be used in the flow package"""
        return json.load(open(self.data_paths.g002_frequency_measures))

    def get_frequency_measures_g003(self) -> list[dict[str, str]]:
        """The frequency measures that can be used in the flow package"""
        return json.load(open(self.data_paths.g003_frequency_measures))

    def get_hto_gates(self) -> pd.DataFrame:
        """Get the HTO gates for the CSO part of the pipeline"""
        return pd.read_csv(self.data_paths.hto_gates, index_col=0)

    def get_g003_hto_gates(self) -> pd.DataFrame:
        """Get the HTO gates for the G003 CSO part of the pipeline"""
        return pd.read_csv(self.data_paths.g003_hto_gates).set_index("biolegend_name")

    def set_cellranger_path(self, path: str) -> None:
        """set path to cellranger if user specifies"""
        if Path(path).exists():
            self.cellranger_path = str(Path(path).absolute())
        else:
            raise ValueError(f"{path} does not exist")

    def set_genome_reference(self, path: str) -> None:
        """set path to genome reference if user specifies"""
        if not Path(path).exists():
            raise ValueError(f"{path} does not exist")
        self.genome_reference = Path(path).absolute()

    def get_human_genome_ref(self) -> str:
        if not self.genome_reference:
            self.genome_reference = Path("/usr/local/cellranger-6.1.2/refdata-gex-GRCh38-2020-A")
        if not self.genome_reference.exists():
            raise FileNotFoundError(
                f"Genome reference not found {self.genome_reference}, set with flags --genome-reference"
            )
        return str(self.genome_reference)

    def get_cellranger_path(self) -> str | None:
        """get the path to cellranger"""
        if not self.cellranger_path:
            self.cellranger_path = find_cellranger()
        return self.cellranger_path

    def get_bcl2fastq_path(self) -> str | None:
        """
        find bcl2fastq for insurance purposes

        we can't actually set this anywhere but we have to make sure its in path
        """
        bcl2fastq_path: str | None = find_bcl2fastq()
        return bcl2fastq_path

    def get_vdj_path(self) -> str:
        """get reference to the path for vdj"""
        return str(self.data_paths.vdj_reference)

    def get_g002_pubids_lookup(self) -> dict[str, str]:
        """From VISC, get the ptid mapped to pubids"""
        # from visc
        order_id = pd.read_excel(self.data_paths.g002_pub_ids_path)

        # correct VISC by adding G002
        # order_id["ptid"] = "G002" + order_id["ptid"].astype(str)
        order_id["ptid"] = order_id["pubID"]

        # turn into lookup map
        lookup_ids = dict(zip(order_id["ptid"].to_list(), order_id["pubID"].to_list()))
        return lookup_ids

    # def get_g003_pubids_lookup(self) -> dict[str, str]:
    #     """From VISC, get the ptid mapped to pubids"""
    #     # from visc
    #     df = pd.read_excel(self.data_paths.g003_pub_ids_path, dtype=str)

    #     # correct VISC by adding G003
    #     # df["ptid"] = "G003-" + df["ptid"].str[:2] + "-" + df["ptid"].str[2:]
    #     df['ptid'] = df['pubID']

    #     # turn into lookup map
    #     lookup_ids = df.set_index("ptid").to_dict()["pubID"]
    #     return lookup_ids

    def get_vh12_reference_airr_table(self) -> pd.DataFrame:
        """get the path to the reference airr table"""
        return pd.read_feather(self.data_paths.vh12_reference_mabs_path)

    def get_cotrell_focus(self) -> dict[str, list[str]]:
        """get the cotrell focus"""
        return json.load(open(self.data_paths.cotrell_focus_path))

    def get_personalized_vh12(self) -> pd.DataFrame:
        """get the personalized VH12"""
        return pd.read_csv(self.data_paths.personalized_vh12)

    def get_g003_personalized_vh12(self) -> pd.DataFrame:
        """get the personalized VH12"""
        return pd.read_csv(self.data_paths.g003_personalized_vh12)

    def get_long_name_sort_order(self) -> list[str]:
        return [
            "run_purpose",
            "run_date",
            "pubID",
            "ptid",
            "group",
            "weeks",
            "visit_id",
            "probe_set",
            "sample_type",
            "Lymphocytes",
            "B cells",
            "Singlets",
            "Dump-",
            "IgD+ B cells",
            "IgD+/Antigen++ B cells",
            "IgD+/Antigen++/KO- B cells",
            "IgD+/KO- B cells",
            "IgD+/KO-/Antigen++ B cells",
            "IgD- B cells",
            "IgD-/Antigen++ B cells",
            "IgD-/Antigen++/KO- B cells",
            "IgD-/IgG-/IgA- B cells",
            "IgD-/IgG-/IgA-/IgM+ B cells",
            "IgD-/IgG-/IgA-/IgM+/Antigen++ B cells",
            "IgD-/IgG-/IgA-/IgM+/Antigen++/KO- B cells",
            "IgD-/IgG-/IgA-/IgM+/KO- B cells",
            "IgD-/IgG-/IgA-/IgM+/KO-/Antigen++ B cells",
            "IgD-/IgG-/IgM- B cells",
            "IgD-/IgG-/IgM-/IgA+ B cells",
            "IgD-/IgG-/IgM-/IgA+/Antigen++ B cells",
            "IgD-/IgG-/IgM-/IgA+/Antigen++/KO- B cells",
            "IgD-/IgG-/IgM-/IgA+/KO- B cells",
            "IgD-/IgG-/IgM-/IgA+/KO-/Antigen++ B cells",
            "IgD-/IgM-/IgA- B cells",
            "IgD-/IgM-/IgA-/IgG+/Antigen++ B cells",
            "IgD-/IgM-/IgA-/IgG+/Antigen++/KO- B cells",
            "IgD-/IgM-/IgA-/IgG+/KO- B cells",
            "IgD-/IgM-/IgA-/IgG+/KO-/Antigen++ B cells",
            "IgD-/KO-/Antigen++ (sorted) B cells",
            "IgD-IgM-/IgA-/IgG+ B cells",
            "IgD-KO- B cells",
            "CD19+ B cells",
            "CD19+/CD20- B cells",
            "CD19+/CD20-/CD27+/CD38+ PB cells",
            "CD19+/CD20-/CD27+/CD38+/Antigen++ PB cells",
            "CD19+/CD20-/CD27+/CD38+/Antigen++/KO- PB cells",
            "CD19+/CD20-/CD27+/CD38+/KO- PB cells",
            "CD19+/CD20-/CD27+/CD38+/KO-/Antigen++ PB cells",
            "Number of IGD- sequences that are VRC01-class",
            "Number of IGHA sequences that are VRC01-class",
            "Number of IGHA1*01 sequences that are VRC01-class",
            "Number of IGHA2*01 sequences that are VRC01-class",
            "Number of IGHD sequences that are VRC01-class",
            "Number of IGHD*02 sequences that are VRC01-class",
            "Number of IGHG sequences that are VRC01-class",
            "Number of IGHG1*01 sequences that are VRC01-class",
            "Number of IGHG2*01 sequences that are VRC01-class",
            "Number of IGHG3*01 sequences that are VRC01-class",
            "Number of IGHG4*01 sequences that are VRC01-class",
            "Number of IGHM sequences that are VRC01-class",
            "Number of IGHM*01 sequences that are VRC01-class",
            "Number of None sequences that are VRC01-class",
            "Number of sequences",
            "Percent IgA^{+}KO^- among Ag^{++}",
            "Percent IgD^{-}KO^{-} among Ag^{++}",
            "Percent IgG^{+}KO^- among Ag^{++}",
            "Percent IgM{+}KO^- among Ag^{++}",
            "Percent PB{+}KO^- among Ag^{++}",
            "Percent antigen-specific among IgD^-",
            "Percent antigen-specific among IgG^{+}",
            "Percent antigen-specific among IgM",
            "Percent epitope-specific (KO^-Ag^{++}) among IgA^{+}",
            "Percent epitope-specific (KO^-Ag^{++}) among IgD^-",
            "Percent epitope-specific (KO^-Ag^{++}) among IgG^{+}",
            "Percent epitope-specific (KO^-Ag^{++}) among IgM",
            "Percent epitope-specific (KO^-Ag^{++}) among PB",
            "Percent of antigen-specific among IgA^{+}",
            "Percent of antigen-specific among PB",
            "Percent of IGD- sequences that are VRC01-class",
            "Percent of IGHA sequences that are VRC01-class",
            "Percent of IGHA1*01 sequences that are VRC01-class",
            "Percent of IGHA2*01 sequences that are VRC01-class",
            "Percent of IGHD sequences that are VRC01-class",
            "Percent of IGHD*02 sequences that are VRC01-class",
            "Percent of IGHG sequences that are VRC01-class",
            "Percent of IGHG1*01 sequences that are VRC01-class",
            "Percent of IGHG2*01 sequences that are VRC01-class",
            "Percent of IGHG3*01 sequences that are VRC01-class",
            "Percent of IGHG4*01 sequences that are VRC01-class",
            "Percent of IGHM sequences that are VRC01-class",
            "Percent of IGHM*01 sequences that are VRC01-class",
            "Percent of None sequences that are VRC01-class",
            "Percent of VRC01-class sequences among IgA",
            "Percent of VRC01-class sequences among IgD-",
            "Percent of VRC01-class sequences among IgG",
            "Percent of VRC01-class sequences among IgM",
        ]

    def get_g003_visit_id_2_week(self) -> dict[str, int]:
        """Get the visit id to week mapping"""
        return json.load(open(self.data_paths.g003_visit_id_2_week, encoding="utf-8"))

    def get_g003_ptid_prefix_2_group(self) -> dict[str, str]:
        """Get the ptid to group mapping"""
        return json.load(open(self.data_paths.g003_ptid_prefix_2_group, encoding="utf-8"))
