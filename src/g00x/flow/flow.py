"""ABV - always be validating"""
import dataclasses
import json
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from g00x.data import Data
from g00x.validations.flow_validation import validate_g00x_box

# these are the unique indexable columns for each flow data file for a single Sorting experiment
index_flow_cols: list[str] = [
    "run_purpose",
    "run_date",
    "sort_id",
    "ptid",
    "group",
    "weeks",
    "visit_id",
    "probe_set",
    "sample_type",
    "sort_software_dv",
    "sort_file_type",
    "sample_tube",
    "gate",
    "phenotype",
    "value_type",
    "extention",
    "sort_pool",
    "hashtag",
]

# these are the unique indexable columns for each flow data file for a pre_screen experiment
index_flow_cols_ps: list[str] = [
    "run_purpose",
    "run_date",
    "sort_id",
    "ptid",
    "group",
    "weeks",
    "visit_id",
    "probe_set",
    "sample_type",
    "sort_software_dv",
    "sort_file_type",
    "sample_tube",
    "gate",
    "phenotype",
    "value_type",
    "extention",
]

sort_unique_columns: list[str] = [
    "run_purpose",
    "run_date",
    "sort_id",
    "ptid",
    "group",
    "weeks",
    "visit_id",
    "probe_set",
    "sample_type",
    "sort_software_dv",
    "sort_file_type",
    "sample_tube",
    "gate",
    "phenotype",
    "value_type",
    "sort_pool",
    "hashtag",
]


presort_unique_columns: list[str] = [
    "run_purpose",
    "run_date",
    "sort_id",
    "ptid",
    "group",
    "weeks",
    "visit_id",
    "probe_set",
    "sample_type",
    "sort_software_dv",
    "sort_file_type",
    "sample_tube",
    "gate",
    "phenotype",
    "value_type",
]

# has na columns only for clinical
has_na_columns: list[str] = ["branch", "easy_name", "notes"]


def sum_up_file_subsets(sub_df: pd.DataFrame) -> pd.Series:
    """
    A groupby function to return summed up gates from, _a,_b_c..etc files

    When the flow samples stop in the middle of the run,there are multiple .csv files output for each "sub" run. We need to groupby and sum up the gates to get the total values

    Parameters
    ----------
    sub_df : pd.DataFrame
        A grouped dataframe

    Returns
    -------
    pd.Series

    Raises
    ------
    ValueError
        If the non-grouped gates are different than expected
    ValueError
        If the gates don't have same branches
    ValueError
        If the gates don't have same notes
    """

    remaning: list[str] = ["branch", "easy_name", "file_path", "file_subset", "notes", "value"]
    if sorted(list(sub_df.columns.difference(index_flow_cols))) != sorted(remaning):  # type: ignore
        raise ValueError("Columns are not as expected")

    # Add the data from the file paths and aggregating the file subsets
    data: dict[str, list[str] | str | int | None] = {
        "file_path": sub_df["file_path"].to_list(),
        "file_subset": sorted(sub_df["file_subset"].to_list()),
        "value": sub_df["value"].sum(),
    }
    # we have to deal with the columns that could potentially contain NAs seperately since they can't join on an index
    for col in has_na_columns:
        value = list(set(sub_df[col]))
        data[col] = value[0]

        if len(value) > 1:
            error_df = sub_df[sub_df[col].isin(value)]
            raise ValueError(f"more than one {col} in aggreagetion {value} from rows {error_df.to_dict('records')}")

    return pd.Series(data=data)


def find_skip_rows(file: str | Path) -> int:
    """Find the number of rows to skip in a file until we get at longform data"""
    for i, x in enumerate(open(file).readlines()):
        if x.split(",")[0].lower() == "population":
            return i
    raise ValueError("No rows found that start with 'Population/population'")


@cache
def read_csv(file: Path | str, skiprows: int) -> pd.DataFrame:
    return pd.read_csv(file, skiprows=skiprows)  # pyright: reportUnknownMemberType=false


@dataclasses.dataclass
class SortPopulationQuery:
    """Data class used as a single flow gate query defined by the gate name and the parent name"""

    gate: str
    phenotype: str
    branch: str
    value_type: str = "count"
    easy_name: str | None = None
    notes: str | None = None

    def assign(self, file: Path | str, skip_rows: int) -> pd.Series:
        """Change the datafile into a pandas series"""
        population_csv = read_csv(file=file, skiprows=skip_rows)  # pyright: reportUnknownMemberType=false
        parsable_series = population_csv.iloc[:, :3].set_index("Population")["#Events"]
        value = parsable_series.loc[self.gate]
        series = pd.Series(self.__dict__)
        series["value"] = value
        return series


@dataclasses.dataclass
class SortPopulationFrequency:
    """Data class used as a single flow frequency that the division of two populations to make a frequency"""

    gate_numerator: str
    gate_parent_numerator: str
    gate_denominator: str
    gate_parent_denominator: str
    easy_name: str
    verbose_name: str
    value_type: str = "frequency"

    def assign(self, file: Path | str, skip_rows: int) -> pd.Series:
        """Change the datafile into a pandas series"""
        population_csv = read_csv(file=file, skiprows=skip_rows)  # pyright: reportUnknownMemberType=false
        parsable_series = population_csv.iloc[:, :3].set_index(["Population", "Parent Name"])["#Events"]
        value_numerator = parsable_series.loc[(self.gate_numerator, self.gate_parent_numerator)]
        value_denominator = parsable_series.loc[(self.gate_denominator, self.gate_parent_denominator)]
        series = pd.Series(
            {
                "gate": "/".join([self.gate_numerator, self.gate_denominator]),
                "gate_parent": "/".join([self.gate_parent_numerator, self.gate_parent_denominator]),
                "easy_name": self.easy_name,
                "verbose_name": self.verbose_name,
                "value_type": self.value_type,
                "value": value_numerator / value_denominator,
            }
        )
        return series


@dataclasses.dataclass
class Gates:
    """A class to hold all the gates from a flow experiment. Is easily serializable to to and from json"""

    gates: list[SortPopulationQuery]

    def to_json(self, file_name: str | Path) -> bool:
        """Dump to json"""
        output_dict = dataclasses.asdict(self)
        json.dump(output_dict, open(file_name, "w"), indent=True, sort_keys=True)
        return True

    @staticmethod
    def from_json(file_name: str | Path) -> "Gates":
        """Load from json"""
        initial_dict = json.load(open(file_name))
        gates = []
        for key in initial_dict:
            gates.append(SortPopulationQuery(**key))
        return Gates(
            gates=gates,
        )


@dataclasses.dataclass
class Frequencies:
    """A class to hold all the frequencies we intend to measure. Easily serializable"""

    frequencies: list[SortPopulationFrequency]

    def to_json(self, file_name: str | Path) -> bool:
        """Dump to json"""
        output_dict = dataclasses.asdict(self)
        json.dump(output_dict, open(file_name, "w"), indent=True, sort_keys=True)
        return True

    @staticmethod
    def from_json(file_name: str | Path) -> "Frequencies":
        """Load from json"""
        initial_dict = json.load(open(file_name))
        frequencies = []
        for key in initial_dict:
            frequencies.append(SortPopulationFrequency(**key))
        return Frequencies(
            frequencies=frequencies,
        )


@dataclasses.dataclass
class SortPopulationQueries:
    """Group of SortPopulations and the initial sort dataframe that will be expanded. This class will actually parse csv and hold the data"""

    individual_sort_queries: Gates
    initial_sort_dataframe: pd.DataFrame
    preclinical: bool

    class Config:
        arbitrary_types_allowed = True

    def get_count_dataframe(self) -> pd.DataFrame:
        """For each file and each gate, build them backup into a combined dataframe"""
        new_structure: list[pd.DataFrame] = []

        def assign_gate(x: Any, gate: SortPopulationQuery | SortPopulationFrequency) -> Any:
            """Assign a gate to a file"""
            skip_rows = find_skip_rows(x)
            return gate.assign(x, skip_rows)

        # @todo make this multiprocessing
        for gate in self.individual_sort_queries.gates:
            new_structure.append(
                self.initial_sort_dataframe.join(
                    self.initial_sort_dataframe["file_path"].apply(assign_gate, args=(gate,))
                )
            )

        return pd.concat(new_structure).reset_index(drop=True).astype({"file_path": str})


def parse_flow_data(data: Data, folder: str | Path) -> pd.DataFrame:
    """Main parse function for flow data for g002

    Parameters
    ----------
    folder : str | Path
        The Box/G002 folder to parse

    Returns
    -------
    pd.DataFrame
        A fully parsed dataframe. Each row represents a flow measurement or frequency with other requisite metadata
    """
    if not Path(folder).exists():
        raise FileNotFoundError(f"Folder {folder} does not exist")

    # must revalidate the folder because it will have our folders parsed in a method
    validation = validate_g00x_box(Path(folder))

    # get all the clinical population sort files
    all_population_sort_files = validation.get_population_sort_files()

    # get prescreen too
    pre_screen_sort_files = validation.get_prescreen_population_sort_files()

    # load gates from json
    clinical_gates_pbmc = Gates.from_json(data.data_paths.g002_pbmc_gates)

    # these are the lfna gates and probably need to be updated
    clinical_gates_lfna = Gates.from_json(data.data_paths.g002_lfna_gates)

    # change them into a dataframe
    sorting_dataframe = all_population_sort_files.get_dataframe()
    prescreen_dataframe = pre_screen_sort_files.get_dataframe()
    sorting_dataframe_pbmc = sorting_dataframe.query("sample_type == 'PBMC'")
    sorting_dataframe_lfna = sorting_dataframe.query("sample_type == 'LFNA'")

    if not prescreen_dataframe.empty:
        prescreen_dataframe_pbmc = prescreen_dataframe.query("sample_type == 'PBMC'")
        prescreen_dataframe_lfna = prescreen_dataframe.query("sample_type == 'LFNA'")
    else:
        prescreen_dataframe_pbmc = pd.DataFrame()
        prescreen_dataframe_lfna = pd.DataFrame()

    ### PBMC Clinical ###
    if not sorting_dataframe_pbmc.empty:
        flow_model_pbmc = SortPopulationQueries(
            individual_sort_queries=clinical_gates_pbmc,
            initial_sort_dataframe=sorting_dataframe_pbmc,
            preclinical=False,
        )
        clinical_count_df_pbmc = flow_model_pbmc.get_count_dataframe()
        summed_clinical_count_df_pbmc = (
            clinical_count_df_pbmc.groupby(index_flow_cols).apply(sum_up_file_subsets).reset_index()
        )

        # make sure we have no duplicate entries
        if summed_clinical_count_df_pbmc.groupby(sort_unique_columns).size().max() > 1:
            raise ValueError(
                f"Multiple files per sort entries {summed_clinical_count_df_pbmc.groupby(sort_unique_columns).size().max()}"
            )

    else:
        summed_clinical_count_df_pbmc = pd.DataFrame()

    ### PBMC Prescreen ###
    if not prescreen_dataframe_pbmc.empty:
        # perhaps a different flow model
        prescreen_flow_model_pbmc = SortPopulationQueries(
            individual_sort_queries=clinical_gates_pbmc,
            initial_sort_dataframe=prescreen_dataframe_pbmc,
            preclinical=True,
        )
        # combines dataframe from inital flow with all the values we need - We have not assigned frequencies yet
        prescreen_count_df_pbmc = prescreen_flow_model_pbmc.get_count_dataframe()
        try:
            summed_preclinical_count_df_pbmc = (
                prescreen_count_df_pbmc.groupby(index_flow_cols_ps).apply(sum_up_file_subsets).reset_index()
            )

            # make sure we have no duplicate entries
            if summed_preclinical_count_df_pbmc.groupby(presort_unique_columns).size().max() > 1:
                raise ValueError(
                    f"Multiple files per presort entries {summed_preclinical_count_df_pbmc.groupby(presort_unique_columns).size().max()}"
                )
        except ValueError as e:
            summed_preclinical_count_df_pbmc = pd.DataFrame()

    else:
        # make an empty one so we can combine it
        summed_preclinical_count_df_pbmc = pd.DataFrame()

    ### LFNA Clinical ###
    if not sorting_dataframe_lfna.empty:
        flow_model_lfna = SortPopulationQueries(
            individual_sort_queries=clinical_gates_lfna,
            initial_sort_dataframe=sorting_dataframe_lfna,
            preclinical=False,
        )
        clinical_count_df_lfna = flow_model_lfna.get_count_dataframe()
        summed_clinical_count_df_lfna = (
            clinical_count_df_lfna.groupby(index_flow_cols).apply(sum_up_file_subsets).reset_index()
        )

        # make sure we have no duplicate entries
        if summed_clinical_count_df_lfna.groupby(sort_unique_columns).size().max() > 1:
            raise ValueError(
                f"Multiple files per sort entries {summed_clinical_count_df_lfna.groupby(sort_unique_columns).size().max()}"
            )

    else:
        summed_clinical_count_df_lfna = pd.DataFrame()

    ### LFNA Prescreen ###
    if not prescreen_dataframe_lfna.empty:
        # perhaps a different flow model
        prescreen_flow_model_lfna = SortPopulationQueries(
            individual_sort_queries=clinical_gates_lfna,
            initial_sort_dataframe=prescreen_dataframe_lfna,
            preclinical=True,
        )
        # combines dataframe from inital flow with all the values we need - We have not assigned frequencies yet
        prescreen_count_df_lfna = prescreen_flow_model_lfna.get_count_dataframe()
        summed_preclinical_count_df_lfna = (
            prescreen_count_df_lfna.groupby(index_flow_cols_ps).apply(sum_up_file_subsets).reset_index()
        )

        # make sure we have no duplicate entries
        if summed_preclinical_count_df_lfna.groupby(presort_unique_columns).size().max() > 1:
            raise ValueError(
                f"Multiple files per presort entries {summed_preclinical_count_df_lfna.groupby(presort_unique_columns).size().max()}"
            )

    else:
        # make an empty one so we can combine it
        summed_preclinical_count_df_lfna = pd.DataFrame()

    # now combine all dataframes
    combined_dataframe = (
        pd.concat(
            [
                summed_preclinical_count_df_pbmc,
                summed_clinical_count_df_pbmc,
                summed_clinical_count_df_lfna,
                summed_preclinical_count_df_lfna,
            ]
        )
        .reset_index(drop=True)
        .replace([None], [np.nan])
    )

    return combined_dataframe
