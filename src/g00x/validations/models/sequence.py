import datetime
import re
from functools import cache
from pathlib import Path
from typing import Literal

import pandas as pd
from pydantic import BaseModel, validator


@cache
def get_expected_sample_manifest_columns() -> list[str]:
    EXPECTED_SAMPLE_MANIFEST_COLUMNS = [
        "pool_number",
        "sorted_date",
        "vdj_sequencing_replicate",
        "cso_sequencing_replicate",
        "vdj_library_replicate",
        "cso_library_replicate",
        "bio_replicate",
        "vdj_index",
        "feature_index",
        "vdj_run_id",
        "cso_run_id",
    ]
    return EXPECTED_SAMPLE_MANIFEST_COLUMNS


class ModelError(ValueError):
    def __init__(self, accepted: str, given: str) -> None:
        self.accepted = accepted
        self.given = given

    def __str__(self) -> str:
        return f"\tGiven: {self.given}\n\tExpected: {self.accepted}"


class SampleManifestModel(BaseModel):
    illuimna_folder_name: list[str]
    path: Path

    @validator("path")
    def validate_path(cls, v: Path, values: dict[str, list[str]]) -> Path:
        illumina_folder_name: list[str] = values["illuimna_folder_name"]
        if v.name != "sample_manifest.csv":
            raise ValueError(f"{v} is not sample_manifest.csv...it must be named sample manifest csv no exceptions")

        # read in dataframe on validation
        _dataframe = pd.read_csv(v)

        # make sure we have just the expected sample manifest columns in order
        _expected_manifest_columns = get_expected_sample_manifest_columns()
        if _dataframe.columns.to_list() != _expected_manifest_columns:
            raise ValueError(
                f"{v} does not have the correct columns...needs to be {_expected_manifest_columns} in exact order"
            )

        # make sure pool number is int
        if _dataframe["pool_number"].dtype != int:
            raise ValueError(f"{v} pool_number column is not an int")

        # make sure we don't have more than 4 pools
        if (_dataframe["pool_number"] > 10).any():
            raise ValueError(f"{v} pool_number column has a value greater than 5")

        # Make sure the sorted date is before now and is in the right format
        try:
            if (
                pd.to_datetime(_dataframe["sorted_date"], format="%y%m%d", errors="raise") > datetime.datetime.now()  # type: ignore
            ).any():  # type: ignore
                raise ValueError(f"{v} sorted_date column has a date in the future")
        except ValueError:
            raise ValueError(f"{v} sorted_date column has a date in the future and has to meet the format %y%m%d")

        if not _dataframe["vdj_run_id"].isin(illumina_folder_name).all():
            raise ValueError(
                f"{_dataframe['vdj_run_id']} vdj_run_id column has a value that is not in the illumina folder name {illumina_folder_name}"
            )
        if not _dataframe["cso_run_id"].isin(illumina_folder_name).all():
            raise ValueError(
                f"{_dataframe['cso_run_id']} cso_run_id column has a value that is not in the illumina folder name {illumina_folder_name}"
            )
        return v

    def get_dataframe(self) -> pd.DataFrame:
        _df = pd.read_csv(self.path)
        _df["sorted_date"] = pd.to_datetime(_df["sorted_date"], format="%y%m%d")
        _df["pool_number"] = _df["pool_number"].apply(lambda x: "P" + str(x).zfill(2))

        # insert at first position the path of the sequence run_dir
        _df.insert(0, "run_dir_path", str(self.path.parent.absolute()))
        return _df


class IlluminaModel(BaseModel):
    run_dir: Path

    @validator("run_dir")
    def validate_run_dir(cls, v: Path) -> Path:
        if not v.is_dir() or not v.exists():
            raise ValueError(f"{v} is not a directory")
        if len(str(v).split("_")) != 4:
            raise ValueError(f"{v} does not have 4 parts seperated by _")
        if datetime.datetime.strptime(str(v.name).split("_")[0], "%y%m%d") > datetime.datetime.now():
            raise ValueError(f"{str(v).split('_')[0]} has a date in the future")
        return v


class RunModel(BaseModel):
    """The run0001 to run00019 model"""

    name: str

    @validator("name", pre=True)
    def validate_name(cls, v: str) -> str:
        if re.fullmatch(r"run000[0-9]|run001[0-9]|run002[0-9]|run003[0-9]", v):
            return v
        raise ModelError(given=v, accepted="run0001 to run0029")


class GlobusEndpointModel(BaseModel):
    name: Literal["globus_endpoint"]
    run: RunModel | None
