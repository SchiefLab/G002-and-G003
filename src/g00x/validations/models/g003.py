import re
from datetime import date, datetime
from typing import Any, Literal, Optional

from pydantic import BaseModel, Field, root_validator, validator


class ModelError(ValueError):
    def __init__(self, accepted: str, given: str) -> None:
        self.accepted = accepted
        self.given = given

    def __str__(self) -> str:
        return f"\tGiven: {self.given}\n\tExpected: {self.accepted}"


class Fields(BaseModel):
    """Base Model for fields for files. If needed to be modified, inherit from this class"""

    run_purpose: str
    run_date: date = Field(..., description="YYMMDD")
    sort_id: str = Field(regex=r"[A-Z]\d{1}")
    ptid: str = Field(regex="G00[1-5]-*[0-9]{2}($|-*[0-9]{3})")
    visit_id: str = Field(regex=r"V\d{1,3}")
    probe_set: Literal["eODGT8"]
    sample_type: Literal["PBMC", "FNA"]
    sort_software_dv: Literal["Chorus", "FlowJo"]
    sort_file_type: Literal["Primary", "Capture", "SortRpt", "Summary", "Data", "Stats"]
    sample_tube: str = Field(regex=r"na|T\d{1}")
    sort_pool_file_subset: str = Field(regex="Control|Presort|FNA|P[1-9][a-z]|P[1-9]")
    extention: Literal[".fcs", ".xlsx", ".png", ".pdf", ".jpg"]

    @root_validator(pre=True)
    def extract_dates(cls, values: dict[str, Any]) -> dict[str, Any]:
        for date_field in ["run_date"]:
            try:
                values[date_field] = datetime.strptime(values[date_field], "%y%m%d").date()
            except:
                values[date_field] = datetime.strptime(values[date_field], "%Y%m%d").date()
        return values

    # class Config:
    #     validate_assignment = True
    #     error_msg_templates = {"value_error.str.regex": "Given: {extra}"}


class DataFilesFromFlowJo(BaseModel):
    run_purpose: str
    run_date: date = Field(..., description="YYMMDD")
    sort_id: str
    ptid: str = Field(regex="G00[1-5]-[0-9]{2}($|-[0-9]{3})")
    visit_id: str = Field(regex=r"V\d{2,3}")
    extention: Literal[".wsp"]

    @root_validator(pre=True)
    def extract_dates(cls, values: dict[str, Any]) -> dict[str, Any]:
        for date_field in ["run_date"]:
            try:
                values[date_field] = datetime.strptime(values[date_field], "%y%m%d").date()
            except:
                values[date_field] = datetime.strptime(values[date_field], "%Y%m%d").date()
        return values


class DataFilesFromMelody(Fields):
    extention: Literal[".fcs"]  # type: ignore


class DataStats(Fields):
    extention: Literal[".xlsx", ".csv"]  # type: ignore


class ScreenshotsCounts(Fields):
    extention: Literal[".jpg", ".JPG", ".png", ".PNG"]  # type: ignore


class ScreenShotsMelodyStats(Fields):
    extention: Literal[".jpg", ".JPG", ".png", ".PNG"]  # type: ignore


class SortReports(Fields):
    extention: Literal[".pdf"]  # type: ignore


class ClinicalSamples(BaseModel):
    name: Literal["ClinicalSamples"]
    datafiles_from_flowjo: DataFilesFromFlowJo
    datafiles_from_melody: DataFilesFromMelody
    datastats: DataStats
    screenshots_counts: ScreenshotsCounts
    screenshots_melody_stats: ScreenShotsMelodyStats
    sortreports: SortReports


class Sort(BaseModel):
    name: str
    run_date: Optional[date] = None
    upload_date: Optional[date] = None

    @root_validator(pre=True)
    def extract_dates(cls, values: dict[str, Any]) -> dict[str, Any]:
        if not re.fullmatch(
            r"Sort_RunDate\d{2}([0][1-9]|[1][0-2])([0-3][0-9])_UploadDate\d{2}([0][1-9]|[1][0-2])([0-3][0-9])",
            values["name"],
        ):
            raise ModelError(
                given=values["name"], accepted="Sort name must be in the format Sort_RunDateYYMMDD_UploadDateYYMMDD"
            )
        try:
            run_date_string = values["name"].split("_")[1].split("RunDate")[1]
            try:
                values["run_date"] = datetime.strptime(run_date_string, "%y%m%d").date()
            except:
                values["run_date"] = datetime.strptime(run_date_string, "%Y%m%d").date()
            upload_date_string = values["name"].split("_")[2].split("UploadDate")[1]
            values["upload_date"] = datetime.strptime(upload_date_string, "%y%m%d").date()
        except (IndexError, ValueError):
            raise ModelError(
                given=values["name"], accepted="Sort name must be in the format Sort_RunDateYYMMDD_UploadDateYYMMDD"
            )
        return values


class Sorts(BaseModel):
    name: Literal["Sorts"]
    sort: Sort


class G00X(BaseModel):
    name: str
    sorts: Optional[Sorts] = None

    @validator("name", pre=True)
    def extract_scheme(cls, v: str) -> str:
        if re.fullmatch("[gG]00[1-5x](_Scheme_Example)?", v):
            return v
        raise ModelError(given=v, accepted="Root folder name must named G001-G005")
