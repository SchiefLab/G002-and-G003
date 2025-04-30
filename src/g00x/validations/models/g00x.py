import json
import re
from datetime import date, datetime
from pathlib import Path
from typing import Any, Literal, Optional

import pandas as pd
from pydantic import BaseModel, Field, root_validator, validator


class VisitIDs:
    """A simple class to assign available visit IDs to based on groups."""

    group_1_available_timepoints: list[str] = ["91", "150", "160", "200", "250", "270", "290", "320"]
    group_1_available_weeks: list[str] = ["-5", "3", "4", "8", "11", "16", "19", "24"]
    group_2_available_timepoints: list[str] = ["91", "150", "160", "200", "250", "270", "320"]
    group_2_available_weeks: list[str] = ["-5", "3", "4", "8", "11", "16", "24"]
    group_3_available_timepoints: list[str] = ["91", "150", "160", "401", "500", "550", "570"]
    group_3_available_weeks: list[str] = ["-5", "3", "4", "8", "16", "19", "24"]
    group_4_available_timepoints: list[str] = ["91", "150", "160", "600"]
    group_4_available_weeks: list[str] = ["-5", "3", "4", "8"]
    group_lookup: dict[int, list[str]] = {
        1: group_1_available_timepoints,
        2: group_2_available_timepoints,
        3: group_3_available_timepoints,
        4: group_4_available_timepoints,
    }
    group_week_lookup: dict[int, dict[str, str]] = {
        1: {f"V{str(i).zfill(3)}": k for i, k in zip(group_1_available_timepoints, group_1_available_weeks)},
        2: {f"V{str(i).zfill(3)}": k for i, k in zip(group_2_available_timepoints, group_2_available_weeks)},
        3: {f"V{str(i).zfill(3)}": k for i, k in zip(group_3_available_timepoints, group_3_available_weeks)},
        4: {f"V{str(i).zfill(3)}": k for i, k in zip(group_4_available_timepoints, group_4_available_weeks)},
    }

    def get_timepoints_by_group(self, group: int) -> list[str]:
        "Given a group, return the available timepoints with attached 'V', e.g. 'V91'"
        return [f"V{str(i).zfill(3)}" for i in self.group_lookup[group]]

    def get_week_by_timepoint_and_group(self, group: int, timepoint: str) -> str:
        "Given a group and a timepoint, return the week"
        return self.group_week_lookup[group][timepoint]


class ParticipantIDs:
    def __init__(self) -> None:
        self.file_path: Path = Path(__file__).parent.parent / Path("data/enrollment.json")
        self.enrollment_data: list[dict[str, int | str]] = json.load(open(self.file_path))

    def get_enrollment_data(self) -> list[dict[str, int | str]]:
        return self.enrollment_data

    def get_ptids(self) -> list[str | int]:
        return list(pt["ptid"] for pt in self.enrollment_data)

    def is_ptid_in_enrollment(self, ptid: str) -> bool:
        return any(ptid == pt["ptid"] for pt in self.enrollment_data)

    def get_group_by_ptid(self, ptid: str | int) -> int | str:
        for pt in self.enrollment_data:
            if pt["ptid"] == ptid:
                return pt["group"]
        raise ValueError(f"Participant {ptid} not found in enrollment data.")


class ModelError(ValueError):
    def __init__(self, accepted: str, given: str) -> None:
        self.accepted = accepted
        self.given = given

    def __str__(self) -> str:
        return f"\tGiven: {self.given}\n\tExpected: {self.accepted}"


class PrescreenSharedFileFieldsModel(BaseModel):
    run_purpose: Literal["PreS", "Sort"]
    run_date: date | str = Field(..., description="YYMMDD")
    sort_id: str
    ptid: str
    visit_id: str
    probe_set: Literal["eODGT8", "Cg28v2"]
    sample_type: Literal["PBMC", "LFNA"]
    sort_software_dv: Literal["DV"]
    sort_file_type: Literal["Primary", "Capture", "SortRpt", "Summary"]
    sample_tube: str
    file_subset: str
    extention: Literal[".fcs", ".csv", ".png", ".PNG", ".pdf"]

    # should not assign this and let the validator do it
    group: Optional[int] = None
    weeks: Optional[str] = None

    @validator("sort_id", pre=True)
    def validate_sort_id(cls, v: str) -> str:
        if re.fullmatch(r"^[a-zA-Z0-9]{3}$", v):
            return v
        raise ModelError(given=v, accepted="Sort ID must be 3 characters/numbers")

    @validator("ptid", pre=True)
    def validate_trial_id_site_id_donor_id(cls, v: str) -> str:
        """Validate the ptid from the enrollment file."""
        participant_ids_instance = ParticipantIDs()
        if participant_ids_instance.is_ptid_in_enrollment(v):
            return v
        else:
            raise ModelError(given=v, accepted=",".join(map(str, participant_ids_instance.get_ptids())))

    @validator("visit_id", pre=True)
    def validate_visit_id(cls, v: str, values: dict[str, int | str]) -> str:
        """Validate visit ID based on group

        First look up the already validated ptid. Then look at the assosciated group.

        This function will also assign the group based on the ptid.

        Then check if the visit ID falls in that group"""
        if "ptid" not in values:
            raise ValueError("ptid must be validated before visit_id")
        ptid = values["ptid"]
        associated_group = int(ParticipantIDs().get_group_by_ptid(ptid))
        values["group"] = associated_group
        timepoints_for_group = VisitIDs().get_timepoints_by_group(associated_group)
        if v in timepoints_for_group:
            values["weeks"] = VisitIDs().get_week_by_timepoint_and_group(associated_group, v)
            return v
        raise ModelError(
            given=v,
            accepted=f"Visit ID must be one of the following: {timepoints_for_group} for group {associated_group}",
        )

    @validator("file_subset", pre=True)
    def validate_file_subset(cls, v: str) -> str:
        if re.fullmatch("[a-j]", v):
            return v
        else:
            raise ModelError(given=v, accepted="file_subset must be a single letter between 'a' and 'j'; i.e. a")

    @validator("run_date", pre=True)
    def extract_date(cls, v: str) -> date:
        "Extract the date from the file name as a date object"
        return datetime.strptime(v, "%y%m%d").date()

    # @validator("group", always=True)
    # def validate_group(cls, v: int, values: dict[str, int | str]) -> int | str:
    #     if v is None:
    #         print(values)
    #         return values["group"]
    #     else:
    #         return v

    # @validator("weeks", always=True)
    # def validate_weeks(cls, v: int, values: dict[str, int | str]) -> int | str:
    #     if v is None:
    #         return values["weeks"]
    #     else:
    #         return v


class ClinicalSharedFileFieldsModel(BaseModel):
    run_purpose: Literal["PreS", "Sort"]
    run_date: date = Field(..., description="YYMMDD")
    sort_id: str
    ptid: str
    visit_id: str
    probe_set: Literal["eODGT8", "Cg28v2"]
    sample_type: Literal["PBMC", "LFNA"]
    hashtag: str
    sort_software_dv: Literal["DV"]
    sort_file_type: Literal["Primary", "Capture", "SortRpt", "Summary"]
    sample_tube: str
    sort_pool: str
    file_subset: str
    extention: Literal[".fcs", ".csv", ".png", ".PNG", ".pdf"]

    # should not assign this and let the validator do it
    group: Optional[int] = None
    weeks: Optional[str] = None

    @validator("sort_id", pre=True)
    def validate_sort_id(cls, v: str) -> str:
        if re.fullmatch(r"^[a-zA-Z0-9]{3}$", v):
            return v
        raise ModelError(given=v, accepted="Sort ID must be 3 characters/numbers")

    @validator("ptid", pre=True)
    def validate_trial_id_site_id_donor_id(cls, v: str) -> str:
        """Validate the ptid from the enrollment file."""
        participant_ids_instance = ParticipantIDs()
        if participant_ids_instance.is_ptid_in_enrollment(v):
            return v
        else:
            raise ModelError(given=v, accepted=",".join(map(str, participant_ids_instance.get_ptids())))

    @validator("visit_id", pre=True)
    def validate_visit_id(cls, v: str, values: dict[str, int | str]) -> str:
        """Validate visit ID based on group

        First look up the already validated ptid. Then look at the assosciated group.

        This function will also assign the group based on the ptid.

        Then check if the visit ID falls in that group"""
        if "ptid" not in values:
            raise ModelError(given=v, accepted="ptid does not exist")
        ptid = values["ptid"]
        associated_group = int(ParticipantIDs().get_group_by_ptid(ptid))
        values["group"] = associated_group
        timepoints_for_group = VisitIDs().get_timepoints_by_group(associated_group)
        if v in timepoints_for_group:
            values["weeks"] = VisitIDs().get_week_by_timepoint_and_group(associated_group, v)
            return v
        raise ModelError(
            given=v,
            accepted=f"Visit ID must be one of the following: {timepoints_for_group} for group {associated_group}",
        )

    @validator("hashtag", pre=True)
    def validate_hashtag(cls, v: str) -> str:
        "Allow only HT01-HT15"
        if v not in [f"HT{i:02d}" for i in range(1, 16)] and v != "NA":
            raise ModelError(
                given=v, accepted="hashtag be one of HT01, HT02, HT03, HT04, HT05, HT06, HT07, HT08, HT09, HT10"
            )
        return v

    @validator("sample_tube", pre=True)
    def validate_sample_tube(cls, v: str) -> str:
        if re.fullmatch("T[1-9]", v):
            return v
        else:
            raise ModelError(given=v, accepted="sample_tube must start with T and be 1 digit long; i.e. T1")

    @validator("sort_pool", pre=True)
    def validate_sort_pool(cls, v: str) -> str:
        "Only allow P01, P02, P03, P04, P05, P06, P07, P08, P09, P10 or NA"
        if re.fullmatch("(P[0-9]{2}|NA)", v):
            return v
        else:
            raise ModelError(given=v, accepted="sort_pool must start with P and be 2 digits long; i.e. P01")

    @validator("file_subset", pre=True)
    def validate_file_subset(cls, v: str) -> str:
        if re.fullmatch("[a-j]", v):
            return v
        else:
            raise ModelError(given=v, accepted="file_subset must be a single letter between 'a' and 'j'; i.e. a")

    @validator("run_date", pre=True)
    def extract_date(cls, v: str) -> date:
        "Extract the date from the file name as a date object"
        return datetime.strptime(v, "%y%m%d").date()

    @validator("group", always=True)
    def validate_group(cls, v: int, values: dict[str, int | str]) -> int | str:
        if v is None:
            """Only return the group if it is not None. Otherwise return the value of weeks"""
            if "group" not in values:
                raise ModelError(given=v, accepted="group does not exist")
            return values["group"]
        else:
            return v

    @validator("weeks", always=True)
    def validate_weeks(cls, v: int, values: dict[str, int | str]) -> int | str:
        if v is None:
            if "weeks" not in values:
                raise ModelError(given=v, accepted="weeks does not exist")
            return values["weeks"]
        else:
            return v


class PackingSharedFileFieldsModel(BaseModel):
    run_purpose: Literal["PreS", "Sort"]
    run_date: date = Field(..., description="YYMMDD")
    sort_id: str
    leuokopack_id: str
    sample_type: Literal["PBMC", "LFNA"]
    hashtag: Literal["HT10"]
    sort_software_dv: Literal["DV"]
    sort_file_type: Literal["Primary", "Capture", "SortRpt", "Summary"]
    sample_tube: str
    sort_pool: str
    file_subset: str
    extention: Literal[".fcs", ".csv", ".png", ".PNG", ".pdf"]

    @validator("sort_id", pre=True)
    def validate_sort_id(cls, v: str) -> str:
        if re.fullmatch(r"^[a-zA-Z0-9]{3}$", v):
            return v
        raise ModelError(given=v, accepted="Sort ID must be 3 characters/numbers")

    @validator("sample_tube", pre=True)
    def validate_sample_tube(cls, v: str) -> str:
        if re.fullmatch("T[1-9]", v):
            return v
        else:
            raise ModelError(given=v, accepted="sample_tube must start with T and be 1 digit long; i.e. T1")

    @validator("sort_pool", pre=True)
    def validate_sort_pool(cls, v: str) -> str:
        "Only allow P01, P02, P03, P04, P05, P06, P07, P08, P09, P10 or NA"
        if re.fullmatch("(P[0-9]{2}|NA)", v):
            return v
        else:
            raise ModelError(given=v, accepted="sort_pool must start with P and be 2 digits long; i.e. P01")

    @validator("file_subset", pre=True)
    def validate_file_subset(cls, v: str) -> str:
        if re.fullmatch("[a-j]", v):
            return v
        else:
            raise ModelError(given=v, accepted="file_subset must be a single letter between 'a' and 'j'; i.e. a")

    @validator("run_date", pre=True)
    def extract_date(cls, v: str) -> date:
        "Extract the date from the file name as a date object"
        if isinstance(v, date):
            return v
        return datetime.strptime(v, "%y%m%d").date()


class ControlSharedFileFieldsModel(BaseModel):
    run_purpose: Literal["PreS", "Sort"]
    sort_date: date = Field(..., description="YYMMDD")
    sort_id: str
    sample_id: str
    hashtag: str
    sort_software_dv: Literal["DV"]
    sort_file_type: Literal["Primary", "Capture", "SortRpt", "Summary"]
    sample_tube: str
    sort_pool: str
    file_subset: str
    extention: Literal[".fcs", ".csv", ".png", ".PNG", ".pdf"]

    @validator("sort_date", pre=True)
    def extract_date(cls, v: str) -> date:
        return datetime.strptime(v, "%y%m%d").date()

    @validator("sort_id", pre=True)
    def validate_sort_id(cls, v: str):
        if re.fullmatch(r"^[a-zA-Z0-9]{3}$", v):
            return v
        raise ModelError(given=v, accepted="Sort ID must be 3 characters/numbers")

    @validator("sample_id", pre=True)
    def validate_sample_id(cls, v: str) -> str:
        if re.fullmatch("[a-zA-Z]{2,}", v):
            return v
        else:
            raise ModelError(given=v, accepted="Sample ID must be 2 characters long; i.e. AA")

    @validator("hashtag", pre=True)
    def validate_hashtag(cls, v: str) -> str:
        "Allow only HT01, HT02, HT03, HT04, HT05, HT06, HT07, HT08, HT09, HT10"
        if v not in [f"HT{i:02d}" for i in range(1, 11)] and v != "NA":
            raise ModelError(
                given=v, accepted="hashtag be one of HT01, HT02, HT03, HT04, HT05, HT06, HT07, HT08, HT09, HT10"
            )
        return v

    @validator("sample_tube", pre=True)
    def validate_sample_tube(cls, v: str) -> str:
        if re.fullmatch("T[1-9]", v):
            return v
        else:
            raise ModelError(given=v, accepted="sample_tube must start with T and be 1 digit long; i.e. T1")

    @validator("sort_pool", pre=True)
    def validate_sort_pool(cls, v: str) -> str:
        if re.fullmatch("(P[0-9]{2}|NA)", v):
            return v
        else:
            raise ModelError(given=v, accepted="sort_pool must start with P and be 2 digits long; i.e. P01")

    @validator("file_subset", pre=True)
    def validate_file_subset(cls, v: str) -> str:
        if re.fullmatch("[a-j]", v):
            return v
        else:
            raise ModelError(given=v, accepted="file_subset must be a single letter between 'a' and 'j'; i.e. a")


class ControlDataFilesModel(ControlSharedFileFieldsModel):
    sort_file_type: Literal["Primary"]
    extention: Literal[".fcs"]


class ControlPopulationSummaryFilesModel(ControlSharedFileFieldsModel):
    sort_file_type: Literal["Summary"]
    extention: Literal[".csv"]


class ControlScreenshotsModel(ControlSharedFileFieldsModel):
    sort_file_type: Literal["Capture"]
    extention: Literal[".png", ".PNG"]


class ControlSortReportsModel(ControlSharedFileFieldsModel):
    sort_file_type: Literal["SortRpt"]
    extention: Literal[".csv", ".pdf"]


class PackingDataFilesModel(PackingSharedFileFieldsModel):
    sort_file_type: Literal["Primary"]
    extention: Literal[".fcs"]


class PackingSortReportsModel(PackingSharedFileFieldsModel):
    sort_file_type: Literal["SortRpt"]
    extention: Literal[".csv", ".pdf"]


class PackingPopulationSummaryFilesModel(PackingSharedFileFieldsModel):
    sort_file_type: Literal["Summary"]
    extention: Literal[".csv"]


class PackingScreenshotsModel(PackingSharedFileFieldsModel):
    sort_file_type: Literal["Capture"]
    extention: Literal[".png", ".PNG"]


class PrescreenDataFilesModel(PrescreenSharedFileFieldsModel):
    sort_file_type: Literal["Primary"]
    extention: Literal[".fcs"]


class PrescreenPopulationSummaryFilesModel(PrescreenSharedFileFieldsModel):
    sort_file_type: Literal["Summary"]
    extention: Literal[".csv"]
    file_path: Optional[str]

    def get_series(self) -> pd.Series:
        return pd.Series(self.__dict__)

    @validator("file_path")
    def validate_file_path(cls, v: str) -> str:
        if not Path(v).exists():
            raise ModelError(given=v, accepted="file_path must exist")
        # population_csv = pd.read_csv(v, skiprows=6)  # pyright: reportUnknownMemberType=false
        # parsable_df = population_csv.iloc[:, :3]
        # population_list = [
        #     "All Events",
        #     "P1",
        #     "P2",
        #     "P3",
        #     "P4",
        #     "P5",
        #     "P6",
        #     "P7",
        #     "P11",
        #     "P12",
        #     "P13",
        #     "P14",
        #     "P15",
        #     "P16",
        #     "P17",
        #     "P18",
        #     "P19",
        #     "P8",
        #     "P9",
        #     "P10",
        # ]
        # parent_list = [
        #     "All Events",
        #     "P1",
        #     "P2",
        #     "P3",
        #     "P4",
        #     "P5",
        #     "P6",
        #     "P5",
        #     "P11",
        #     "P12",
        #     "P13",
        #     "P12",
        #     "P15",
        #     "P15",
        #     "P5",
        #     "P18",
        #     "P5",
        #     "P8",
        #     "P8",
        # ]
        # if parsable_df["Population"].tolist() != population_list:
        #     raise ModelError(
        #         given=v, accepted=f"file_path must contain the correct populations in order {population_list}"
        #     )
        # if parsable_df["Parent Name"].tolist()[1:] != parent_list:
        #     raise ModelError(given=v, accepted="file_path must contain the correct populations in order {parent_list}")
        return v


class PrescreenScreenshotsModel(PrescreenSharedFileFieldsModel):
    sort_file_type: Literal["Capture"]
    extention: Literal[".png", ".PNG"]


class ClinicalDataFilesModel(ClinicalSharedFileFieldsModel):
    """G00x -> Sorts -> Sort_RunDate..../ClinicalSamples/DataFilesFromDV/{Fields}.csv"""

    sort_file_type: Literal["Primary"]
    extention: Literal[".fcs"]


class ClinicalPopulationSummaryFilesModel(ClinicalSharedFileFieldsModel):
    """G00x -> Sorts -> Sort_RunDate..../ClinicalSamples/PopulationSummarFilesFromDV/{Fields}.csv"""

    sort_file_type: Literal["Summary"]
    extention: Literal[".csv"]
    file_path: Optional[str]

    def get_series(self) -> pd.Series:
        return pd.Series(self.__dict__)

    @validator("file_path")
    def validate_file_path(cls, v: str) -> str:
        if not Path(v).exists():
            raise ModelError(given=v, accepted="file_path must exist")
        # population_csv = pd.read_csv(v, skiprows=4)  # pyright: reportUnknownMemberType=false
        # parsable_df = population_csv.iloc[:, :3]
        # population_list = [
        #     "All Events",
        #     "P1",
        #     "P2",
        #     "P3",
        #     "P4",
        #     "P5",
        #     "P6",
        #     "P7",
        #     "P11",
        #     "P12",
        #     "P13",
        #     "P14",
        #     "P15",
        #     "P16",
        #     "P17",
        #     "P18",
        #     "P19",
        #     "P8",
        #     "P9",
        #     "P10",
        # ]
        # parent_list = [
        #     "All Events",
        #     "P1",
        #     "P2",
        #     "P3",
        #     "P4",
        #     "P5",
        #     "P6",
        #     "P5",
        #     "P11",
        #     "P12",
        #     "P13",
        #     "P12",
        #     "P15",
        #     "P5",
        #     "P18",
        #     "P5",
        #     "P8",
        #     "P8",
        # ]
        # if parsable_df["Population"].tolist() != population_list:
        #     raise ModelError(
        #         given=v,
        #         accepted=f"file_path must contain the correct populations in order {population_list}",
        #     )
        # if parsable_df["Parent Name"].tolist()[1:] != parent_list:
        #     raise ModelError(given=v, accepted=f"file_path must contain the correct populations in order {parent_list}")
        return v


class ClinicalScreenshotsModel(ClinicalSharedFileFieldsModel):
    """G00x -> Sorts -> Sort_RunDate..../ClinicalSamples/ScreenShotsFromDV/{Fields}.csv"""

    sort_file_type: Literal["Capture"]
    extention: Literal[".png", ".PNG"]


class ClinicalSortReportsModel(ClinicalSharedFileFieldsModel):
    """G00x -> Sorts -> Sort_RunDate..../ClinicalSamples/SortReportsFromDV/{Fields}.csv"""

    sort_file_type: Literal["SortRpt"]
    extention: Literal[".csv", ".pdf"]


class XMLModel(BaseModel):
    run_date: date = Field(..., description="YYMMDD")
    sorter_id_run_id: str
    experimenter_initials: str
    file_subset: str
    extension: Literal[".xml"]

    @validator("run_date", pre=True)
    def extract_date(cls, v: str):
        return datetime.strptime(v, "%y%m%d").date()

    @validator("sorter_id_run_id", pre=True)
    def validate_sorter_id_run_id(cls, v: str) -> str:
        if re.fullmatch("S[0-9][A-Z][0-9]{2}", v):
            return v
        else:
            raise ModelError(
                given=v, accepted="sorter_id_run_id must start with S and be 4 char/digits long; i.e. S6C01"
            )

    @validator("experimenter_initials", pre=True)
    def validate_experimenter_initials(cls, v: str) -> str:
        if re.fullmatch("[A-Z]{2}", v):
            return v
        else:
            raise ModelError(given=v, accepted="experimenter_initials must be 2 char long; i.e. MP")

    @validator("file_subset", pre=True)
    def validate_file_subset(cls, v: str) -> str:
        if re.fullmatch("[a-j]", v):
            return v
        else:
            raise ModelError(given=v, accepted="file_subset must be a single letter between 'a' and 'j'; i.e. a")


class ClinicalSamplesModel(BaseModel):
    name: Literal["ClinicalSamples"]
    xmlfile: Optional[XMLModel]
    datafiles: Optional[ClinicalDataFilesModel]
    population_summary_files: Optional[ClinicalPopulationSummaryFilesModel]
    screenshots: Optional[ClinicalScreenshotsModel]
    sortreports: Optional[ClinicalSortReportsModel]


class ControlSamplesModel(BaseModel):
    name: Literal["ControlSamples"]
    xmlfile: Optional[XMLModel]
    datafiles: Optional[ControlDataFilesModel]
    population_summary_files: Optional[ControlPopulationSummaryFilesModel]
    screenshots: Optional[ControlScreenshotsModel]
    sortreports: Optional[ControlSortReportsModel]


class FlowExtrasModel(BaseModel):
    """Including manifest, counts and flags .xlsx"""

    run_date: date = Field(..., description="YYMMDD")
    name: Literal["FlowManifest", "PBMCCounts", "Flags", "LFNACounts"]
    extension: Literal[".xlsx"]

    @validator("run_date", pre=True)
    def extract_date(cls, v: str) -> date:
        return datetime.strptime(v, "%y%m%d").date()


class PrescreenModel(BaseModel):
    name: str
    run_date: Optional[date] = None
    upload_date: Optional[date] = None
    flow_manifest: Optional[FlowExtrasModel] = None
    xmlfile: Optional[XMLModel] = None
    datafiles: Optional[PrescreenDataFilesModel] = None
    population_summary_files: Optional[PrescreenPopulationSummaryFilesModel] = None
    screenshots: Optional[PrescreenScreenshotsModel] = None

    @root_validator(pre=True)
    def extract_dates(cls, values: dict[str, Any]) -> dict[str, Any]:
        if not re.fullmatch(
            r"Prescreen_RunDate\d{2}([0][1-9]|[1][0-2])([0-3][0-9])_UploadDate\d{2}([0][1-9]|[1][0-2])([0-3][0-9])",
            values["name"],
        ):
            raise ModelError(
                given=values["name"],
                accepted="Prescreen name must be in the format Prescreen_RunDateYYMMDD_UploadDateYYMMDD",
            )
        try:
            run_date_string = values["name"].split("_")[1].split("RunDate")[1]
            values["run_date"] = datetime.strptime(run_date_string, "%y%m%d").date()
            upload_date_string = values["name"].split("_")[2].split("UploadDate")[1]
            values["upload_date"] = datetime.strptime(upload_date_string, "%y%m%d").date()
        except (IndexError, ValueError):
            raise ModelError(
                given=values["name"], accepted="Prescreens name must be in the format Prescreen_RunDateYYMMDD"
            )
        return values


class PrescreensModel(BaseModel):
    name: Literal["Prescreens"]
    prescreen: Optional[PrescreenModel]


class SortModel(BaseModel):
    name: str
    run_date: Optional[date] = None
    upload_date: Optional[date] = None
    flow_manifest: Optional[FlowExtrasModel] = None
    clinical_sample: Optional[ClinicalSamplesModel] = None

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
            values["run_date"] = datetime.strptime(run_date_string, "%y%m%d").date()
            upload_date_string = values["name"].split("_")[2].split("UploadDate")[1]
            values["upload_date"] = datetime.strptime(upload_date_string, "%y%m%d").date()
        except (IndexError, ValueError):
            raise ModelError(
                given=values["name"], accepted="Sort name must be in the format Sort_RunDateYYMMDD_UploadDateYYMMDD"
            )
        return values


class SortsModel(BaseModel):
    name: Literal["Sorts"]
    sort: SortModel


class G00XModel(BaseModel):
    name: str
    sorts: Optional[SortsModel] = None
    prescreens: Optional[PrescreensModel] = None

    @validator("name", pre=True)
    def extract_scheme(cls, v: str) -> str:
        if re.fullmatch("[gG]00[1-5x](_Scheme_Example)?", v):
            return v
        raise ModelError(given=v, accepted="Root folder name must named G001-G005")


class NIHBoxModel(BaseModel):
    name: Literal["NIH_Box"]
    g00x: Optional[G00XModel]


class PopulationSortFile(BaseModel):
    """represents a single population sort data file .csv"""

    data: ClinicalPopulationSummaryFilesModel | PrescreenPopulationSummaryFilesModel
    file_path: Path

    def get_series(self) -> pd.Series:
        s = self.data.get_series()
        s["file_path"] = self.file_path
        return s


class PreScreenPopulationSortFile(BaseModel):
    """represents a single population sort data file .csv"""

    data: PrescreenPopulationSummaryFilesModel | ClinicalPopulationSummaryFilesModel
    file_path: Path

    def get_series(self) -> pd.Series:
        s = self.data.get_series()
        s["file_path"] = self.file_path
        return s


class PopulationSortFiles(BaseModel):
    """Represents a collection of population sort data files .csv"""

    data: list[PopulationSortFile | PreScreenPopulationSortFile] = []

    def get_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame([x.get_series() for x in self.data])
