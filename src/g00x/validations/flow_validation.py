"""
Validate G00x file path and internal data.
"""
import logging
from glob import glob
from pathlib import Path
from typing import Any

from g00x.validations.models.g00x import (
    ClinicalDataFilesModel,
    ClinicalPopulationSummaryFilesModel,
    ClinicalScreenshotsModel,
    ClinicalSortReportsModel,
    ControlDataFilesModel,
    ControlPopulationSummaryFilesModel,
    ControlScreenshotsModel,
    ControlSortReportsModel,
    FlowExtrasModel,
    G00XModel,
    PackingDataFilesModel,
    PackingPopulationSummaryFilesModel,
    PackingScreenshotsModel,
    PackingSortReportsModel,
    PopulationSortFile,
    PopulationSortFiles,
    PrescreenDataFilesModel,
    PrescreenModel,
    PreScreenPopulationSortFile,
    PrescreenPopulationSummaryFilesModel,
    PrescreenScreenshotsModel,
    SortModel,
    XMLModel,
)

logger = logging.getLogger("FlowValidation")


class ValidateG00X:
    def __init__(self) -> None:
        self.population_files: list[PopulationSortFile | PreScreenPopulationSortFile] = []
        self.prescreen_population_files: list[PreScreenPopulationSortFile | PopulationSortFile] = []

    def get_model(
        self, pydantic_model: Any, path: Path, naming_parts: list[str]
    ) -> ClinicalPopulationSummaryFilesModel | PrescreenPopulationSummaryFilesModel:
        """
        Populate model inferring variables by order of args

        Parameters
        ----------
        pydantic_model : BaseModel
            Pydantic model to populate
        naming_parts : list
            List of variables to populate model with. Must be in order of model definition.

        Returns
        -------
        BaseModel
            Pydantic model populated with variables
        """
        field_num = len([i for i in pydantic_model.__fields__.values() if i.required])

        if len(naming_parts) < field_num:
            raise ValueError(f"file {path.stem} has missing parts in its name")
        elif len(naming_parts) > field_num:
            raise ValueError(f"In {path.parent} \n\tfile {path.stem} has extra parts in its name")

        # these need a file path
        if (
            pydantic_model == PrescreenPopulationSummaryFilesModel
            or pydantic_model == ClinicalPopulationSummaryFilesModel
        ):  # isinstance doesn't work here
            add_dict = dict(zip(pydantic_model.__fields__, naming_parts))
            add_dict["file_path"] = str(path)
            return pydantic_model(**add_dict)
        return pydantic_model(**{k: v for k, v in zip(pydantic_model.__fields__, naming_parts)})

    def check_filename_len(self, path: Path, naming_parts: list[str]) -> None:
        """
        Check if file name and file type are the same length

        Parameters
        ----------
        path : str
            File path
        naming_parts : list
            Filename split by "_"
        """
        if len(naming_parts) < 14:
            raise ValueError(f"file {path.stem} has missing parts in its name")
        elif len(naming_parts) > 14:
            raise ValueError(f"In {path.parent} \n\tfile {path.stem} has extra parts in its name")

    def validate_scheme(
        self,
        root_folder: str,
    ) -> None:
        """
        Validate G00X folder structure scheme

        Parameters
        ----------
        root_folder : str
            Root folder to start validation
        """
        if not Path(root_folder).exists():
            raise ValueError(f"Folder {root_folder} does not exist")
        for path_str in glob(f"{root_folder}/**", recursive=True):
            path = Path(path_str)
            parts = path.relative_to(Path(root_folder).parent).parts
            naming_parts = path.stem.split("_") + [path.suffix]
            logger.debug(f"checking path: {parts}")
            logger.debug(f"naming parts:{ path}, {naming_parts}")
            match parts:
                # G00X
                case [g00x]:
                    G00XModel(name=g00x)
                # G00X -> Prescreens
                case [g00x, "Prescreens"]:
                    continue
                # G00X -> Prescreens -> Prescreen
                case [g00x, "Prescreens", prescreen]:
                    PrescreenModel(name=prescreen)
                # G00X -> Prescreens -> Prescreen -> .xml ; xmlFile only
                case [g00x, "Prescreens", prescreen, file] if file.endswith(".xml"):
                    self.get_model(XMLModel, path, naming_parts)
                # G00X -> Prescreens -> Prescreen -> .xlsx ; FlowManifest only
                case [g00x, "Prescreens", prescreen, file] if file.endswith(".xlsx"):
                    # could be flags, counts or  manifest
                    self.get_model(FlowExtrasModel, path, naming_parts)
                # G00X -> Prescreens -> Prescreen -> DataFiles
                case [g00x, "Prescreens", prescreen, "DataFilesFromDV"]:
                    continue
                # G00X -> Prescreens -> Prescreen -> PopulationSummaryFilesFromDV
                case [g00x, "Prescreens", prescreen, "PopulationSummaryFilesFromDV"]:
                    continue
                # G00X -> Prescreens -> Prescreen -> ScreenshotsFromDV
                case [g00x, "Prescreens", prescreen, "ScreenshotsFromDV"]:
                    continue
                # G00X -> Prescreens -> Prescreen -> DataFiles -> .fcs
                case [g00x, "Prescreens", prescreen, "DataFilesFromDV", _]:
                    self.get_model(PrescreenDataFilesModel, path, naming_parts)
                # G00X -> Prescreens -> Prescreen -> PopulationSummary -> .csv
                case [g00x, "Prescreens", prescreen, "PopulationSummaryFilesFromDV", _]:
                    model = self.get_model(PrescreenPopulationSummaryFilesModel, path, naming_parts)
                    self.prescreen_population_files.append(PreScreenPopulationSortFile(data=model, file_path=path))
                # G00X -> Prescreens -> Prescreen -> Screenshots -> .png
                case [g00x, "Prescreens", prescreen, "ScreenshotsFromDV", _]:
                    self.get_model(PrescreenScreenshotsModel, path, naming_parts)
                # G00X -> Sorts
                case [g00x, "Sorts"]:
                    continue
                # G00X -> Sorts -> Sort
                case [g00x, "Sorts", sort]:
                    SortModel(name=sort)
                # G00X -> Sorts -> Sort -> .xlsx
                case [g00x, "Sorts", sort, sort_file] if Path(sort_file).suffix == ".xlsx":
                    self.get_model(FlowExtrasModel, path, naming_parts)
                # G00X -> Sorts -> Sort -> ClinicalSamples
                case [g00x, "Sorts", sort, "ClinicalSamples"]:
                    continue
                # G00X -> Sorts -> Sort -> ClinicalSamples -> DataFiles
                case [g00x, "Sorts", sort, "ClinicalSamples", "DataFilesFromDV"]:
                    continue
                # G00X -> Sorts -> Sort -> ClinicalSamples -> PopulationSummaryFilesFromDV
                case [g00x, "Sorts", sort, "ClinicalSamples", "PopulationSummaryFilesFromDV"]:
                    continue
                # G00X -> Sorts -> Sort -> ClinicalSamples -> ScreenshotsFromDV
                case [g00x, "Sorts", sort, "ClinicalSamples", "ScreenshotsFromDV"]:
                    continue
                # G00X -> Sorts -> Sort -> ClinicalSamples -> SortReportsFromDV
                case [g00x, "Sorts", sort, "ClinicalSamples", "SortReportsFromDV"]:
                    continue
                # G00X -> Sorts -> Sort -> ClinicalSamples -> .xml
                case [g00x, "Sorts", sort, "ClinicalSamples", file] if file.endswith(".xml"):
                    self.get_model(XMLModel, path, naming_parts)
                # G00X -> Sorts -> Sort -> ClinicalSamples -> .xlsx (internal manifest)
                case [g00x, "Sorts", sort, "ClinicalSamples", file] if file.endswith(".xlsx"):
                    self.get_model(FlowExtrasModel, path, naming_parts)
                # G00X -> Sorts -> Sort -> ClinicalSamples -> DataFiles ->.fcs
                case [g00x, "Sorts", sort, "ClinicalSamples", "DataFilesFromDV", _]:
                    self.get_model(ClinicalDataFilesModel, path, naming_parts)
                # G00X -> Sorts -> Sort -> ClinicalSamples -> PopulationSummary -> .csv
                case [g00x, "Sorts", sort, "ClinicalSamples", "PopulationSummaryFilesFromDV", _]:
                    model = self.get_model(ClinicalPopulationSummaryFilesModel, path, naming_parts)
                    self.population_files.append(PopulationSortFile(data=model, file_path=path))
                # G00X -> Sorts -> Sort -> ClinicalSamples -> Screenshots -> .png
                case [g00x, "Sorts", sort, "ClinicalSamples", "ScreenshotsFromDV", _]:
                    self.get_model(ClinicalScreenshotsModel, path, naming_parts)
                # G00X -> Sorts -> Sort -> ClinicalSamples -> SortReports -> .csv | .pdf
                case [g00x, "Sorts", sort, "ClinicalSamples", "SortReportsFromDV", _]:
                    self.get_model(ClinicalSortReportsModel, path, naming_parts)
                # G00X -> Sorts -> Sort -> ControlSamples
                case [g00x, "Sorts", sort, "ControlSamples"]:
                    continue
                # G00X -> Sorts -> Sort -> ControlSamples -> DataFiles
                case [g00x, "Sorts", sort, "ControlSamples", "DataFilesFromDV"]:
                    continue
                # G00X -> Sorts -> Sort -> ControlSamples -> PopulationSummaryFilesFromDV
                case [g00x, "Sorts", sort, "ControlSamples", "PopulationSummaryFilesFromDV"]:
                    continue
                # G00X -> Sorts -> Sort -> ControlSamples -> ScreenshotsFromDV
                case [g00x, "Sorts", sort, "ControlSamples", "ScreenshotsFromDV"]:
                    continue
                # G00X -> Sorts -> Sort -> ControlSamples -> SortReportsFromDV
                case [g00x, "Sorts", sort, "ControlSamples", "SortReportsFromDV"]:
                    continue
                # G00X -> Sorts -> Sort -> ControlSamples -> .xml
                case [g00x, "Sorts", sort, "ControlSamples", file] if file.endswith(".xml"):
                    self.get_model(XMLModel, path, naming_parts)
                # G00X -> Sorts -> Sort -> ControlSamples -> .xlsx (internal manifest)
                case [g00x, "Sorts", sort, "ControlSamples", file] if file.endswith(".xlsx"):
                    self.get_model(FlowExtrasModel, path, naming_parts)
                # G00X -> Sorts -> Sort -> ControlSamples -> DataFiles ->.fcs
                case [g00x, "Sorts", sort, "ControlSamples", "DataFilesFromDV", _]:
                    self.get_model(ControlDataFilesModel, path, naming_parts)
                # G00X -> Sorts -> Sort -> ControlSamples -> PopulationSummary -> .csv
                case [g00x, "Sorts", sort, "ControlSamples", "PopulationSummaryFilesFromDV", _]:
                    model = self.get_model(ControlPopulationSummaryFilesModel, path, naming_parts)
                # G00X -> Sorts -> Sort -> ControlSamples -> Screenshots -> .png
                case [g00x, "Sorts", sort, "ControlSamples", "ScreenshotsFromDV", _]:
                    self.get_model(ControlScreenshotsModel, path, naming_parts)
                # G00X -> Sorts -> Sort -> ControlSamples -> SortReports -> .csv | .pdf
                case [g00x, "Sorts", sort, "ControlSamples", "SortReportsFromDV", _]:
                    self.get_model(ControlSortReportsModel, path, naming_parts)

                # G00X -> Sorts -> Sort -> PackingSamples
                case [g00x, "Sorts", sort, "PackingSamples"]:
                    continue
                # G00X -> Sorts -> Sort -> PackingSamples -> PopulationSummaryFilesFromDV
                case [g00x, "Sorts", sort, "PackingSamples", "PopulationSummaryFilesFromDV"]:
                    continue
                # G00X -> Sorts -> Sort -> PackingSamples -> ScreenshotsFromDV
                case [g00x, "Sorts", sort, "PackingSamples", "ScreenshotsFromDV"]:
                    continue
                # G00X -> Sorts -> Sort -> PackingSamples -> SortReportsFromDV
                case [g00x, "Sorts", sort, "PackingSamples", "SortReportsFromDV"]:
                    continue
                # G00X -> Sorts -> Sort -> PackingSamples -> DataFilesFromDV
                case [g00x, "Sorts", sort, "PackingSamples", "DataFilesFromDV"]:
                    continue

                # G00X -> Sorts -> Sort -> PackingSamples -> SortReports -> .csv | .pdf
                case [g00x, "Sorts", sort, "PackingSamples", "SortReportsFromDV", _]:
                    self.get_model(PackingSortReportsModel, path, naming_parts)

                case [g00x, "Sorts", sort, "PackingSamples", "PopulationSummaryFilesFromDV", _]:
                    model = self.get_model(PackingPopulationSummaryFilesModel, path, naming_parts)

                # G00X -> Sorts -> Sort -> PackingSamples -> .xml
                case [g00x, "Sorts", sort, "PackingSamples", file] if file.endswith(".xml"):
                    self.get_model(XMLModel, path, naming_parts)

                # G00X -> Sorts -> Sort -> PackingSamples -> ScreenshotsFromDV -> .png
                case [g00x, "Sorts", sort, "PackingSamples", "ScreenshotsFromDV", _]:
                    self.get_model(PackingScreenshotsModel, path, naming_parts)

                # ignore readme or docs or fail
                case [g00x, "Sorts", sort, "PackingSamples", "DataFilesFromDV", _]:
                    self.get_model(PackingDataFilesModel, path, naming_parts)
                case _:
                    if naming_parts[-1] in [".docx", ".pdf", ".doc", ".md", ".txt"]:
                        continue
                    raise ValueError(f"{path} does not match any known scheme pattern")

    def get_population_sort_files(self) -> PopulationSortFiles:
        """Get the population sort files as combined dataframe"""
        return PopulationSortFiles(data=self.population_files)

    def get_prescreen_population_sort_files(self) -> PopulationSortFiles:
        """Get the population sort file prescreens as a combined dataframe"""
        return PopulationSortFiles(data=self.prescreen_population_files)


# This is the main function that will be called by the CLI
def validate_g00x_box(folder: Path) -> ValidateG00X:
    """Validate the folder structure of a G00X from Box.

    Parameters
    ----------
    folder : Path
        Path to the G00X folder
    debug : bool
        The debug flag
    print_scheme : bool
        The print scheme flag
    """
    validate_g00x = ValidateG00X()
    validate_g00x.validate_scheme(root_folder=str(folder))
    print("Schema validation passed \u2713")
    return validate_g00x
