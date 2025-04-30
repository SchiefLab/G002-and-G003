"""
Validate G00x file path and internal data.
"""
import logging
from glob import glob
from pathlib import Path
from typing import Any

from g00x.validations.models import g003 as models

# from pydantic import BaseModel


logger = logging.getLogger("G003 Validation for Folder Structure")


class ValidateG003:
    """Perform validation of G00X folder structure and file naming scheme."""

    def __init__(self) -> None:
        self.models = models
        self.data_stats: list[tuple[models.Fields, Path, Path]] = []

    def get_model(self, pydantic_model: Any, path: Path, naming_parts: list[str]) -> Any:
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

        # Presorts always suffix and removes the need for sample tube field
        if naming_parts[-2].lower() in ["control", "presort"]:
            naming_parts = naming_parts[:-2] + ["na"] + naming_parts[-2:]

        if len(naming_parts) < field_num:
            raise ValueError(f"file {path.stem} has missing parts in its name")
        elif len(naming_parts) > field_num:
            raise ValueError(f"In {path.parent} \n\tfile {path.stem} has extra parts in its name")

        try:
            return pydantic_model(**{k: v for k, v in zip(pydantic_model.__fields__, naming_parts)})
        except ValueError as e:
            raise ValueError(f"In {path.parent} \n\tfile {path.stem} has invalid parts in its name {e}")

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
                    self.models.G00X(name=g00x)
                case [g00x, "Presorts", *_]:  # ignore presorts for validation; manual view only
                    continue
                # G00X -> Sorts -> Sort
                case [g00x, "Sorts", sort]:
                    self.models.Sort(name=sort)
                # All static directorys will be indirectly validated by the below cases
                case [*_] if path.is_dir():
                    continue
                # G00X -> Sorts -> Sort -> ClinicalSamples -> DataFilesFromFlowJo -> .wsp
                case [g00x, "Sorts", sort, "ClinicalSamples", "DataFilesFromFlowJo", _]:
                    self.get_model(self.models.DataFilesFromFlowJo, path, naming_parts)
                # G00X -> Sorts -> Sort -> ClinicalSamples -> DataFiles -> .fcs
                case [g00x, "Sorts", sort, "ClinicalSamples", "DataFilesFromMelody", _]:
                    continue  # ingore this folder, its not important to validate.
                    self.get_model(self.models.DataFilesFromMelody, path, naming_parts)
                # G00X -> Sorts -> Sort -> ClinicalSamples -> DataStats -> .xlsx
                case [g00x, "Sorts", sort, "ClinicalSamples", "DataStats", _]:
                    model = self.get_model(self.models.DataStats, path, naming_parts)
                    self.data_stats.append((model, Path(path_str), Path(root_folder)))
                # G00X -> Sorts -> Sort -> ClinicalSamples -> ScreenshotsCounts -> .jpg
                case [g00x, "Sorts", sort, "ClinicalSamples", "ScreenshotsCounts", _]:
                    self.get_model(self.models.ScreenshotsCounts, path, naming_parts)
                # G00X -> Sorts -> Sort -> ClinicalSamples -> ScreenshotsCountsFullImage -> .jpg
                case [g00x, "Sorts", sort, "ClinicalSamples", "ScreenshotsCountsFullImage", _]:
                    self.get_model(self.models.ScreenshotsCounts, path, naming_parts)
                # G00X -> Sorts -> Sort -> ClinicalSamples -> ScreenshotsMelodyStats -> .jpg
                case [g00x, "Sorts", sort, "ClinicalSamples", "ScreenshotsMelodyStats", _]:
                    self.get_model(self.models.ScreenShotsMelodyStats, path, naming_parts)
                # G00X -> Sorts -> Sort -> ClinicalSamples -> SortReports -> .pdf
                case [g00x, "Sorts", sort, "ClinicalSamples", "SortReports", _]:
                    self.get_model(self.models.SortReports, path, naming_parts)
                # ignore readme or docs or fail
                case _:
                    if naming_parts[-1] in [".docx", ".pdf", ".doc", ".md", ".txt"]:
                        continue
                    raise ValueError(f"{path} does not match any known scheme pattern")


# This is the main function that will be called by the CLI
def validate_g003_sorting(folder: Path) -> ValidateG003:
    """Validate the folder structure of a G00X from Box.

    Parameters
    ----------
    folder : Path
        Path to the G00X folder
    """
    validate = ValidateG003()
    validate.validate_scheme(root_folder=str(folder))
    print("Schema validation passed \u2713")
    return validate
