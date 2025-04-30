"""
Validate G00x file path and internal data.
"""
import logging
import sys
from glob import glob
from pathlib import Path

import pandas as pd

from g00x.validations.models.g003_sequence import (
    IlluminaModel,
    RunModel,
    SequenceManifestModel,
)

sys.tracebacklimit = 1
logger = logging.getLogger("")


class ValidateSequencing:
    def __init__(self):
        self.sequencing_csvs: list[SequenceManifestModel] = []

    def validate_scheme(self, root_folder: str, debug: bool = False) -> None:
        """
        Validate G00X folder structure scheme

        Parameters
        ----------
        root_folder : str
            Root folder to start validation
        """
        for run_dir in glob(f"{root_folder}/*", recursive=False):
            run_dir_path = Path(run_dir)
            # check run000x first
            RunModel(name=run_dir_path.name)
            illumina_model_name: list[str] = []
            for run_dir_compoenents in sorted(glob(f"{root_folder}/{run_dir_path.name}/*")):
                run_dir_compoenents_path = Path(run_dir_compoenents)
                if run_dir_compoenents_path.name == "working_directory":
                    logger.info(f"Skipping working directory...in {run_dir_compoenents_path}")
                    continue
                elif run_dir_compoenents_path.is_dir():
                    model = IlluminaModel(run_dir=run_dir_compoenents_path)
                    illumina_model_name.append(model.run_dir.name)
                else:
                    self.sequencing_csvs.append(
                        SequenceManifestModel(
                            path=run_dir_compoenents_path,
                            illuimna_folder_name=illumina_model_name,
                        )
                    )

    def get_parsed_sampled_manifests(self) -> pd.DataFrame:
        """Get all parsed sample manifests"""
        return pd.concat([i.get_dataframe() for i in self.sequencing_csvs]).reset_index(drop=True)


# This is the main function that will be called by the CLI
def validate_g003_sequencing(folder: Path) -> pd.DataFrame:
    """Validate the folder structure of the globus_endpoint folder for Illumina seq data.

    Parameters
    ----------
    folder : Path
        Path to the G003 folder
    debug : bool
        The debug flag
    """
    validate_sequence = ValidateSequencing()
    validate_sequence.validate_scheme(root_folder=str(folder))
    logger.info("Sequence validation passed \u2713")
    return validate_sequence.get_parsed_sampled_manifests()
