import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from g00x.data import Data
from g00x.flow.flow import parse_flow_data
from g00x.validations.sequencing_validation import validate_sequencing


def merge_flow_and_sequencing(data: Data, flow_path: Path, sequencing_path: Path) -> pd.DataFrame:
    # validate seqencing path
    sequencing_manifest = validate_sequencing(sequencing_path)

    # validate flow path
    flow_manifest = parse_flow_data(data, flow_path)

    # merge the data
    unique_cols: list[str] = [
        "ptid",
        "group",
        "weeks",
        "visit_id",
        "probe_set",
        "sample_type",
        "run_date",
        "sort_pool",
        "hashtag",
    ]

    # only get a single flow manifest
    single_flow = flow_manifest.query("run_purpose=='Sort'").groupby(unique_cols).head(1)[unique_cols]

    outer_merge = single_flow.astype({"run_date": str}).merge(
        sequencing_manifest.astype({"sorted_date": str}),
        left_on=["run_date", "sort_pool"],
        right_on=["sorted_date", "pool_number"],
        how="outer",
    )

    # change POSIX object to string
    outer_merge["run_dir_path"] = outer_merge["run_dir_path"].astype(str).replace("nan", np.nan)
    if outer_merge["ptid"].isna().any():
        raise ValueError("There are missing ptids in the merge. This means flow data is missing")

    # if there are vdj run ids missing, warn the user that sequencing data is missing but proceed.
    if outer_merge["vdj_run_id"].isna().any():
        warnings.warn(
            "There are missing vdj_run_ids in the merge. This means the sequencing data is missing", UserWarning
        )

    return outer_merge
