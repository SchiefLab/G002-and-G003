"""ABV - always be validating"""
from typing import Any

import pandas as pd

from g00x.validations.g003_flow_validation import ValidateG003


def pull_flow_from_validation(
    validation: ValidateG003,
    ptid2pubid: dict[str, str],
    ptid_prefix2group: dict[str, str],
    visit_id2week: dict[str, int],
) -> pd.DataFrame:
    """Pull flow data from validation object creating from parsing sort folder structure

    Parameters
    ----------
    validation : ValidateG003
        Validation object created from parsing sort folder structure
    ptid2group : dict
        Dictionary mapping ptid to group
    visit_id2week : dict
        Dictionary mapping visit_id to week
    Returns
    -------
    pd.DataFrame
        Dataframe containing flow data
    """
    data_stats = validation.data_stats

    records = []

    for model, file_path, root_folder in data_stats:
        # print(file_path)
        data_stats_df = pd.read_excel(file_path) if file_path.suffix == ".xlsx" else pd.read_csv(file_path)

        home = "/" + "/".join(root_folder.parts[1:-5])  # type: ignore
        relative_file_path = file_path.relative_to(home)  # type: ignore

        if len(model.visit_id) < 3:
            model.visit_id = "V0" + model.visit_id[-1]

        record_template = {
            "run_purpose": model.run_purpose,
            "run_date": model.run_date,
            "sort_id": model.sort_id,
            "ptid": model.ptid,
            "pubID": model.ptid,  # ptid2pubid[model.ptid] if not model.ptid.startswith('G001') else model.ptid,
            "group": ptid_prefix2group[model.ptid.split("-")[1]] if not model.ptid.startswith("G001") else 0,
            "weeks": visit_id2week[model.visit_id] if not model.ptid.startswith("G001") else model.visit_id,
            "visit_id": model.visit_id,
            "probe_set": model.probe_set,
            "sample_type": model.sample_type,
            "sort_software_dv": model.sort_software_dv,
            "sort_file_type": model.sort_file_type,
            "sample_tube": model.sample_tube,
            "gate": None,
            "phenotype": None,
            "value_type": "count",  # TODO: hardcoded?
            "extention": model.extention,
            "file_path": str(relative_file_path),  # type: ignore
            "file_subset": model.sort_pool_file_subset[-1],
            "value": None,
            "branch": None,
            "easy_name": None,
            "notes": None,  # TODO: no longer needed?
            "sort_pool": model.sort_pool_file_subset[:-1],
            "hashtag": None,  # TODO: no longer needed?
        }
        for _, row in data_stats_df.iterrows():
            record: dict[str, Any] = record_template.copy()

            for k, v in row.items():
                if isinstance(k, str):
                    row[k.strip()] = v

            if row.get("Gate_ID"):
                row["Gate Name"] = row["Gate_ID"]
            if row.get("Count"):
                row["Events Count"] = row["Count"]
            if row.get("Gate_Short_Name"):
                row["Gate Short Name"] = row["Gate_Short_Name"]

            if row["Gate Name"].strip().lower() == "all events":
                continue

            value = str(row["Events Count"]).strip().replace(",", "").replace("-", "")
            if value == "":
                value = None
            else:
                value = int(value)
            record.update(
                {
                    "gate": row["Gate Name"],
                    "branch": row["Gate Short Name"],
                    "phenotype": row[0],
                    "value": value,
                    "easy_name": row[0],
                }
            )
            records.append(record)

    flow_df = pd.DataFrame(records)

    group_key = [
        "run_purpose",
        # "run_date",
        # "sort_id",
        "ptid",
        "group",
        "weeks",
        "visit_id",
        "probe_set",
        "sample_type",
        "gate",
        # "sort_pool",
        # "hashtag",
    ]
    agg_fields = {
        "value": "sum",
        "run_date": lambda x: list(x),
        "sort_id": lambda x: list(x),
        "file_subset": lambda x: list(x),
        "file_path": lambda x: sorted(list(x)),
        "sort_pool": lambda x: list(x),
    }
    flow_df = (
        (
            flow_df.groupby(group_key)
            .agg(agg_fields | {k: "first" for k in flow_df.columns if k not in agg_fields})  # type: ignore
            .reset_index(drop=True)
        )[
            record_template  # type: ignore
        ]
        .explode(["run_date", "sort_id", "sort_pool"])
        .reset_index()
    )

    return flow_df
