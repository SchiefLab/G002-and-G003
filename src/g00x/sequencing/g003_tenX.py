import hashlib
import logging
import os
import shutil
import subprocess
import warnings
from pathlib import Path
from typing import Any

import pandas as pd

from g00x.data import Data
from g00x.tools.path import cd, pathing

logger = logging.getLogger("G00x")


def get_hash_digest(df: pd.DataFrame) -> str:
    digestable: list[Any] = df["Sample"].to_list()
    h: str = hashlib.md5("".join(digestable).encode("utf-8")).hexdigest()
    return h


def g003_run_demultiplex(
    data: Data,
    merged_dataframe: pd.DataFrame,
    out: Path,
    overwrite: bool = False,
) -> pd.DataFrame:
    logger.info("Begining Demultiplexing...")

    # Remove RAMOS entries, which are only used as control to confirm vdj recovery efficiency?
    # NOTE: may not be needed after G002
    merged_dataframe = merged_dataframe.query("ptid!='RAMOS'").reset_index(drop=True)
    if merged_dataframe["run_dir_path"].isna().any():
        warnings.warn("There are missing run_dir_paths in the merge", UserWarning)

        # I really hate this but need to compensate for missing data.
        merged_dataframe = merged_dataframe[~merged_dataframe["run_dir_path"].isna()].reset_index(drop=True)

    merged_dataframe["vdj_run_dir_path"] = (
        merged_dataframe["run_dir_path"] + "/" + merged_dataframe["run_id"].astype(str)
    )
    merged_dataframe["cso_run_dir_path"] = (
        merged_dataframe["run_dir_path"] + "/" + merged_dataframe["run_id"].astype(str)
    )
    all_run_dir_paths: list[Path] = sorted(
        map(
            pathing,  # type: ignore
            list(
                set(
                    list(merged_dataframe["vdj_run_dir_path"].unique())  # type: ignore
                    + list(merged_dataframe["cso_run_dir_path"].unique())  # type: ignore
                )
            ),
        )
    )

    demux_cmd = []
    for run_path in all_run_dir_paths:
        working_dir = out / run_path.parent.stem
        demultiplexed_dir = working_dir / "demultiplexed"

        if not run_path.parent.exists():
            raise ValueError(f"Run directory {run_path.parent} does not exist")

        if not working_dir.exists():
            logger.info(f"Creating working directory {working_dir}")
            working_dir.mkdir(parents=True)
        else:
            logger.info(f"Working directory {working_dir} already exists. Skipping")
        logger.info(f"Demultiplexing from {run_path} to {working_dir}")

        # make demultiplexed directory
        if not demultiplexed_dir.exists():
            logger.info(f"Creating demultiplexedc directory {demultiplexed_dir}")
            demultiplexed_dir.mkdir()
        else:
            logger.info(f"Demultiplex directory {demultiplexed_dir} already exists. Skipping...")

        vdj_run_id_dataframe = merged_dataframe.query(f"vdj_run_dir_path == '{run_path}'").copy()
        cso_run_id_dataframe = merged_dataframe.query(f"cso_run_dir_path == '{run_path}'").copy()

        # VDJ indexes
        vdj_run_id_dataframe["vdj-sample_name"] = "vdj-" + vdj_run_id_dataframe["vdj_index"]
        vdj_indexes = sorted(list(vdj_run_id_dataframe["vdj_index"].unique()))
        vdj_sample_names = ["vdj-" + i for i in vdj_indexes]
        vdj_csv = pd.DataFrame(
            {
                "Lane": ["*"] * len(vdj_indexes),
                "Sample": vdj_sample_names,
                "index": vdj_indexes,
            }
        )

        # CSO
        cso_run_id_dataframe["cso-sample_name"] = "cso-" + cso_run_id_dataframe["cso_index"]
        cso_indexes = sorted(list(cso_run_id_dataframe["cso_index"].unique()))
        cso_sample_names = ["cso-" + i for i in cso_indexes]
        cso_csv = pd.DataFrame(
            {
                "Lane": ["*"] * len(cso_indexes),
                "Sample": cso_sample_names,
                "index": cso_indexes,
            }
        )

        # write out sample sheet
        combined_csv = pd.concat([vdj_csv, cso_csv]).reset_index().drop("level_0", axis=1)
        hash_output = get_hash_digest(combined_csv)
        hash_running_dir = demultiplexed_dir / Path(hash_output)

        if overwrite and hash_running_dir.exists():
            logger.info(f"Overwrite is set to True. Removing {hash_running_dir}")
            shutil.rmtree(hash_running_dir)
        if not overwrite and hash_running_dir.exists():
            vdj_indexes_to_update = vdj_run_id_dataframe.index
            cso_indexes_to_update = cso_run_id_dataframe.index
            merged_dataframe.loc[vdj_indexes_to_update, "vdj_fastq_dir"] = str(
                hash_running_dir / Path("outs/fastq_path")
            )
            merged_dataframe.loc[vdj_indexes_to_update, "vdj_sample_name"] = (
                "vdj-" + merged_dataframe.loc[vdj_indexes_to_update, "vdj_index"]
            )
            merged_dataframe.loc[cso_indexes_to_update, "cso_fastq_dir"] = str(
                hash_running_dir / Path("outs/fastq_path")
            )
            merged_dataframe.loc[cso_indexes_to_update, "cso_sample_name"] = (
                "cso-" + merged_dataframe.loc[cso_indexes_to_update, "feature_index"]
            )
            continue

        csv_output = demultiplexed_dir / Path(f"sample_sheet.csv")
        logger.info(f"Writing sample sheet to {csv_output}")
        combined_csv.to_csv(csv_output, index=False)

        # where are we before we start
        logger.info(f"Saving dir {Path('.').absolute()}")
        current_working_dir = Path(".").absolute()

        # go to working dir
        logger.info(f"Changing working directory to {demultiplexed_dir}")
        os.chdir(demultiplexed_dir)

        try:
            demux_cmd = [
                data.get_cellranger_path(),
                "mkfastq",
                "--csv",
                csv_output,
                "--run",
                run_path,
                "--id",
                hash_output,
                "--localcores=48",
                "--uiport=40575",
                "--jobmode=local",
            ]
            command = subprocess.run(demux_cmd)
            if command.returncode != 0:
                raise ValueError(f"Demultiplexing failed with command {demux_cmd}")
            vdj_indexes_to_update = vdj_run_id_dataframe.index
            cso_indexes_to_update = cso_run_id_dataframe.index
            merged_dataframe.loc[vdj_indexes_to_update, "vdj_fastq_dir"] = str(
                hash_running_dir / Path("outs/fastq_path")
            )
            merged_dataframe.loc[vdj_indexes_to_update, "vdj_sample_name"] = (
                "vdj-" + merged_dataframe.loc[vdj_indexes_to_update, "vdj_index"]
            )
            merged_dataframe.loc[cso_indexes_to_update, "cso_fastq_dir"] = str(
                hash_running_dir / Path("outs/fastq_path")
            )
            merged_dataframe.loc[cso_indexes_to_update, "cso_sample_name"] = (
                "cso-" + merged_dataframe.loc[cso_indexes_to_update, "cso_index"]
            )
        except Exception as e:
            logger.error(f"Failed to run {demux_cmd}")
            raise e
        finally:
            os.chdir(current_working_dir)
        logger.info(f"Saving merged dataframe to {demultiplexed_dir}")

    return merged_dataframe


def ensure_singleton(df: pd.DataFrame, column: str) -> bool:
    if len(df[column].unique()) != 1:
        return False
    return True


def g003_run_vdj(
    data: Data,
    demux_dataframe: pd.DataFrame,
    out: Path,
    overwrite: bool = False,
) -> pd.DataFrame:
    """Run the VDJ 10x pipeline for G003

    Parameters
    ----------
    data : Data
        G00x data pathways
    demux_dataframe : pd.DataFrame
        The demux dataframe from the demultiplex pipelien
    out : Path
        output to save the dataframe
    Returns
    -------
    pd.DataFrame
        The updated dataframe with vdj paths
    Raises
    ------
    ValueError
        Non singletons in the dataframe
    """
    logger.info("Running VDJ")

    # ensure that vdj_fastq_dir is not null
    if demux_dataframe["vdj_fastq_dir"].isnull().any():
        logger.error("vdj_fastq_dir is null")
        raise ValueError("vdj_fastq_dir is null")
    # ensure that vdj_fastq_sample is not null
    if demux_dataframe["vdj_sample_name"].isnull().any():
        logger.error("vdj_sample_name is null")
        raise ValueError("vdj_sample_name is null")

    groupby = demux_dataframe.groupby(["vdj_fastq_dir", "vdj_sample_name"])
    enumerate_groupby = enumerate(groupby)
    for numerator, (index, group_df) in enumerate_groupby:
        # first get the fastq path which will be first argument of gropuby index
        fastq_path = index[0]

        # then get sample name which will be second argument of groupby index
        sample_name = index[1]

        if not ensure_singleton(group_df, "vdj_fastq_dir"):
            error = f"vdj_fastq_dir is not singleton {group_df['vdj_fastq_dir'].unique()}"
            logger.error(error)
            raise ValueError(error)

        if not ensure_singleton(group_df, "run_dir_path"):
            error = f"run_dir_path is not singleton {group_df['run_dir_path'].unique()}"
            logger.error(error)
            raise ValueError(error)

        # we can make a working dir in run000x/vdj
        print("here", out)
        working_dir = out / Path(group_df["run_dir_path"].unique()[0]).stem / Path("vdj")

        if working_dir.exists():
            logger.info(f"Skipping creating{working_dir} as it already exists")
        else:
            working_dir.mkdir(parents=True)
            logger.info(f"Creating {working_dir}")

        # log current direcotry
        current_dir = os.getcwd()
        logger.info(f"Current directory is {current_dir}")

        # the actual output will be in vdj_output_000N
        vdj_output = working_dir / f"vdj_output_{str(numerator).zfill(4)}"

        if vdj_output.exists() and overwrite:
            logger.info(f"Removing {vdj_output} as overwrite is set to True")
            shutil.rmtree(vdj_output)
        if vdj_output.exists() and not overwrite:
            logger.info(f"{vdj_output} already exists. Skipping and adding path to manifest.")
            demux_dataframe.loc[group_df.index, "vdj_output"] = str(vdj_output)
            continue
        try:
            logger.info(f"Changing directory to {working_dir}")
            os.chdir(working_dir)
            vdj_cmd = [
                data.get_cellranger_path(),
                "vdj",
                "--id",
                vdj_output.name,  # unique_name
                "--sample",
                sample_name,
                "--reference",
                data.get_vdj_path(),
                "--fastqs",
                fastq_path,
                "--localcores=48",
                "--uiport=40575",
                "--jobmode=local",
            ]
            command_string = " ".join(vdj_cmd)
            logger.info(f"Running {command_string}")
            f = subprocess.run(vdj_cmd)
            if f.returncode != 0:
                logger.error(f"Failed to run {command_string}")
                raise ValueError(f"Failed to run {command_string}")
            demux_dataframe.loc[group_df.index, "vdj_output"] = str(vdj_output)
        finally:
            logger.info(f"changing directory back to {current_dir}")
            os.chdir(current_dir)
    # if not Path(vdj_frame_output).parent.exists():
    #     Path(vdj_frame_output).parent.mkdir()
    #     logger.info(f"Created {Path(vdj_frame_output).parent}")

    # logger.info(f"Saving vdj dataframe to {Path(vdj_frame_output).stem}")
    # demux_dataframe.to_feather(f"{Path(vdj_frame_output).parent}/{Path(vdj_frame_output).stem}.feather")
    # demux_dataframe.to_csv(str(out) + ".csv", index=False)
    return demux_dataframe


def get_hashtaglookup(data: Data, hto: str) -> str:
    hashtag_lookup = data.get_g003_hto_gates()
    hto_seq: str = hashtag_lookup.loc[hto, "seq"]  # type: ignore
    return hto_seq


def get_single_entry(data: Data, row: pd.Series) -> pd.Series:
    return pd.Series(
        {
            "id": row["hto"],
            "name": row["cso_sample_name"] + "_" + row["hto"],
            "read": "R2",
            "pattern": "5PNNNNNNNNNN(BC)",
            "sequence": get_hashtaglookup(data, row["hto"]),
            "feature_type": "Antibody Capture",
        }
    )


def get_library_df(fastq_path: str, sample: str):
    return pd.DataFrame(
        [
            {
                "fastqs": fastq_path,
                "sample": sample,
                "library_type": "Antibody Capture",
            }
        ]
    )


def get_feature_dataframe(data: Data, group_df: pd.DataFrame) -> pd.DataFrame:
    assert len(group_df["hto"]) == len(group_df)
    dfs: list[pd.Series] = []
    for _, row in group_df.iterrows():
        df = get_single_entry(data, row)
        dfs.append(df)
    return pd.DataFrame(dfs)


def g003_run_cso(
    data: Data,
    demux_dataframe: pd.DataFrame,
    out: Path,
    genome_reference: Path,
    overwrite: bool = False,
) -> pd.DataFrame:
    """Run the feature barcode 10x pipeline. It should be run after vdj, but that is your call bro.

    Parameters
    ----------
    data : Data
        G00x data pathways
    demux_dataframe : pd.DataFrame
        The demux dataframe from the demultiplex pipelien
    out : Path
        output to save the dataframe
    Returns
    -------
    pd.DataFrame
        The updated dataframe with vdj paths
    Raises
    ------
    ValueError
        Non singletons in the dataframe
    """
    logger.info("Running CSO")
    # ensure that vdj_fastq_dir is not null
    if demux_dataframe["cso_fastq_dir"].isnull().any():
        raise ValueError("vdj_fastq_dir is null")
    # ensure that vdj_fastq_sample is not null
    if demux_dataframe["cso_sample_name"].isnull().any():
        raise ValueError("cso_sample_name is null")

    groupby = demux_dataframe.groupby(["cso_fastq_dir", "cso_sample_name"])
    enumerate_groupby = enumerate(groupby)
    index: tuple[str, str]
    for numerator, (index, group_df) in enumerate_groupby:
        # first get the fastq path which will be first argument of gropuby index
        fastq_path = Path(index[0])

        # then get sample name which will be second argument of groupby index
        sample_name = str(index[1])

        if not ensure_singleton(group_df, "cso_fastq_dir"):
            raise ValueError(f"cso_fastq_dir is not singleton {group_df['vdj_fastq_dir'].unique()}")

        if not ensure_singleton(group_df, "run_dir_path"):
            raise ValueError(f"run_dir_path is not singleton {group_df['run_dir_path'].unique()}")

        # we can make a working dir in run000x/cso
        working_dir = out / Path(group_df["run_dir_path"].unique()[0]).stem / Path("cso")

        if working_dir.exists():
            logger.info(f"Skipping creating{working_dir} as it already exists")
        else:
            working_dir.mkdir(parents=True)
            logger.info(f"Creating {working_dir}")

        # log current direcotry
        current_dir = os.getcwd()
        logger.info(f"Current directory is {current_dir}")

        # the actual output will be in vdj_output_000N
        cso_output = working_dir / f"cso_output_{str(numerator).zfill(4)}"
        if cso_output.exists() and overwrite:
            logger.info(f"Removing {cso_output} as overwrite is set")
            shutil.rmtree(cso_output)
        if cso_output.exists() and not overwrite:
            logger.info(f"{cso_output} already exists. Skipping and adding path to manifest.")
            demux_dataframe.loc[group_df.index, "cso_output"] = str(cso_output)
            continue
        try:
            logger.info(f"Changing directory to {working_dir}")
            os.chdir(working_dir)
            feature_df = get_feature_dataframe(data, group_df)
            library_df = get_library_df(str(fastq_path), sample_name)
            feature_csv_name = f"feature_frame_{str(numerator).zfill(4)}.csv"
            library_csv_name = f"library_df_{str(numerator).zfill(4)}.csv"
            feature_df.to_csv(feature_csv_name, index=False)
            library_df.to_csv(library_csv_name, index=False)
            cso_cmd: list[str] = [
                str(data.get_cellranger_path()),
                "count",
                "--id",
                cso_output.name,
                "--feature-ref",
                feature_csv_name,
                "--libraries",
                library_csv_name,
                "--transcriptome",
                str(genome_reference),
                "--localcores=48",
                "--uiport=40576",
                "--jobmode=local",
            ]
            command_string: str = " ".join(cso_cmd)
            logger.info(f"Running {command_string}")
            print(command_string)
            cso_process = subprocess.run(command_string, shell=True)
            if cso_process.returncode != 0:
                raise ValueError(f"CSO failed with return code {cso_process.returncode}")
            demux_dataframe.loc[group_df.index, "cso_output"] = str(cso_output)
        finally:
            logger.info(f"changing directory back to {current_dir}")
            os.chdir(current_dir)
    # logger.info(f"Saving cso dataframe to {Path(cso_frame_output).stem}")
    # demux_dataframe.to_feather(f"{Path(cso_frame_output).parent}/{Path(cso_frame_output).stem}.feather")
    return demux_dataframe
