import csv
import gzip
import logging
import re
from pathlib import Path
from typing import Any

import pandas as pd
import scipy
from sadie.airr import Airr
from sadie.airr.airrtable import AirrTable, LinkedAirrTable
from sadie.airr.methods import run_igl_assignment, run_mutational_analysis
from sadie.cluster import Cluster
from sadie.reference import Reference, References

from g00x.data import Data

logger = logging.getLogger("Airr")


def personalize(airr_df: pd.DataFrame, data: Data) -> pd.DataFrame:
    """
    Split the data, unused..precalculated
    """
    haplotype_lookup = data.get_personalized_vh12()

    # our baseline references is all available names
    baseline_human = References().get_dataframe().query("name=='human'")

    # take out VH12
    baseline_no_vh12 = baseline_human[baseline_human["gene"].str.split("*").str.get(0) != "IGHV1-2"].copy()

    references: References = References()
    haplotype_df = []

    # group by allele1 and 2
    for allele_group, allele_group_df in haplotype_lookup.groupby(["allele_1", "allele_2"]):
        # the haplotype name will be ex allele1_allele2, eg 02_02
        name = "_".join(list(map(lambda x: x.split("*")[-1], allele_group)))

        # only get the common
        need_alleles = list(set(list(allele_group)))

        # create a baseline ref with no vh12
        reference = Reference().from_dataframe(baseline_no_vh12.drop(["name", "imgt.sequence_aa"], axis=1))

        # add each allele one at a time
        for needed_allele in need_alleles:
            reference.add_gene({"species": "human", "gene": needed_allele, "source": "imgt"})
        print(f"adding_references:{name}")

        # add the reference to references with name
        references.add_reference(name=name, reference=reference)
        allele_group_df["haplotype"] = name
        haplotype_df.append(allele_group_df)

    haplotype_df = pd.concat(haplotype_df).reset_index(drop=True)
    personalized_df = []
    for haplo, haplot_df in haplotype_df.groupby("haplotype"):
        print("gettings", haplo)
        # from IPython import embed; embed()
        # only get ptids of that group
        pub_ids = haplot_df["ptid"].to_list()

        # sub df is unblind dataframe with just the subjects of interest
        sub_df = airr_df[airr_df["ptid"].isin(pub_ids)]
        if sub_df.empty:
            continue

        # make an api call to get the allele
        haplo = str(haplo)
        print(references, "references")
        airr_api = Airr(haplo, adaptable=True, references=references)
        heavy_airr_df = airr_api.run_dataframe(sub_df, "cellid", "sequence_heavy")
        light_airr_df = airr_api.run_dataframe(sub_df, "cellid", "sequence_light")
        # make some temporary dataframes to update
        h_ = heavy_airr_df.rename({"sequence_id": "cellid"}, axis=1).set_index("cellid")
        l_ = light_airr_df.rename({"sequence_id": "cellid"}, axis=1).set_index("cellid")

        # suffix with _heavy and light to update
        h_.columns = [str(i) + "_heavy" for i in h_.columns]
        l_.columns = [str(i) + "_light" for i in l_.columns]
        sub_df_cellid = sub_df.set_index("cellid")
        sub_df_cellid.update(h_)
        sub_df_cellid.update(l_)
        personalized_df.append(sub_df_cellid)

    # concat them all
    before_df_len = len(airr_df)
    working_dataframe = pd.concat(personalized_df).reset_index()
    if before_df_len != len(working_dataframe):
        raise ValueError(f"personalized {len(working_dataframe):,} != {before_df_len:,} before")
    logger.info(f"Personalized {len(working_dataframe):,} rows")
    return working_dataframe


def find_centroid(
    cluster_df: LinkedAirrTable | AirrTable,
    cluster_n: int = 5,
    lookup: list[str] | None = None,
) -> LinkedAirrTable | AirrTable:
    """Find the centroid of each cluster and mark it as such"""
    cluster_df["is_centroid"] = False
    if lookup is None:
        lookup = [
            "cdr1_aa_heavy",
            "cdr2_aa_heavy",
            "cdr3_aa_heavy",
            "cdr1_aa_light",
            "cdr2_aa_light",
            "cdr3_aa_light",
        ]
    for _, g_df in cluster_df.groupby(["cluster"]):
        if len(g_df) > 1:
            logger.info(f"Cluster df columns: {g_df.columns}")
            cluster_api = Cluster(
                LinkedAirrTable(g_df, key_column="cellid"),  # type: ignore
                linkage="average",
                lookup=lookup,
                pad_somatic=True,
                groupby=[
                    "v_call_top_heavy",
                    "v_call_top_light",
                    "junction_aa_length_heavy",
                    "top_c_call",
                    "is_vrc01_class",
                ],
            )
            _ = cluster_api.cluster(cluster_n)
            d_df = pd.DataFrame(cluster_api.distance_df)

            d_df.columns = g_df.index
            d_df.index = g_df.index

            centroid_index = d_df.mean().sort_values().head(1).index[0]
            cluster_df.loc[centroid_index, "is_centroid"] = True
        else:
            cluster_df.loc[g_df.index, "is_centroid"] = True
    return cluster_df


def cluster(
    dataframe: LinkedAirrTable, cluster_n: int, cluster_heavy_only: bool = False
) -> LinkedAirrTable | AirrTable:
    """Run a simple agglomorative clustering on the cdr columns that are pre-grouped by ptid, vrc01 class and isotype"""
    if cluster_heavy_only:
        lookup = [
            "cdr3_aa_heavy",
            "cdr3_aa_light",
        ]
    else:
        lookup = [
            "cdr1_aa_heavy",
            "cdr2_aa_heavy",
            "cdr3_aa_heavy",
            "cdr1_aa_light",
            "cdr2_aa_light",
            "cdr3_aa_light",
        ]
    cluster_api = Cluster(
        dataframe,
        linkage="average",  # TODO: change?
        lookup=lookup,
        # groupby=[
        #     "hcdr3_len",
        #     "lcdr3_len",
        #     "is_vrc01_class",
        #     "top_c_call",
        # ],
        pad_somatic=True,
        groupby=["v_call_top_heavy", "v_call_top_light", "junction_aa_length_heavy", "top_c_call", "is_vrc01_class"],
    )
    return find_centroid(cluster_api.cluster(cluster_n), cluster_n=cluster_n, lookup=lookup)


def find_100b(row: str) -> bool:
    "If we have a tryptophan at the purported 100b position, then return true"
    if not isinstance(row, str) or not row:
        return False
    elif len(row) < 6:
        return False
    elif row[-6] == "W":
        return True
    return False


def get_iter_kabat(x: str) -> list[str]:
    """return list of string up to kabat number 93"""
    return_list = []
    if isinstance(x, str):
        x = eval(x)
    for y in x:
        number = re.findall(r"\d+", y)[0]
        if int(number) > 93:
            continue
        return_list.append(y)
    return return_list


def add_mutational_sets(data: Data, dataframe: LinkedAirrTable):
    """Add which mutations are present in the heavy chain.

    There is cottrell_focused_v_common_score score
    which is the specific mutations we pay attention to in the VRC01 class


    Parameters
    ----------
    data : Data
        The data object
    dataframe : LinkedAirrTable
        A linked Airrtable
    """
    if not isinstance(dataframe, LinkedAirrTable):
        raise TypeError(f"{dataframe} must be type of LinkedAirrTable")
    if "is_vrc01_class" not in dataframe.columns:
        raise ValueError(f"is_vrc01_class field not found in {dataframe.columns}")
    if "mutations_heavy" not in dataframe.columns:
        raise ValueError(
            f"""mutations_heavy field not found in {dataframe.columns}.
            Run sadie.airr methods `run_mutational_analysis` to generate"""
        )
    vh12_reference_airr_table = data.get_vh12_reference_airr_table()
    cottrell_super_focus = data.get_cotrell_focus()
    cottrell_super_focus_positive = cottrell_super_focus["positive_set"]
    cottrell_super_focus_negative = cottrell_super_focus["negative_set"]

    cottrell_mabs = "N6 VRC27 VRC01 12A12 PCIN63_71I VRC-PG20".split()
    jardine_mabs = "12A12 3BNC60 VRC-PG04 VRC-PG20 VRC-CH31 VRC01".split()
    cotrell_focus = vh12_reference_airr_table[vh12_reference_airr_table["sequence_id"].isin(cottrell_mabs)]
    jardine_focus = vh12_reference_airr_table[vh12_reference_airr_table["sequence_id"].isin(jardine_mabs)]

    cotrell_focus_heavy_sets = set([item for sublist in cotrell_focus["mutations_heavy"].to_list() for item in sublist])
    cotrell_focus_light_sets = set([item for sublist in cotrell_focus["mutations_light"].to_list() for item in sublist])

    jardine_focus_heavy_sets = set([item for sublist in jardine_focus["mutations_heavy"].to_list() for item in sublist])
    jardine_focus_light_sets = set([item for sublist in jardine_focus["mutations_light"].to_list() for item in sublist])

    dataframe["cottrell_focused_v_common_heavy_positive"] = dataframe["mutations_heavy"].apply(
        lambda x: list(set(get_iter_kabat(x)).intersection(cottrell_super_focus_positive))
    )
    dataframe["cottrell_focused_v_common_heavy_negative"] = dataframe["mutations_heavy"].apply(
        lambda x: list(set(map(lambda y: y[1:-1], x)).intersection(cottrell_super_focus_negative))
    )
    dataframe["cottrell_focused_v_common_score"] = dataframe["cottrell_focused_v_common_heavy_positive"].apply(
        lambda x: len(x)
    ) - dataframe["cottrell_focused_v_common_heavy_negative"].apply(lambda x: len(x))

    dataframe["100bW"] = dataframe["junction_aa_heavy"].apply(find_100b)
    dataframe["cottrell_focused_v_common_score"] += dataframe["100bW"].apply(lambda x: {True: 1, False: 0}[x])

    dataframe["cottrell_v_common_heavy"] = dataframe["mutations_heavy"].apply(
        lambda x: list(set(get_iter_kabat(x)).intersection(cotrell_focus_heavy_sets))
    )
    dataframe["cottrell_v_common_heavy_score"] = dataframe["cottrell_v_common_heavy"].apply(lambda x: len(x))

    dataframe["cottrell_v_common_light"] = dataframe["mutations_light"].apply(
        lambda x: list(set(get_iter_kabat(x)).intersection(cotrell_focus_light_sets))
    )
    dataframe["cottrell_v_common_light_score"] = dataframe["cottrell_v_common_light"].apply(lambda x: len(x))

    dataframe["jardine_v_common_heavy"] = dataframe["mutations_heavy"].apply(
        lambda x: list(set(get_iter_kabat(x)).intersection(jardine_focus_heavy_sets))
    )
    dataframe["jardine_v_common_heavy_score"] = dataframe["jardine_v_common_heavy"].apply(lambda x: len(x))
    dataframe["jardine_v_common_light"] = dataframe["mutations_light"].apply(
        lambda x: list(set(get_iter_kabat(x)).intersection(jardine_focus_light_sets))
    )
    dataframe["jardine_v_common_light_score"] = dataframe["jardine_v_common_light"].apply(lambda x: len(x))
    return dataframe


def determine_if_vrc01(df: pd.DataFrame):
    # has vh12
    df["has_vh12_02_04"] = False
    df.loc[df[df["v_call_top_heavy"].str.contains(r"IGHV1-2\*0[24]")].index, "has_vh12_02_04"] = True

    # has five leng
    df["has_5_len"] = False
    df.loc[df[df["cdr3_aa_light"].str.len() == 5].index, "has_5_len"] = True

    # has both
    df["is_vrc01_class"] = False
    df.loc[df.query("has_vh12_02_04").query("has_5_len").index, "is_vrc01_class"] = True
    return df


def get_keyed_cso_file(grouped_df: pd.DataFrame, drop_lt: int = 100, drop_lt_percentile: float = 0.95) -> pd.DataFrame:
    """Get the CSO file keyed by the HTO

    Parameters
    ----------
    grouped_df : pd.DataFrame
        We are grouping by ptid and timepoint and looking at the cso output
    drop_lt : int, optional
        If a timepoint has less than this number of counts, we will drop it, by default 100
    drop_lt_percentile : float, optional
        If a sample does not have at least this many percentiles, we will drop it because its a doublet

    Returns
    -------
    pd.DataFrame
        The dataframe with the HTO as the index and the barcodes as the columns

    """
    cso_output: pd.Series = grouped_df["cso_output"]
    if len(cso_output.unique()) > 1:
        raise ValueError("cso_output has multiple values")
    else:
        cso_output_path: str = cso_output.iloc[0]

    matrix_path = Path(cso_output_path) / Path("outs/filtered_feature_bc_matrix/matrix.mtx.gz")
    logger.info(f"Reading in matrix from {matrix_path}")

    if not matrix_path.exists():
        raise ValueError(f"{matrix_path} does not exist")

    matrix = scipy.io.mmread(matrix_path)

    # now get the features and barcodes
    features_path = Path(cso_output_path) / Path("outs/filtered_feature_bc_matrix/features.tsv.gz")
    if not features_path.exists():
        raise ValueError(f"{features_path} does not exist")

    feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]

    barcodes_path = features_path = Path(cso_output_path) / Path("outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
    if not barcodes_path.exists():
        raise ValueError(f"{barcodes_path} does not exist")

    barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode="rt"), delimiter="\t")]

    # now combine the matrix with the barcodes and features
    matrix = pd.DataFrame.sparse.from_spmatrix(matrix)
    matrix.columns = barcodes
    matrix.columns = barcodes

    matrix.insert(loc=0, column="feature_id", value=feature_ids)
    matrix = matrix.set_index("feature_id")

    # just get the count matrix
    count_matrix = matrix.transpose()

    # we will throw out anything less than
    indexes_lt = count_matrix[(count_matrix > drop_lt).sum(axis=1) == 0].index
    logger.info(f"Dropping {len(indexes_lt)} the following indexes because they have less than {drop_lt} counts")
    logger.info(f"{indexes_lt}")

    # normal matrix
    nomral_matrix = count_matrix.apply(lambda x: (x / x.sum()), axis=1)

    indexes_no_dominant = nomral_matrix[(nomral_matrix > drop_lt_percentile).sum(axis=1) == 0].index
    logger.info(f"Dropping {len(indexes_no_dominant)} the following indexes because they have no dominant HTO")
    logger.info(f"{indexes_no_dominant}")

    # keep indexes to drop
    indexes_to_drop = indexes_lt.append(indexes_no_dominant)

    # anything that survived
    survival = nomral_matrix.drop(indexes_to_drop)

    def find_true(row: pd.Series) -> Any:
        # find the column thats true
        _index = [i for i, j in enumerate(row) if j][0]
        return row.index[_index]

    if survival.empty:
        logger.warn("Survival is empty, returning empty dataframe")
        return pd.DataFrame({"HTO": []})
    keyed_df = (survival > drop_lt_percentile).apply(find_true, axis=1).to_frame().rename({0: "HTO"}, axis=1)
    logger.info(f"Kept {len(keyed_df)} indexes")
    return keyed_df


def get_pairing(df: pd.DataFrame) -> pd.DataFrame:
    """Get the pairing of the heavy and light chain

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe with the airr data

    Returns
    -------
    pd.DataFrame
        A paired dataframe
    """
    df.insert(0, "cellhash", df["sequence_id"].str.split("_").str.get(0))
    complete_productive = df.query("productive and complete_vdj").reset_index(drop=True)
    complete_productive["locus"] = pd.Categorical(
        complete_productive["locus"], categories=["IGH", "IGK", "IGL"], ordered=True
    )

    # have to change this to a dataframe
    complete_productive = pd.DataFrame(complete_productive)

    # 1-1 pariring hashes
    good_hashes = (
        complete_productive.sort_values("locus")
        .groupby("cellhash")
        .apply(lambda x: "_".join(x["locus"].to_list()))
        .to_frame()
        .rename({0: "pair_type"}, axis=1)
        .reset_index()
        .query("pair_type=='IGH_IGK' or pair_type=='IGH_IGL'")["cellhash"]
        .to_list()
    )
    complete_pairing_candidates = complete_productive[complete_productive["cellhash"].isin(good_hashes)]
    complete_productive_pairing_candidates_heavy = complete_pairing_candidates.query("locus=='IGH'")
    complete_productive_pairing_candidates_light = complete_pairing_candidates.query("locus!='IGH'")
    paired = complete_productive_pairing_candidates_heavy.merge(
        complete_productive_pairing_candidates_light,
        on="cellhash",
        how="inner",
        suffixes=["_heavy", "_light"],  # type: ignore
    )
    return paired


def run_airr(
    data: Data,
    vdj_dataframe: pd.DataFrame,
    cso_dataframe: pd.DataFrame,
    output: Path | str,
    overwrite: bool,
    skip_mutation: bool,
    cluster_n: int = 5,
    cluster_heavy_only: bool = False,
) -> pd.DataFrame:
    """Run AIRR on the vdj files and demultiplex them with the CSO files"""
    logger.info("Running AIRR")
    difference = vdj_dataframe.columns.symmetric_difference(cso_dataframe.columns)
    logger.info(f"Columns in vdj but not cso: {difference}")

    # columns we want to merge the CSO on
    mergable_columns = [
        "ptid",
        "group",
        "weeks",
        "visit_id",
        "probe_set",
        "sample_type",
        "sort_pool",
        "hashtag",
        "sorted_date",
        "run_date",
        "vdj_index",
    ]

    if (vdj_dataframe.groupby(mergable_columns).size() > 1).any():
        raise ValueError("vdj has multiple rows for the same sample")

    combined_df = vdj_dataframe.merge(
        cso_dataframe[mergable_columns + ["feature_index", "cso_output"]],
        on=mergable_columns,
    )

    logger.info("Removing LNFA from analysis")
    combined_df = combined_df.query("sample_type=='PBMC'").reset_index(drop=True)

    airr_api = Airr("human", adaptable=True)
    complete_df = []
    for g, g_df in combined_df.groupby("vdj_output"):
        # declare these files to see if they exist

        # lets use filtered even though it should not matter since our pairing algorighm essentially gets the same thing
        contig_path = Path(str(g)) / Path("outs/filtered_contig.fasta")
        # contig_path = Path(str(g)) / Path("outs/all_contig.fasta")
        airr_out = Path(str(g)) / Path("outs/sadie_airr.feather")
        paired_airr_out = Path(str(g)) / Path("outs/paired_sadie_airr.feather")
        logging.info(f"VDJ path: {contig_path}")
        if airr_out.exists():
            if overwrite:
                logger.info(f"Overwriting {airr_out}\n")
                airr_file = airr_api.run_fasta(contig_path)
                airr_file.to_feather(airr_out)
            else:
                logger.info(f"{airr_out} exists\n")
                airr_file = pd.read_feather(airr_out)

        else:
            airr_file = airr_api.run_fasta(contig_path)
            airr_file.to_feather(airr_out)

        if paired_airr_out.exists():
            if overwrite:
                logger.info(f"pairing file {paired_airr_out} exists but overwrite was passed")
                paired_airr_file = get_pairing(airr_file)
                paired_airr_file.to_feather(paired_airr_out)
            else:
                logger.info(f"Skipping pairing because {paired_airr_out} exists\n")
                paired_airr_file = pd.read_feather(paired_airr_out)
        else:
            paired_airr_file = get_pairing(airr_file)
            if paired_airr_file.empty:
                paired_airr_file.reset_index().to_feather(paired_airr_out)
            paired_airr_file.to_feather(paired_airr_out)

        # locate those in dataframe
        combined_df.loc[g_df.index, "sadie_airr_path"] = str(airr_out)
        combined_df.loc[g_df.index, "paired_sadie_airr_path"] = str(paired_airr_out)
        g_df["sadie_airr_path"] = str(airr_out)
        g_df["paired_sadie_airr_path"] = str(paired_airr_out)

        # get the hash lookup
        keyed_file = get_keyed_cso_file(g_df)

        # join the airrtable with the with the CSO table
        with_cso = paired_airr_file.set_index("cellhash").join(keyed_file)

        # drop the ones without CSO
        no_cso = with_cso[with_cso["HTO"].isna()].index
        logger.info(f"Found {len(no_cso)} without CSO")
        with_cso = with_cso.drop(no_cso).reset_index()

        # merge the with the meta airr file with the hashtag (not cellhash) and HTO
        final_df = g_df.merge(with_cso, left_on="hashtag", right_on="HTO")
        complete_df.append(final_df)

    airr_df = pd.concat(complete_df).reset_index(drop=True)

    lookup_maps = data.get_g002_pubids_lookup()

    # Insert the pubid at 0 column for pubid identification
    # airr_df.insert(0, "pubID", airr_df["ptid"].map(lookup_maps))
    airr_df['pubID'] = airr_df['ptid']

    # Insert a unique cellid
    cellid_cols = ["pubID", "group", "weeks", "probe_set", "sort_pool", "cellhash"]
    airr_df.insert(
        0,
        "cellid",
        airr_df[cellid_cols].apply(lambda x: "_".join([str(i) for i in x]), axis=1),
    )
    if airr_df["cellid"].duplicated().any():
        raise ValueError(f"cellid is not unique {airr_df[airr_df['cellid'].duplicated()]['cellid']}")

    logger.info("Personalizing....")
    airr_df = personalize(airr_df, data)

    logger.info("Converting table to linked airrtable")
    airr_df_lat = LinkedAirrTable(airr_df, key_column="cellid")  # type: ignore

    logger.info("Running personalized VH12 analysis")

    # turn it back into a dataframe
    airr_df_lat = LinkedAirrTable(airr_df_lat, key_column="cellid")  # type: ignore

    logger.info("Running mutational analysis")
    airr_df_lat = run_mutational_analysis(airr_df_lat, "kabat")

    logger.info("Running iGL assignment")
    did_not_run: bool = False
    try:
        airr_df_lat = run_igl_assignment(airr_df_lat)
    except ValueError:
        did_not_run = True
        logger.error("IGL DID NOT RUN, MUST RUN MANUALLY, continuing")

    logger.info("Finding VRC01 Class")
    airr_df_lat = determine_if_vrc01(airr_df_lat)

    logger.info("Finding Mutational Sets")
    airr_df_lat = LinkedAirrTable(airr_df_lat, key_column="cellid")  # type: ignore
    airr_df_lat = add_mutational_sets(data, airr_df_lat)

    logger.info("Adding hcdr3_len, lcdr3_len, and top_c_call")
    airr_df_lat["hcdr3_len"] = airr_df_lat["cdr3_aa_heavy"].str.len()
    airr_df_lat["lcdr3_len"] = airr_df_lat["cdr3_aa_light"].str.len()
    airr_df_lat["top_c_call"] = airr_df_lat["c_call_heavy"].str[0:4].fillna("")

    logger.info(f"Clustering {len(airr_df_lat)} rows")
    logger.info(f"Airr df columns: {airr_df_lat.columns}")
    airr_df_lat = cluster(airr_df_lat, cluster_n, cluster_heavy_only)

    # these are the columns we will take from LinkedAirrTable
    mergable_cols: list[str] = [
        "cellid",
        # "iGL_aa_heavy",
        # "iGL_aa_light",
        "mutations_heavy",
        "mutations_light",
        "100bW",
        "is_vrc01_class",
        "cottrell_focused_v_common_heavy_positive",
        "cottrell_focused_v_common_heavy_negative",
        "cottrell_focused_v_common_score",
        "hcdr3_len",
        "lcdr3_len",
        "top_c_call",
        "cluster",
        "is_centroid",
    ]

    logger.info("Merging AIRR file with mutational analysis, iGL and mutational assignment")
    # merge back with LinkedAirrTable with subselected columns. This preserves our meta data
    airr_df = airr_df.merge(airr_df_lat[mergable_cols], on="cellid")
    logger.info(f"Saving AIRR file to {str(output)}.feather/csv.gz")

    # write out save in function so we can use it as an API call
    airr_df.to_feather(str(output) + ".feather")
    airr_df.to_csv(str(output) + ".csv.gz")
    return airr_df_lat
