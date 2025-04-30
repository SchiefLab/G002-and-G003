import logging

import pandas as pd

from g00x.data import Data

logger = logging.getLogger("report")


def get_frequency_df(data: Data, combined_flow_df: pd.DataFrame) -> pd.DataFrame:
    """Get the flow frequency dataframe part

    Parameters
    ----------
    data : Data
        The data object housing G00x data
    combined_flow_df : pd.DataFrame
        The flow dataframe which combines multiple sampling pools. Usually only 1
    """

    # These frequency measures are a json which we can define
    frequency_measures = data.get_frequency_measures()

    # These are the count measures
    count_measures = data.get_pbmc_gates()

    # These are the indexes we will keep in the flow dataframe
    sort_index: list[str] = [
        "run_purpose",
        "run_date",
        "pubID",
        "ptid",
        "group",
        "weeks",
        "visit_id",
        "probe_set",
        "sample_type",
        "value_type",
    ]

    # get rid of these columns
    v = list(combined_flow_df.columns.difference(sort_index))

    # don't want to remove value though
    v.remove("value")
    logger.info(f"Frequency dataframe will drop {v}")

    # Create a new dataframe object
    new_dfs: list[pd.DataFrame] = []

    # Do pre-sort and sort in seperate instances
    for q in ["PreS", "Sort"]:
        # first do for run purpose
        sub_df = combined_flow_df.query(f"run_purpose=='{q}'")
        # now do for each count, which will yield a dataframe
        for count in count_measures:
            # print(freq)
            gate = count["gate"]
            long_name = count["easy_name"]

            gate_count: pd.DataFrame = sub_df.query(f"gate=='{gate}'").set_index(sort_index)["value"].reset_index()

            # Add some extra columns
            gate_count["value_type"] = "count"
            gate_count["short_name"] = long_name
            gate_count["long_name"] = long_name

            # Add the gate measures, e.g. P8 etc
            gate_count["pbmc_gate_expression"] = gate

            new_dfs.append(gate_count)

        # now do for each frequency, which will yield a dataframe
        for freq in frequency_measures:
            # print(freq)
            num_gate = freq["numerator"]
            dom_gate = freq["denominator"]
            num: pd.Series = sub_df.query(f"gate=='{num_gate}'").set_index(sort_index)["value"]
            denom: pd.Series = sub_df.query(f"gate=='{dom_gate}'").set_index(sort_index)["value"]
            freq_df = ((num / denom) * 100).reset_index()

            # Add some extra columns
            freq_df["value_type"] = "frequency"
            freq_df["short_name"] = freq["axis_name"]
            freq_df["long_name"] = freq["long_name"]

            # Add the gate measures, e.g. P8/P5 etc
            freq_df["pbmc_gate_expression"] = f"{num_gate}/{dom_gate}"

            new_dfs.append(freq_df)
    freq_df = pd.concat(new_dfs).reset_index(drop=True)
    return freq_df


def calculate_frequency_dataframe(data: Data, flow_df: pd.DataFrame) -> pd.DataFrame:
    # flow_df["pubID"] = flow_df["ptid"].map(data.get_g002_pubids_lookup())
    flow_df["pubID"] = flow_df["ptid"]

    # only work with these columns
    combine_columns: list[str] = [
        "run_purpose",
        "run_date",
        "pubID",
        "sort_id",
        "ptid",
        "group",
        "weeks",
        "visit_id",
        "probe_set",
        "sample_type",
        "gate",
        "phenotype",
        "value_type",
        "branch",
        "easy_name",
        "value",
        "sort_pool",
    ]

    logger.info(
        f"Removing {list(flow_df.columns.difference(combine_columns))} from flow dataframe before frequency calculation\n"
    )

    # get minimal columns
    flow_df_min = flow_df[combine_columns]

    # here is what we will groupby to collapse sort pools into one value
    gb_columns: list[str] = [
        "run_purpose",
        "run_date",
        "sort_id",
        "ptid",
        "pubID",
        "group",
        "weeks",
        "visit_id",
        "probe_set",
        "sample_type",
        "gate",
        "phenotype",
        "value_type",
        "branch",
        "easy_name",
    ]
    logger.info(f"Grouping by {gb_columns} to combine {list(flow_df_min.columns.difference(gb_columns))}")

    def get_combined_apply(group: pd.DataFrame):
        """groupby_combine will put the sort pools into a list and sum up the values"""
        return pd.Series({"sort_pool": group["sort_pool"].to_list(), "value": group["value"].sum()})

    # combine the sort pools with fancy groupby statement
    # drop na must be false here since we have nan values in `branch` field...but just sometimes..
    combined_sort_pools = flow_df_min.groupby(gb_columns, dropna=False).apply(get_combined_apply).reset_index()

    # now can calculate frequency
    frequency_df = get_frequency_df(data, combined_sort_pools)

    # need to convert to string to combine with sequencing
    frequency_df["run_date"] = frequency_df["run_date"].astype(str)
    return frequency_df


def calculate_sequence_frequency_dataframe(sequence_df: pd.DataFrame) -> pd.DataFrame:
    """Calulate the frequency of vrc01 sequences in the sequence dataframe

    Parameters
    ----------
    sequence_df : pd.DataFrame
        The sequence dataframe from the pipeline

    Returns
    -------
    pd.DataFrame
        A dataframe with the frequency of vrc01 sequences with meta data
    """

    # we have to add this to add it to the flow
    sequence_df["run_purpose"] = "Sort"
    sequence_df["top_c_allele"] = sequence_df["c_call_heavy"].str.split(",").str.get(0)

    gb_columns: list[str] = [
        "run_purpose",
        "run_date",
        "pubID",
        "ptid",
        "group",
        "weeks",
        "visit_id",
        "probe_set",
        "sample_type",
    ]

    # go into global
    gb = sequence_df.groupby(gb_columns)

    def get_combined_groupby(list_of_isotypes: list[str], column: str) -> list[pd.Series]:
        # We will collect all the series in a empty array
        collect_series: list[pd.Series] = []

        # iterate through the groupby
        for g, g_df in gb:
            for normalize in [False, True]:
                for isotype in list_of_isotypes:
                    # Create an empty series from the indexing group g
                    series = pd.Series(dict(zip(gb_columns, g)))

                    # Get the break down of the vrc01 class by isotype
                    break_down = (
                        g_df.query(f"{column}=='{isotype}'")["is_vrc01_class"]
                        .value_counts(normalize=normalize)
                        .reindex([True, False])
                        .fillna(0)
                    )

                    # if we normalize, it will be a percentage
                    if normalize:
                        sn = f"percent_{isotype}_vrc01_class_sequences"
                        ln = f"Percent of {isotype} sequences that are VRC01-class"
                        break_down = break_down * 100

                    # else it will be a count
                    else:
                        sn = f"num_{isotype}_vrc01_class_sequences"
                        ln = f"Number of {isotype} sequences that are VRC01-class"

                    # Add short and long name
                    series["short_name"] = sn
                    series["long_name"] = ln

                    # True is VRC01 class, false is non-vrc01 class
                    series["value"] = break_down[True]

                    # other meta data we can merge
                    series["value_type"] = "sequence measure"
                    series["pbmc_gate_expression"] = ""
                    collect_series.append(series)
        return collect_series

    list_of_isotypes = ["IGHM", "IGHG", "IGHA", "IGHD"]
    list_of_alleles = sequence_df["top_c_allele"].unique().tolist()
    collect_series = get_combined_groupby(list_of_isotypes, "top_c_call")
    collect_series += get_combined_groupby(list_of_alleles, "top_c_allele")

    for g, g_df in gb:
        for normalize in [False, True]:
            #:Need one outside of isotype for all IGHD-
            series = pd.Series(dict(zip(gb_columns, g)))

            # Get the break down of the vrc01 class by isotype
            break_down = (
                g_df["is_vrc01_class"].value_counts(normalize=normalize).reindex([True, False]).fillna(0)
            )  # add zero here so we can distinguish between 0 and nan

            # if we normalize, it will be a percentage
            if normalize:
                sn = f"percent_igdneg_vrc01_class_sequences"
                ln = f"Percent of IGD- sequences that are VRC01-class"
                break_down = break_down * 100

            # else it will be a count
            else:
                sn = f"num_igdneg_vrc01_class_sequences"
                ln = f"Number of IGD- sequences that are VRC01-class"

            # Add short and long name
            series["short_name"] = sn
            series["long_name"] = ln

            # True is VRC01 class, false is non-vrc01 class
            series["value"] = break_down[True]

            # other meta data we can merge
            series["value_type"] = "sequence measure"
            series["pbmc_gate_expression"] = ""
            collect_series.append(series)

    remade = pd.DataFrame(collect_series).reset_index(drop=True)

    # add one final series with all the sequences counted
    logger.info("Adding number of sequences")
    number_of_seqs = (
        gb.size()
        .reset_index()
        .rename({0: "value"}, axis=1)
        .assign(
            long_name="Number of sequences",
            short_name="num_sequences",
            value_type="sequence measure",
            pbmc_gate_expression="",
        )
    )
    remade = pd.concat([remade, number_of_seqs]).reset_index(drop=True)
    return remade


def combine_seq_and_flow(
    data: Data, seq_dataframe: pd.DataFrame, flow_dataframe: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Combine the sequence and flow dataframes into one dataframe with the requisite columns"""

    # Get frequency dataframe
    frequency_df = calculate_frequency_dataframe(data, flow_dataframe)

    # Sequence dataframe
    sequence_df = calculate_sequence_frequency_dataframe(seq_dataframe)

    combined_df = pd.concat([frequency_df, sequence_df]).reset_index(drop=True)
    row_index: list[str] = [
        "run_purpose",
        "run_date",
        "pubID",
        "ptid",
        "group",
        "weeks",
        "visit_id",
        "probe_set",
        "sample_type",
    ]

    def add_combined_frequencies(
        df: pd.DataFrame, x_value: str, y_value: str, short_name: str, long_name: str
    ) -> pd.DataFrame:
        x = df.query(f"short_name=='{x_value}'").set_index(row_index)["value"]
        y = df.query(f"short_name=='{y_value}'").set_index(row_index)["value"]
        x = x / 100
        y = y / 100
        new_df = ((x * y) * 100).reset_index()
        new_df["short_name"] = short_name
        new_df["value_type"] = "frequency"
        new_df["long_name"] = long_name
        return new_df

    # Add percent VRC01 class among IgG
    combined_df = combined_df.append(
        add_combined_frequencies(
            combined_df,
            "percent_ep_among_igg",
            "percent_IGHG_vrc01_class_sequences",
            "percent_vrc01_among_igg",
            "Percent of VRC01-class sequences among IgG",
        )
    )  # type: ignore

    # Add percent VRC01 class among IgM
    combined_df = combined_df.append(
        add_combined_frequencies(
            combined_df,
            "percent_ep_among_igm",
            "percent_IGHM_vrc01_class_sequences",
            "percent_vrc01_among_igm",
            "Percent of VRC01-class sequences among IgM",
        )
    )

    # Add percent VRC01 class among IgA
    combined_df = combined_df.append(
        add_combined_frequencies(
            combined_df,
            "percent_ep_among_iga",
            "percent_IGHA_vrc01_class_sequences",
            "percent_vrc01_among_iga",
            "Percent of VRC01-class sequences among IgA",
        )
    )

    # Add percent VRC01 class among IgD-
    combined_df = combined_df.append(
        add_combined_frequencies(
            combined_df,
            "percent_ep_among_igd_neg",
            "percent_igdneg_vrc01_class_sequences",
            "percent_vrc01_among_igd_neg",
            "Percent of VRC01-class sequences among IgD-",
        )
    )

    combined_df["value_type"] = pd.Categorical(combined_df["value_type"], ["count", "sequence measure", "frequency"])
    combined_df = combined_df.sort_values(["value_type", "short_name"]).reset_index(drop=True)
    combined_pivot = combined_df.pivot(row_index, "short_name", "value").reset_index()  # type: ignore
    combined_pivot.columns.name = ""

    # With long names
    combined_pivot_long = combined_df.pivot(row_index, "long_name", "value").reset_index()  # type: ignore
    combined_pivot_long.columns.name = ""

    # if you are missing columns, add them to data.get_long_name_sort_order()
    logger.info(
        f"Missing columns won't be in final: {set(data.get_long_name_sort_order()) - set(combined_pivot_long.columns)}"
    )
    combined_pivot_long = combined_pivot_long[data.get_long_name_sort_order()]

    return combined_pivot, combined_pivot_long, combined_df
