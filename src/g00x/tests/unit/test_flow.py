import pandas as pd
from click.testing import CliRunner
from conftest import GeneralFixture
from pandas.testing import assert_frame_equal

from g00x.cli import g00x
from g00x.data import Data
from g00x.flow.flow import parse_flow_data

# from g00x.flow.frequency import get_frequency_df


def test_flow_workflow(fixture_setup: GeneralFixture) -> None:
    valid_box_structure = fixture_setup.get_valid_box_data_structure()
    data = Data()
    output_dataframe = parse_flow_data(data, valid_box_structure)
    output_dataframe["file_subset"] = output_dataframe["file_subset"].apply(lambda x: "_".join(x)).astype(str)
    output_dataframe["file_path"] = output_dataframe["file_path"].apply(lambda x: "_".join(x)).astype(str)
    output_dataframe = output_dataframe.sort_values(
        ["run_purpose", "ptid", "run_date", "group", "visit_id", "probe_set", "file_subset"]
    ).reset_index(drop=True)

    valid_dataframe = pd.read_feather(fixture_setup.get_valid_flow_data())
    valid_dataframe["file_subset"] = valid_dataframe["file_subset"].apply(lambda x: "_".join(x)).astype(str)
    valid_dataframe = valid_dataframe.sort_values(
        ["run_purpose", "ptid", "run_date", "group", "visit_id", "probe_set", "file_subset"]
    ).reset_index(drop=True)

    # need to drop file_subset because it is sorted differently in the CI.
    assert_frame_equal(
        output_dataframe.drop(["file_path", "file_subset"], axis=1),
        valid_dataframe.drop(["file_path", "file_subset"], axis=1),
    )

    # won't do this until analysis now
    # # test frequency
    # output_freq: pd.DataFrame = (
    #     get_frequency_df(data, output_dataframe)
    #     .sort_values(["run_purpose", "ptid", "run_date", "group", "visit_id", "probe_set"])
    #     .reset_index(drop=True)
    # )
    # valide_freq_dataframe: pd.DataFrame = (
    #     pd.read_feather(fixture_setup.get_valid_freq_data())
    #     .sort_values(["run_purpose", "ptid", "run_date", "group", "visit_id", "probe_set"])
    #     .reset_index(drop=True)
    # )
    # assert_frame_equal(output_freq, valide_freq_dataframe)


def test_flow_cli(fixture_setup: GeneralFixture) -> None:
    click_runner = CliRunner()
    path = fixture_setup.get_valid_box_data_structure()
    result = click_runner.invoke(g00x, ["g002", "validate", "flow", str(path)])
    assert result.exit_code == 0
