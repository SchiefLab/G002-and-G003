import pandas as pd
from click.testing import CliRunner
from conftest import GeneralFixture
from pandas.testing import assert_frame_equal

from g00x.cli import g00x
from g00x.data import Data
from g00x.sequencing.merge import merge_flow_and_sequencing
from g00x.validations.sequencing_validation import validate_sequencing


def test_globus_endpoint_validator(fixture_setup: GeneralFixture):
    """Test the globus endpoint validator"""
    valid_box_structure = fixture_setup.get_valid_globus_endpoint_structure()
    validate_sequencing(valid_box_structure)


def test_sequence_merge(fixture_setup: GeneralFixture) -> None:
    """Test the sequence validator from the CLI"""
    click_runner = CliRunner()
    seq_path = fixture_setup.get_valid_globus_endpoint_structure()
    flow_path = fixture_setup.get_valid_box_data_structure()

    cmd = ["g002", "validate", "merge", "-s", str(seq_path), "-f", str(flow_path)]
    print("\nRunning:", " ".join(cmd))
    result = click_runner.invoke(g00x, cmd)
    if result.exit_code != 0:
        raise AssertionError(result.output)

    data = Data()
    output_file = merge_flow_and_sequencing(data, flow_path, seq_path).drop(columns=["run_dir_path"])
    finalized_file = pd.read_feather(fixture_setup.get_valid_merge_data()).drop(columns=["run_dir_path"])
    assert_frame_equal(finalized_file, output_file)
