"""Pytest conftest with all the fixture classes"""
from pathlib import Path

import pytest


class PreAuthFixtures:
    def __init__(self, tmp_path: Path, base_datadir: Path):
        self.base_datadir = base_datadir
        self.tmp_path = tmp_path


class G00XFixtures:
    def __init__(self, tmp_path: Path, base_datadir: Path):
        self.base_datadir = base_datadir
        self.tmp_path = tmp_path
        self.valid_box_dir = self.base_datadir / Path("sorting/G002/")

    def get_valid_box_data_structure(self) -> Path:
        """Get a contrived valid box structure"""
        return self.valid_box_dir


class GlobusFixtures:
    def __init__(self, tmp_path: Path, base_datadir: Path):
        self.base_datadir = base_datadir
        self.tmp_path = tmp_path
        self.globus_dir = self.base_datadir / Path("sequencing/G002")

    def get_valid_globus_endpoint_structure(self) -> Path:
        """Get a contrived valid box structure"""
        return self.globus_dir


class FlowFixtures:
    def __init__(self, tmp_path: Path, base_datadir: Path, fixture_dir: Path):
        self.base_datadir = base_datadir
        self.fixture_dir = fixture_dir
        self.tmp_path = tmp_path
        self.flow_data = self.fixture_dir / Path("flow_dataframe.feather")
        self.frequency_data = self.fixture_dir / Path("flow_dataframe_frequency.feather")
        self.merge_data = self.fixture_dir / Path("flow_dataframe_merged.feather")

    def get_valid_flow_data(self) -> Path:
        """Get a flow dataframe"""
        return self.flow_data

    def get_valid_freq_data(self) -> Path:
        """Get a flow dataframe"""
        return self.frequency_data

    def get_valid_merge_data(self) -> Path:
        """Get a merged flow and sequence dataframe"""
        return self.merge_data


class GeneralFixture(PreAuthFixtures, G00XFixtures, GlobusFixtures, FlowFixtures):
    def __init__(self, tmp_path_factory: pytest.TempPathFactory) -> None:
        tmp_path = tmp_path_factory.mktemp("g00x_fixture")
        base_datadir = Path("~/g002/G002").expanduser()
        fixture_dir = Path("tests/data/fixtures/")

        # then rest of attributes
        self.tmp_path = tmp_path
        self.base_datadir = base_datadir

        # three subclasses need to be initialized this way
        PreAuthFixtures.__init__(self, tmp_path, base_datadir)
        G00XFixtures.__init__(self, tmp_path, base_datadir)
        GlobusFixtures.__init__(self, tmp_path, base_datadir)
        FlowFixtures.__init__(self, tmp_path, base_datadir, fixture_dir)

    def get_card(self) -> Path:
        """get a png file path that has a card image. This is a nonsense file path to test unexpected files"""
        return self.base_datadir / "card.png"


@pytest.fixture(scope="session", autouse=True)
def fixture_setup(tmp_path_factory: pytest.TempPathFactory) -> GeneralFixture:
    return GeneralFixture(tmp_path_factory)
