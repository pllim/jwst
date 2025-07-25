"""Test file utilities"""

import os

import pytest

from jwst.lib.file_utils import pushdir


def test_pushdir(tmp_path):
    """Test for successful change"""
    # Retrieve a temp folder
    current = os.getcwd()
    tmp_dir = str(tmp_path)

    # Temporarily move to the temp folder
    with pushdir(tmp_dir):
        assert tmp_dir == os.getcwd(), "New directory is not what was specified."
    assert current == os.getcwd(), "Directory was not restored"


def test_pusdir_fail():
    """Test for failing changing."""
    current = os.getcwd()
    with pytest.raises(Exception):
        with pushdir("Really_doesNOT-exist"):
            # Nothing should happen here. The assert should never be checked.
            assert False
    assert current == os.getcwd()
