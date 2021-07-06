# -*- coding: utf-8 -*-
"""Tests for the scripts.data_handling.request submodule.

These unit tests are designed to test that genes are passed to the query correctly.
They do not test the API,
as tests of data returned by realworld API queries are best left to integrations tests.
"""

import pandas as pd

from scripts.data_handling.request import lut_check

lut = pd.DataFrame.from_dict({"name": ["abc", "def"], "id": ["a1", "a2"]})


def test_return_none() -> None:
    """It returns None when the gene is not found."""
    result = lut_check("ghi", lut)
    assert result is None


def test_return_str() -> None:
    """It returns the id when the gene is found."""
    result = lut_check("abc", lut)
    assert result == "a1"
