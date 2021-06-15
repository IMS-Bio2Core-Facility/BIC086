# -*- coding: utf-8 -*-
"""Tests for the multithreading.biomart submodule."""

from scripts.multithreading.biomart import concurrent_biomart
import pytest
from ..custom_tmp_file import CustomTempFile

def test_raises_index_error(tmpdir) -> None:
    """It raises an IndexError when there are no transcripts."""
    with pytest.raises(IndexError), CustomTempFile("transcriptId") as file:
        concurrent_biomart(file.filename, str(tmpdir))
