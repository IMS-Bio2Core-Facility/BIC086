# -*- coding: utf-8 -*-
"""Tests for the multithreading.biomart submodule."""

from pathlib import Path

import pytest

from scripts.multithreading.biomart import concurrent_biomart

from ..custom_tmp_file import CustomTempFile


def test_raises_index_error(tmpdir: Path) -> None:
    """It raises an IndexError when there are no transcripts."""
    with pytest.raises(IndexError), CustomTempFile("transcriptId") as file:
        concurrent_biomart(file.filename, str(tmpdir))
