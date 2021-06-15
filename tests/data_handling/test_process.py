# -*- coding: utf-8 -*-
"""Tests for the scripts.data_handling.process submodule.

These unit tests are designed to test that the results of the query
are handled correctly.
They assume the data returned by the API query is formatted correctly,
as tests of data returned by realworld API queries are best left to integrations tests.
As such,
representative data is included in the tests.data module.

Attributes
----------
MANE : pd.DataFrame
    A minimal MANE dataset
"""
from io import StringIO

import pandas as pd

from scripts.data_handling.process import merge_data

from ..custom_tmp_file import (
    BIOMART_CONTENTS,
    GTEX_CONTENTS,
    MANE_CONTENTS,
    CustomTempFile,
)

MANE: pd.DataFrame = pd.read_csv(StringIO(MANE_CONTENTS))


def test_returns_dataframe() -> None:
    """It returns a DataFrame."""
    results = merge_data(
        CustomTempFile(GTEX_CONTENTS).filename,
        CustomTempFile(BIOMART_CONTENTS).filename,
        MANE,
    )
    assert type(results) == pd.DataFrame


def test_results_columns() -> None:
    """Its columns are named correctly ."""
    columns = [
        "gencodeId",
        "geneSymbol",
        "tissueSiteDetailId",
        "transcriptId",
        "median",
        "unit",
        "datasetId",
        "refseq",
        "#NCBI_GeneID",
        "HGNC_ID",
        "name",
        "RefSeq_prot",
        "Ensembl_prot",
        "MANE_status",
        "GRCh38_chr",
        "chr_start",
        "chr_end",
        "chr_strand",
    ]
    results = merge_data(
        CustomTempFile(GTEX_CONTENTS).filename,
        CustomTempFile(BIOMART_CONTENTS).filename,
        MANE,
    )
    assert all(x in columns for x in results.columns), "Found an unexpected column."
    assert len(results.columns) == len(
        columns
    ), f"There should be {len(columns)} columns"


def test_sorted_results() -> None:
    """The results are sorted by median.

    As the ``merge_data`` function technically sorts on "MANE_status" as well,
    it would be ideal to test that sort, too.
    However, it it impossible to know in advance how many will have this status,
    so we cannot check count.
    Additionally, we cannot check the sort as most values are NaN,
    and knowing the correct order would require prior knowledge about the number
    of GTEx transcripts and the number with MANE status.
    """
    results = merge_data(
        CustomTempFile(GTEX_CONTENTS).filename,
        CustomTempFile(BIOMART_CONTENTS).filename,
        MANE,
    )
    assert results["median"].is_monotonic
