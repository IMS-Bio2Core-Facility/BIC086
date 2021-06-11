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
GTEX_PATH : Path
    Path to representative GTEx data.
BM_PATH : Path
    Path to representative BioMart data.
MANE : pd.DataFrame
    Representative MANE data.
"""
from pathlib import Path

import pandas as pd

from scripts.data_handling.process import merge_data

GTEX_PATH = Path("tests", "data", "gtex_message.csv")
BM_PATH = Path("tests", "data", "biomart_message.csv")
MANE: pd.DataFrame = pd.read_csv(Path("tests", "data", "mane_minimal.csv"), index_col=0)


def test_returns_dataframe() -> None:
    """It returns a DataFrame."""
    results = merge_data(GTEX_PATH, BM_PATH, MANE)
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
    results = merge_data(GTEX_PATH, BM_PATH, MANE)
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
    results = merge_data(GTEX_PATH, BM_PATH, MANE)
    assert results["median"].is_monotonic
