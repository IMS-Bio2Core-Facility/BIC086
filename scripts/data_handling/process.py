# -*- coding: utf-8 -*-
"""Data handling for *process* step.

Attributes
----------
Pathlike : TypeVar
    A custom type to unify string and Path
"""
import logging
from pathlib import Path
from typing import TypeVar

import pandas as pd

logger = logging.getLogger(__name__)


Pathlike = TypeVar("Pathlike", Path, str)


def merge_data(
    gtex_path: Pathlike, bm_path: Pathlike, mane: pd.DataFrame
) -> pd.DataFrame:
    """Merge the data from previous pipeline queries.

    Parameters
    ----------
    gtex_path : Pathlike
        Path to the file containing GTEx query data.
    bm_path : Pathlike
        Path to the file containing BioMart query data.
    mane : pd.DataFrame
        A DataFrame containing MANE annotations

    Returns
    -------
    pd.DataFrame
        The merged DataFrame, containing the GTEx query,
        the BioMart query, and the MANE annotations.

    """
    gtex = pd.read_csv(gtex_path, header=0, index_col=None)
    bm = pd.read_csv(bm_path, header=0, index_col=None)
    data = (
        gtex.merge(bm, on=["geneSymbol", "gencodeId", "transcriptId"], how="outer")
        .merge(
            mane,
            on=["geneSymbol", "gencodeId", "transcriptId", "refseq"],
            how="left",
        )
        .sort_values(["median", "MANE_status"])
    )
    return data


def write_data(data: pd.DataFrame, writer: pd.ExcelWriter) -> None:
    """Write a DataFrame to an Excel file.

    Note
    ----
    This function is best used within a ``with`` block, so that:

    #. The ``ExcelWriter`` is already open.
    #. It will be properly closed.

    Parameters
    ----------
    data : pd.DataFrame
        The DataFrame to be written.
    writer : pd.ExcelWriter
        An **open** pandas ExcelWriter.

    """
    gene = data["geneSymbol"].unique()[0]
    data.to_excel(writer, index=False, sheet_name=gene)
    logger.info(f"{gene} add to output file.")
