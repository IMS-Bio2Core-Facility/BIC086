# -*- coding: utf-8 -*-
"""Data handling for *request* step."""
import contextlib
import logging
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)


def lut_check(gene: str, lut: pd.DataFrame) -> Optional[str]:  # type: ignore
    """Check that a gene is found in the Gencode annotations.

    If the gene is found, then it is converted to its Ensembl ID.
    If it is not found, then None is returned to signal to downstream processing
    that the gene is not in Gencode.

    `contextlib.suppress`_ is used rather than the more verbose ``try/except``
    as we are only interested in returning ``none`` for this specific case.
    Any and all other cases should fail fast and hard.

    .. _contextlib.suppress:
       https://docs.python.org/3/library/contextlib.html

    Note
    ----
    It's likely that your gene is in Gencode even if it is not found.
    Common reasons (at least for me!) that a gene might not be found include
    spelling errors and name errors (ie. using NGN2 instead of NEUROG2).

    Parameters
    ----------
    gene : str
        The gene name to be queried
    lut : pd.DataFrame
        The dataframe containing the name-to-id conversion for the genes

    Returns
    -------
    Optional[str]

    Example
    -------
    >>> lut = pd.DataFrame.from_dict({"name": ["ASCL1"], "id": ["ENSG00000139352.3"]})
    >>> lut_check("ASCL1", lut)
    'ENSG00000139352.3'
    >>> lut_check("NotAGene", lut)

    """
    with contextlib.suppress(IndexError):
        return lut.loc[lut["name"] == gene, "id"].values[0]
