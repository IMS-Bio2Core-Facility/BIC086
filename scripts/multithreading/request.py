# -*- coding: utf-8 -*-
"""Functions for concurrent GTEx requests.

.. warning::

   Though every effort has been made to ensure thread safety,
   the concurrency is, as of yet, untested.

.. note::

   This module will likely undergo significant refactoring to generalise
   the concurrency pipeline and move all data handling code to a separate
   model.

Any API query is likely to be an I/O bound process,
particularly if there are many to make.
As this is a single step,
many in, many out process,
concurrency can be easily achieved with a thread local ``requests.session``
and mapping with ``concurrent.futures.ThreadPoolExecutor.map``.
The call to ``concurrent.futures.ThreadPoolExecutor.map`` is handled in the analysis
script.
"""
import logging
import threading
from io import StringIO

import pandas as pd
import requests

thread_local = threading.local()
logger = logging.getLogger(__name__)


def _get_session(region: str) -> requests.Session:
    """Instantiate a thread local session.

    The requests session is not thread safe, per `this thread`_.
    To circumvent this, we create a thread local session. This means each session
    will still make multiple requests but remain isolated to its calling thread.

    .. _this thread:
        https://github.com/psf/requests/issues/2766

    Parameters
    ----------
    region : str
        The GTEx region to query.

    Returns
    -------
    requests.Session
    """
    # session still worth it - re-used by each thread
    if not hasattr(thread_local, "session"):
        thread_local.session = requests.Session()
        thread_local.session.headers.update({"Accept": "text/html"})
        thread_local.session.params.update(  # type: ignore
            {
                "datasetId": "gtex_v8",
                "tissueSiteDetailId": region,
                "format": "tsv",
            }
        )
    return thread_local.session


def gtex_request(region: str, gene: str, output: str) -> None:
    """Make a thead-safe gtex request against medianTranscriptExpression.

    A thread local session is provided by a call to ``_get_session``.
    This allows the reuse of sessions, which, among other things,
    provides significant speed ups.

    Parameters
    ----------
    region : str
        The GTEx region to query.
    gene : str
        The ENSG to query.
    output : str
        Where to save the output file.

    Raises
    ------
    requests.HTTPError
        When the GET request returns an error
    Exception
        Any other errors
    """
    s = _get_session(region)
    response = s.get(
        "https://gtexportal.org/rest/v1/expression/medianTranscriptExpression",
        params={
            "gencodeId": gene,
        },
    )

    try:
        response.raise_for_status()
    except requests.HTTPError:
        logger.exception(
            f"An HTTP Error occurred on gene {gene}. A detailed report follows"
        )
        raise
    except Exception:
        logger.exception(
            f"A non-HTTP Error occurred on gene {gene}. A detailed report follows"
        )
        raise
    else:
        logger.info(
            f"GET request for {gene} successful! Details: \n\t{response.headers}"
        )
        data = pd.read_csv(StringIO(response.text), sep="\t")
        data = data.loc[data["median"] > 0, :].sort_values("median", ascending=False)
        data.loc[:, ["gencodeId", "transcriptId"]] = data.loc[
            :, ["gencodeId", "transcriptId"]
        ].apply(lambda x: x.str.split(".").str.get(0))
        data.to_csv(output, index=False)
