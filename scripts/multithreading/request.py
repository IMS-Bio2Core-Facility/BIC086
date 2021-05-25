# -*- coding: utf-8 -*-
"""Make a GET request to GTEx medianTranscriptExpression."""
import json
import logging
import threading

import requests

# from logs.get_logger import get_logger

thread_local = threading.local()
logger = logging.getLogger(__name__)


def get_session() -> requests.Session:
    """Instantiate a thread local session.

    Returns
    -------
    requests.Session

    """
    # session still worth it - re-used by each thread
    if not hasattr(thread_local, "session"):
        thread_local.session = requests.Session()
        thread_local.session.headers.update({"Accept": "application/json"})
        thread_local.session.params.update(  # type: ignore
            {
                "datasetId": "gtex_v8",
                "tissueSiteDetailId": "Brain_Hypothalamus",
                "hcluster": "False",
                "format": "json",
            }
        )
    return thread_local.session


def gtex_request(gene: str, output: str) -> None:
    """Make a gtex request against medianTranscriptExpression.

    Parameters
    ----------
    gene : str
        gene
    output : str
        output

    Raises
    ------
    requests.HTTPError
        When the GET request returns an error
    Exception
        Any other errors
    """
    s = get_session()
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
        message = response.json()
        with open(output, "w") as file:
            json.dump(message, file)
