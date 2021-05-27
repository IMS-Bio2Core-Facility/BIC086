# -*- coding: utf-8 -*-
"""Make a GET request to GTEx medianTranscriptExpression."""
import logging
import threading

import requests

# from logs.get_logger import get_logger

thread_local = threading.local()
logger = logging.getLogger(__name__)


def get_session() -> requests.Session:
    """Instantiate a thread local session.

    The requests session is not thread safe, per `this`_ thread.
    To circumvent this, we create a thread local session. This means each session
    will still make multiple requests but remain isolated to its calling thread.

    Returns
    -------
    requests.Session

    .. _this:
       https://github.com/psf/requests/issues/2766
    """
    # session still worth it - re-used by each thread
    if not hasattr(thread_local, "session"):
        thread_local.session = requests.Session()
        thread_local.session.headers.update({"Accept": "text/html"})
        thread_local.session.params.update(  # type: ignore
            {
                "datasetId": "gtex_v8",
                "tissueSiteDetailId": "Brain_Hypothalamus",
                "format": "tsv",
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
        with open(output, "w") as file:
            file.write(response.text)
