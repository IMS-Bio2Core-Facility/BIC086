# -*- coding: utf-8 -*-
"""Make a GET request to GTEx medianTranscriptExpression."""
import json

import requests
from logs.get_logger import get_logger

LOGS = snakemake.log[0]  # noqa: F821
PARAMS = snakemake.params  # noqa: F821
OUTS = snakemake.output  # noqa: F821

logger = get_logger(__name__, LOGS)

for gene, output in zip(PARAMS["gene_ids"], OUTS["data"]):
    response = requests.get(
        "https://gtexportal.org/rest/v1/expression/medianTranscriptExpression",
        headers={"Accept": "application/json"},
        params={
            "datasetId": "gtex_v8",
            "gencodeId": gene,
            "tissueSiteDetailId": "Brain_Hypothalamus",
            "hcluster": "False",
            "format": "json",
        },
    )

    # By not raising, we ensure the whole loop runs, as no individual error is fatal
    try:
        response.raise_for_status()
    except requests.HTTPError:
        logger.exception(
            f"An HTTP Error occurred on gene {gene}. A detailed report follows"
        )
    except Exception:
        logger.exception(
            f"A non-HTTP Error occurred on gene {gene}. A detailed report follows"
        )
    else:
        logger.info(f"GET request successful! Details: \n\t{response.headers}")
        message = response.json()
        with open(output, "w") as file:
            json.dump(message, file)
