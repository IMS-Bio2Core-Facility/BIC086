# -*- coding: utf-8 -*-
"""Make a GET request to GTEx medianTranscriptExpression."""
import concurrent.futures

from logs.get_logger import get_logger
from multithreading.request import gtex_request

LOGS = snakemake.log[0]  # noqa: F821
OUTS = snakemake.output  # noqa: F821
PARAMS = snakemake.params  # noqa: F821
THREADS = snakemake.threads  # noqa: F821


logger = get_logger(__name__, LOGS)

with concurrent.futures.ThreadPoolExecutor(max_workers=THREADS) as ex:
    ex.map(gtex_request, PARAMS["gene_ids"], OUTS["data"])
