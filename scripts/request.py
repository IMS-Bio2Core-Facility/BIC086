# -*- coding: utf-8 -*-
"""Make a GET request to GTEx medianTranscriptExpression."""
import concurrent.futures

import pandas as pd
from logs.get_logger import get_logger
from multithreading.request import gtex_request

INPUTS = snakemake.input  # noqa: F821
LOGS = snakemake.log[0]  # noqa: F821
OUTS = snakemake.output  # noqa: F821
PARAMS = snakemake.params  # noqa: F821
THREADS = snakemake.threads  # noqa: F821


logger = get_logger(__name__, LOGS)

lut = pd.read_csv(
    INPUTS["lut"], sep=" ", index_col=None, header=None, names=["id", "name"]
)

genes = [lut.loc[lut["name"] == gene, "id"].values[0] for gene in PARAMS["gene_ids"]]

with concurrent.futures.ThreadPoolExecutor(max_workers=THREADS) as ex:
    ex.map(gtex_request, genes, OUTS["data"])
