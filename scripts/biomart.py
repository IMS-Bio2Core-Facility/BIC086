# -*- coding: utf-8 -*-
"""Make a GET request to GTEx medianTranscriptExpression."""
if __name__ == "__main__":
    import concurrent.futures

    # import pandas as pd
    from logs.get_logger import get_logger
    from multithreading.biomart import concurrent_biomart

    INPUTS = snakemake.input  # noqa: F821
    LOGS = snakemake.log[0]  # noqa: F821
    OUTS = snakemake.output  # noqa: F821
    THREADS = snakemake.threads  # noqa: F821

    logger = get_logger(__name__, LOGS)

    with concurrent.futures.ThreadPoolExecutor(max_workers=THREADS) as ex:
        ex.map(concurrent_biomart, INPUTS["data"], OUTS["data"])
