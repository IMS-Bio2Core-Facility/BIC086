# -*- coding: utf-8 -*-
"""Process GTEx json data to xlsx."""
if __name__ == "__main__":
    import pandas as pd
    from logs.get_logger import get_logger
    from multithreading.process import Pipeline

    INS = snakemake.input  # noqa: F821
    LOGS = snakemake.log[0]  # noqa: F821
    OUTS = snakemake.output  # noqa: F821
    THREADS = snakemake.threads  # noqa: F821

    logger = get_logger(__name__, LOGS)

    mane = pd.read_csv(INS["mane"], index_col=0)

    with pd.ExcelWriter(OUTS["data"], engine="openpyxl") as writer:
        pl = Pipeline(
            gtex=INS["gtex"], bm=INS["bm"], writer=writer, mane=mane, maxworkers=THREADS
        )
        pl.run()
