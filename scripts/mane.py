# -*- coding: utf-8 -*-
"""Make a GET request to GTEx medianTranscriptExpression."""
if __name__ == "__main__":

    import pandas as pd
    from logs.get_logger import get_logger

    LOGS = snakemake.log[0]  # noqa: F821
    OUTS = snakemake.output  # noqa: F821
    PARAMS = snakemake.params  # noqa: F821

    logger = get_logger(__name__, LOGS)

    logger.info(f"Fetching MANE from {PARAMS['url']}")
    mane = pd.read_csv(PARAMS["url"], sep="\t", header=0)
    mane = mane.rename(
        columns={
            "Ensembl_Gene": "gencodeId",
            "symbol": "geneSymbol",
            "Ensembl_nuc": "transcriptId",
            "RefSeq_nuc": "refseq",
        }
    )

    mane.loc[:, ["gencodeId", "transcriptId"]] = mane.loc[
        :, ["gencodeId", "transcriptId"]
    ].apply(lambda x: x.str.split(".").str.get(0))

    mane.to_csv(OUTS["data"])
    logger.info(f"MANE saved to {OUTS['data']}")
