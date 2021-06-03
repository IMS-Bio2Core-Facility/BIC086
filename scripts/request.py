# -*- coding: utf-8 -*-
"""Make a GET request to GTEx medianTranscriptExpression.

This step queries the GTEx API for transcript expression data in the region
specified by the user,
using a user provided list of genes names.
As the GTEx API is quite straightforward,
these queries can be made using the standard `requests.session`_ object.
Data were pulled from the ``gtex_v8`` dataset limited to the
region specified by the user.

There are a number of bugs in the GTEx API that necesitated some work arounds.
The data are returned in a TSV format,
as attempts to return in a CSV format produced a tab-delimited first line
and comma-delimited body.
Additionally,
the ``hclust`` option cannot be passed as it is force-set to ``true``,
no matter the user-passed value.
This is circumvented by not passing the value in the first place.

.. _requests.session: https://docs.python-requests.org/en/master/user/advanced/
"""
if __name__ == "__main__":
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

    genes = [
        lut.loc[lut["name"] == gene, "id"].values[0] for gene in PARAMS["gene_ids"]
    ]

    with concurrent.futures.ThreadPoolExecutor(max_workers=THREADS) as ex:
        ex.map(gtex_request, [PARAMS["region"]] * len(genes), genes, OUTS["data"])
