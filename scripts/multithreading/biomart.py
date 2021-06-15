# -*- coding: utf-8 -*-
"""Classes and functions for concurrent BioMart requests.

.. warning::

   Though every effort has been made to ensure thread safety,
   the concurrency is, as of yet, untested.

.. note::

   This module will likely undergo significant refactoring to generalise
   the concurrency pipeline and move all data handling code to a separate
   model.

Any API query is likely to be an I/O bound process,
particularly if there are many to make.
As this is a many in, many out process,
a call to ``concurrency.futures.ThreadPoolExecutor.map`` is sufficient,
though this call occurs in the analysis script.

To simplify the multi-step process for mapping,
and allow for multiple transcripts to be queried,
a ``BMSession`` class is provided.
"""
import logging
import threading
from dataclasses import dataclass
from io import StringIO
from time import sleep

import pandas as pd
from bioservices import BioMart

thread_local = threading.local()
logger = logging.getLogger(__name__)


def concurrent_biomart(file: str, output: str) -> None:
    """Given an input file, extract the ENST IDs and query BioMart.

    Note
    ----
    For use with ``concurrency.futures.ThreadPoolExecutor.map``.

    This is designed to allow concurrency within the given pipeline,
    as reading files and querying an external API are inherently I/O
    bound processes. Given the nature of this pipeline, the function
    expects the input file to contain a column `transcriptId` that
    contains a list of *version-less* ENST IDs from GTEx.

    Parameters
    ----------
    file : str
        Where to read the data from
    output : str
        Where to write the results to

    Raises
    ------
    IndexError
        If there are no transcipts to process in ``file``
    """
    logger.info(f"Processing file {file}")
    df = pd.read_csv(file, index_col=None)
    if len(transcripts := df.loc[:, "transcriptId"].unique().tolist()) == 0:
        logging.error(f'There are no transcripts for {file}. Raising index error')
        raise IndexError
    s = BMSession(transcripts, output)
    s.biomart_request()


@dataclass()
class BMSession:
    """A pipeline for requesting multiple transcripts from BioMart.

    Parameters
    ----------
    transcripts : list[str]
        A list of *version-less* ENST IDs
    output : str
        Where to save the csv results file
    """

    transcripts: list[str]
    output: str

    def __post_init__(self) -> None:
        """Initialise thread local BioMart session."""
        self._s = self._get_biomart()

    def _get_biomart(self) -> BioMart:
        """Instantiate a thread local BioMart session.

        Returns
        -------
        BioMart

        """
        # session still worth it - re-used by each thread
        if not hasattr(thread_local, "bm"):
            thread_local.bm = BioMart(host="www.ensembl.org", verbose=False, cache=True)
            logger.info(f"Generated BioMart with host {thread_local.bm.host}")
        return thread_local.bm

    def _query_biomart(self, transcript: str) -> pd.DataFrame:
        """Reguest data on a transcript from BioMart.

        A thread-local session is provided during post-initialisation.

        With the BioMart interface, simultaneous requests for multiple transcripts only
        return data for the last transcript. Ergo, this function must be mapped over all
        input transcipts individually.

        Parameters
        ----------
        transcript : str
            A *version-less* ENST ID

        Returns
        -------
        pd.DataFrame
            Containing the geneSymbol, gencodeId, transcriptId, and refseq ID of the
            transcript. If refseq is NA, then there is no corresponding ID in BioMart.

        Raises
        ------
        TimeoutError
            When BioMart cannot be reached after 5 tries

        """
        self._s.new_query()
        self._s.add_dataset_to_xml("hsapiens_gene_ensembl")
        self._s.add_filter_to_xml("ensembl_transcript_id", transcript)
        self._s.add_attribute_to_xml("hgnc_symbol")
        self._s.add_attribute_to_xml("ensembl_gene_id")
        self._s.add_attribute_to_xml("ensembl_transcript_id")
        self._s.add_attribute_to_xml("refseq_mrna")

        logger.info(f"Requesting transcript {transcript}")
        xml = self._s.get_xml()

        i = 0
        # Bioservices does not raise error, so we must check manually.
        # This can probably be improved...
        while "ERROR" in (message := self._s.query(xml)):
            logger.warning(f"Query error for {transcript}. Trying again in 0.1 s.")
            i += 1
            try:
                if i > 5:
                    raise TimeoutError
            except TimeoutError:
                logger.exception(
                    f"Query for {transcript} failed. Likely a BioMart connection issue."
                )
                raise
            else:
                sleep(0.1)
        else:
            data = pd.read_csv(
                StringIO(message),
                sep="\t",
                header=None,
                names=["geneSymbol", "gencodeId", "transcriptId", "refseq"],
            )
        return data

    def biomart_request(self) -> None:
        """Request data on multiple ENST IDs from BioMart."""
        data = pd.concat([self._query_biomart(t) for t in self.transcripts])
        data.to_csv(self.output, index=False)
        logger.info(f"Data saved to {self.output}")
