# -*- coding: utf-8 -*-
"""A simple queue for data processing.

.. warning::

   Though every effort has been made to ensure thread safety,
   the concurrency is, as of yet, untested.

.. note::

   This module will likely undergo significant refactoring to generalise
   the concurrency pipeline and move all data handling code to a separate
   model.

Combining the relevant data is a massively I/O bound process.
That is,
the file reading/writing takes more time than the data processing,
making this an ideal candidate for concurrency.
Additionally,
as this is a "many in, one out" situation,
a queue is a useful structure.
Here,
a common producer-consumer queue is structured as a class to simplify calls
to the queue.
"""
import logging
import queue
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Union

import pandas as pd
from data_handling.process import merge_data, write_data

logger = logging.getLogger(__name__)


@dataclass()
class Pipeline:
    """A class-based implementation of a simple queue.

    This pipeline takes a list of files, performs basic data processing on them,
    and then concatenates the output to a single excel file. As such a process is
    IO bound (a quick check with pprofile suggests the majority of the time is
    writing to Excel), it is an excellent choice for concurrency.

    We take advantage of the ``dataclass __post_init__`` method to initialise the queue
    with user-specified parameters,
    and to mark the end of the files with None.

    Using joining has the interesting problem that that the consumer can finish
    before all items have been processed by the producer. Instead, we trigger
    each to return by marking the last item as None.

    Attributes
    ----------
    gtex : list[str]
        The list of gtex data to iterate over.
    bm : list[str]
        The list of BioMart data to iterate over.
    writer : pd.ExcelWriter
        The *open* pd.ExcelWriter instance.
        This should be opened, preferably with a `with` block before instantiating
        the class.
    mane : pd.DataFrame
        The MANE reference data to be used.
    maxworkers : int
        The maximum number of workers to use.
    maxsize : int, default: 0
        The maximum size of the queue.
        Defaults to no limit.
    """

    gtex: list[Union[str, None]]
    bm: list[Union[str, None]]
    writer: pd.ExcelWriter
    mane: pd.DataFrame
    maxworkers: int
    maxsize: int = 0

    def __post_init__(self) -> None:
        """Initialise queue from given values and end lists with ``None``."""
        self.gtex += [None]
        self.bm += [None]
        self._q: queue.Queue = queue.Queue(maxsize=self.maxsize)

    def _producer(self) -> None:
        """Initialise the queue's producer.

        The producer will return when all files have been iterated over, as indicated
        by the presence of `None` appended in `__post_init__`. This works by popping
        from 0 instead of -1.

        Additionally, the None signal is sent to the consumer to indicate when there are
        no more items to process.
        """
        while (gtex_path := self.gtex.pop(0)) is not None and (
            bm_path := self.bm.pop(0)
        ) is not None:
            data = merge_data(gtex_path, bm_path, self.mane)
            self._q.put(data)
            logger.info(f"Contents of file {gtex_path} added to queue")
        else:
            self._q.put(None)  # Send end signal to consumer
            logger.info("All files added. None signal sent. Producer returns")
            return

    def _consumer(self) -> None:
        """Initialise the queue's consumer.

        The consumer will return when producer sends None. This works as the default
        queue is a FIFO queue.
        """
        while (data := self._q.get()) is not None:
            write_data(data, self.writer)
            self._q.task_done()
        else:
            logging.info("None received. Queue consumed.")
            self._q.task_done()
            return

    def run(self) -> None:
        """Run the queue.

        The producer will return when there are no more files, and the consumer
        will return when the producer sends `None`. This prevent the interesting,
        and entirely common, situation of the consumer returning before the producer
        finishes the known number of inputs. A call to `join` is still used, however,
        to insure the queue locks till done.
        """
        with ThreadPoolExecutor(max_workers=self.maxworkers) as ex:
            ex.submit(self._producer)
            ex.submit(self._consumer)
            self._q.join()
