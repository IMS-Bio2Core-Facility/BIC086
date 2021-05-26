# -*- coding: utf-8 -*-
"""Multi-thread pandas processing.

Writing to files is definitely not thread safe. But here, we have multiple read-ins
that are independent. Additionally, as the bottleneck is definitely the read step,
(ie this is I/O bound), we want concurrency. Here, we implement a simple queue to
allow exactly this.
"""
import json
import logging
import queue
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Union

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass()
class Pipeline:
    """A class-based implementation of a simple queue.

    This pipeline takes a list of files, performs basic data processing on them,
    and then concatenates the output to a single excel file. As such a process is
    IO bound (a quick check with pprofile suggests the majority of the time is
    writing to Excel), it is an excellent choice for concurrency.

    Using joining has the interesting problem that that the consumer can finish
    before all items have been processed by the producer. Instead, we trigger
    each to return by marking the last item as None.

    Attributes
    ----------
    files : list[str]
        The list of files to iterate over.
    writer : pd.ExcelWriter
        The *open* pd.ExcelWriter instance.
        This should be opened, preferably with a `with` block before instantiating
            the class.
    maxworkers : int
        The maximum number of workers to use.
    maxsize : int, default: 0
        The maximum size of the queue.
        Defaults to no limit.
    """

    files: list[Union[str, None]]
    writer: pd.ExcelWriter
    maxworkers: int
    maxsize: int = 0

    def __post_init__(self) -> None:
        """Initialise queue from given values."""
        self.files += [None]
        self._q: queue.Queue = queue.Queue(maxsize=self.maxsize)

    def _producer(self) -> None:
        """Initialise the queue's producer.

        The producer will return when all files have been iterated over, as indicated
        by the presence of `None` appended in `__post_init__`. This works by popping
        from 0 instead of -1.

        Additionally, the None signal is sent to the consumer to indicate when there are
        no more items to process.
        """
        while (path := self.files.pop(0)) is not None:
            with open(path, "r") as file:
                raw = json.load(file)
            data = pd.DataFrame.from_dict(raw["medianTranscriptExpression"])
            data = data.sort_values("median", ascending=False)
            self._q.put(data)
        else:
            self._q.put(path)  # Send end signal to consumer
            return

    def _consumer(self) -> None:
        """Initialise the queue's consumer.

        The consumer will return when producer sends None. This works as the default
        queue is a FIFO queue.
        """
        while (data := self._q.get()) is not None:
            gene = data["geneSymbol"].unique()[0]
            data.to_excel(self.writer, index=False, sheet_name=gene)
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
