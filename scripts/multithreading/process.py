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
from threading import Event

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass()
class Pipeline:
    """A class-based implementation of a simple queue.

    This pipeline takes a list of files, performs basic data processing on them,
    and then concatenates the output to a single excel file. As such a process is
    IO bound (a quick check with pprofile suggests the majority of the time is
    writing to Excel), it is an excellent choice for concurrency.

    The producer will return when there are no more files, and the consumer will
    returen when the queue is joined. The former is tracked with a simple while
    loop and `list.pop`, while the latter is tracked using a `threading.Event`.

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

    files: list[str]
    writer: pd.ExcelWriter
    maxworkers: int
    maxsize: int = 0

    def __post_init__(self) -> None:
        """Initialise queue from given values."""
        self._q: queue.Queue = queue.Queue(maxsize=self.maxsize)

    def _producer(self) -> None:
        """Initialise the queue's producer.

        The producer will return when all files have been iterated over.
        This has the added advantage of terminating immediately if no files are passed.
        """
        while len(self.files) > 0:
            with open(self.files.pop(0), "r") as file:
                raw = json.load(file)
            data = pd.DataFrame.from_dict(raw["medianTranscriptExpression"])
            data = data.sort_values("median", ascending=False)
            self._q.put(data)

    def _consumer(self, event: Event) -> None:
        """Initialise the queue's consumer.

        The consumer will return when event is triggered.
        This is achieved by using `self._q.task_done` in this method and
        `self._q.join` in the run method.

        Parameters
        ----------
        event : Event
            An event marking whether the queue has joined.
        """
        while not event.is_set():
            data = self._q.get()
            gene = data["geneSymbol"].unique()[0]
            data.to_excel(self.writer, index=False, sheet_name=gene)
            self._q.task_done()
        else:
            logging.info("Queue consumed")

    def run(self) -> None:
        """Run the queue.

        The consumer will return when there are no more files, and the producer
        will not return until the queue has joined. The latter is triggered using
        a `threading.Event` afte a call to `queue.Queue.join`.
        """
        is_consumed = Event()
        with ThreadPoolExecutor(max_workers=self.maxworkers) as ex:
            ex.submit(self._producer)
            ex.submit(self._consumer, is_consumed)
            self._q.join()
            is_consumed.set()
