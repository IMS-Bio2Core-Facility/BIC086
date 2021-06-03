# -*- coding: utf-8 -*-
"""Configure logger.

This takes advantage of the parent-child relationship between modules to configure
logging.
This only needs to be called once,
in the top module -
here, the analysis scripts -
for these effects to be propagated to all modules that configure logging using:

.. code-block:: python

   import logging
   logger = logging.getLogger(__name__)

And then log using one of the standard levels.
"""
import logging


def get_logger(module: str, file: str) -> logging.Logger:
    """Configure a file logger for use in a script.

    Parameters
    ----------
    module : str
        The name of the module from which the logger is called
    file : str
        The name of the log file to which the logger will write


    Returns
    -------
    logging.Logger
        The configured logger instance.
    """
    handler = logging.FileHandler(file)
    formatter = logging.Formatter(
        "{asctime} :: {levelname} :: {name} :: {message}", style="{"
    )

    handler.setLevel(logging.INFO)
    handler.setFormatter(formatter)

    logging.basicConfig(level=logging.INFO, handlers=[handler])

    logger = logging.getLogger(module)

    return logger
