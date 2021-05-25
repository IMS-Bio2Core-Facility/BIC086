# -*- coding: utf-8 -*-
"""Configure logger."""
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
