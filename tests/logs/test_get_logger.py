# -*- coding: utf-8 -*-
"""Tests for scripts.logs.get_logger.

Frustratingly, because of the way PyTest, logging, and the parent/child relationship
interact, it is difficult to test the undelying handler. For the time being, we suffice
to test that the logger exists, is enabled, and has the correct name.
"""
import logging
from pathlib import Path

from scripts.logs.get_logger import get_logger


def test_returns_logger(tmpdir: Path) -> None:
    """It returns a logger."""
    logger = get_logger(__name__, str(tmpdir / "tmp.log"))
    assert type(logger) == logging.Logger


def test_logger_enabled(tmpdir: Path) -> None:
    """It returns an enabled logger."""
    logger = get_logger(__name__, str(tmpdir / "tmp.log"))
    assert not logger.disabled


def test_logger_named(tmpdir: Path) -> None:
    """It returns a logger with the correct name."""
    logger = get_logger(__name__, str(tmpdir / "tmp.log"))
    assert logger.name == __name__
