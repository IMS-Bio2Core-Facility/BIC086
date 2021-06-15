"""Pytest fixtures."""

from __future__ import annotations
import os
import pytest
from tempfile import NamedTemporaryFile
from types import TracebackType
from typing import Optional, Type

class CustomTempFile:
    """Create a temporary file with custom content.

    This is best used within a context manager to insure proper file clean up.

    Adapted from a StackOverflow `answer`_.

    .. _answer: https://stackoverflow.com/a/54053967

    Parameters
    ----------
    content : str
        Content of created temporary file.

    Examples
    --------
    >>> with CustomTempFile('Hello World') as tmp:
            with tmp.file as file:
                assert file.read() == 'Hello World'

    """

    def __init__(self, content: str) -> None:
        self.file = NamedTemporaryFile(mode="w", delete=False)
        with self.file as file:
            file.write(content)

    @property
    def filename(self) -> str:
        """Return the filename of the created file.

        Returns
        -------
        str
            The filename.
        """
        return self.file.name

    def __enter__(self) -> CustomTempFile:
        """Call on entry into ``with`` statement.

        Returns
        -------
        CustomTempFile
            Instance of self
        """
        return self

    def __exit__(
        self,
        ex_type: Optional[Type[BaseException]],
        ex_val: Optional[BaseException],
        tb: Optional[TracebackType],
    ) -> bool:
        """Call on exit from with statement.

        Typing via `mypy`_.

        .. _mypy: https://github.com/python/mypy/issues/4885

        Parameters
        ----------
        ex_type : Optional[Type[BaseException]]
            Exception type
        ex_val : Optional[BaseException]
            Exception value
        tb : Optional[TracebackType]
            Traceback

        Returns
        -------
        bool
            True, if successful.
        """
        os.unlink(self.filename)
        return True
