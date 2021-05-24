# -*- coding: utf-8 -*-
"""Nox session configuration."""
from typing import List

import nox
from nox.sessions import Session

PACKAGE: str = "bic083"
LOCATIONS: List[str] = [
    PACKAGE,
    "noxfile.py",
    "tests",
]
VERSIONS: List[str] = ["3.8"]

nox.options.stop_on_first_error = True
nox.options.reuse_existing_virtualenvs = False


@nox.session(python="3.8")
def form(session: Session) -> None:
    """Format code with isort and black."""
    args = session.posargs or LOCATIONS
    session.install("isort==5.8.0", "black==21.4b2")
    session.run("isort", *args)
    session.run("black", *args)


@nox.session(python=VERSIONS)
def lint(session: Session) -> None:
    """Lint files with flake8."""
    args = session.posargs or LOCATIONS[:2]
    session.install(
        "flake8==3.9.1",
        "pyproject-flake8",
        "flake8-annotations==2.6.2",
        "flake8-bandit==2.1.2",
        "flake8-bugbear==21.4.3",
        "flake8-comprehensions==3.4.0",
        "flake8-docstrings==1.6.0",
        "flake8-pytest-style==1.4.1",
        "flake8-spellcheck==0.24.0",
        "darglint==1.8.0",
    )
    session.run("pflake8", *args)


@nox.session(python=VERSIONS)
def type(session: Session) -> None:
    """Type check files with mypy."""
    args = session.posargs or LOCATIONS
    session.install("mypy==0.812")
    session.run(
        "mypy",
        "--ignore-missing-imports",
        "--disable-error-code",
        "name-defined",
        *args
    )


@nox.session(python=VERSIONS)
def security(session: Session) -> None:
    """Check security safety."""
    session.run("safety", "check", external=True)


@nox.session(python=VERSIONS)
def tests(session: Session) -> None:
    """Run the test suite with pytest."""
    args = session.posargs or []
    session.install(
        "coverage==5.5",
        "pytest==6.2.4",
        "pytest-clarity==0.2.0a1",
        "pytest-sugar==0.9.4",
        "pytest-mock==3.6.0",
        "pytest-cov==2.11.1",
        "typeguard==2.12.0",  # Though typing, run best in pytest
        "six==1.15.0",
    )
    session.run("pytest", *args)
