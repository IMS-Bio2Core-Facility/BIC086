# -*- coding: utf-8 -*-
"""Nox session configuration."""
import nox
from nox.sessions import Session

PACKAGE: str = "scripts"
LOCATIONS: list[str] = [
    PACKAGE,
    "noxfile.py",
]
VERSIONS: list[str] = ["3.9"]
CONDA_PARAMS: list[str] = ["-c", "bioconda", "-c", "conda-forge", "--file"]

nox.options.stop_on_first_error = True
nox.options.default_venv_backend = "conda"
nox.options.reuse_existing_virtualenvs = True
nox.options.sessions = ["form", "lint", "type", "security", "tests", "doc_tests"]


@nox.session(python=VERSIONS[0])
def form(session: Session) -> None:
    """Format code with isort and black."""
    args = session.posargs or LOCATIONS
    session.conda_install(*CONDA_PARAMS, "environments/form.txt")
    session.run("isort", *args)
    session.run("black", *args)


@nox.session(python=VERSIONS)
def lint(session: Session) -> None:
    """Lint files with flake8."""
    args = session.posargs or LOCATIONS
    session.conda_install(*CONDA_PARAMS, "environments/lint.txt")
    session.install("-r", "environments/lint_pip.txt", "--no-deps")
    session.run("pflake8", *args)


@nox.session(python=VERSIONS)
def type(session: Session) -> None:
    """Type check files with mypy."""
    args = session.posargs or LOCATIONS
    session.conda_install(*CONDA_PARAMS, "environments/type.txt")
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
    args = session.posargs or []
    session.conda_install(*CONDA_PARAMS, "environments/security.txt")
    session.run("safety", "check", "--bare", *args)


@nox.session(python=VERSIONS)
def tests(session: Session) -> None:
    """Run the test suite with pytest."""
    args = session.posargs or []
    session.conda_install(*CONDA_PARAMS, "environments/tests.txt")
    session.run("pytest", *args)


@nox.session(python=VERSIONS)
def doc_tests(session: Session) -> None:
    """Build the docs."""
    args = session.posargs or []
    session.conda_install(*CONDA_PARAMS, "environments/doc_tests.txt")
    session.run(
        "python",
        "-m",
        "xdoctest",
        "--verbose",
        "2",
        "--report",
        "cdiff",
        "--nocolor",
        PACKAGE,
        *args
    )


@nox.session(python="3.9")
def doc_build(session: Session) -> None:
    """Build the docs."""
    session.conda_install(*CONDA_PARAMS, "environments/doc_build.txt")
    session.install("-r", "environments/doc_build_pip.txt", "--no-deps")
    session.run("sphinx-build", "docs", "docs/_build")
