"""
Global thread limit shared by every PhyKIT command.

The limit comes from the ``--threads`` option (accepted by every command) or
the ``PHYKIT_THREADS`` environment variable. It caps the number of worker
processes PhyKIT starts and the number of threads used by numeric libraries
(BLAS/OpenMP via NumPy, SciPy, scikit-learn, numba). Those libraries read
their environment variables when first imported, which in PhyKIT always
happens after argument parsing because services are loaded lazily.
"""

from __future__ import annotations

import os

from ..errors import PhykitUserError


PHYKIT_THREADS_ENV = "PHYKIT_THREADS"
THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "NUMBA_NUM_THREADS",
)
_OPTION = "--threads"


def _parse_positive_int(value: str | None) -> int | None:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return None
    return parsed if parsed >= 1 else None


def _invalid_value_error(value: str | None) -> PhykitUserError:
    shown = "nothing" if value is None else repr(value)
    return PhykitUserError(
        [f"{_OPTION} requires a positive integer, got {shown}"], code=2
    )


def extract_threads_option(argv: list[str]) -> tuple[int | None, list[str]]:
    """Remove ``--threads N`` / ``--threads=N`` from argv.

    Returns the requested thread count (last occurrence wins) and the
    remaining arguments.
    """
    value = None
    rest = []
    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == _OPTION:
            raw = argv[i + 1] if i + 1 < len(argv) else None
            i += 2
        elif arg.startswith(_OPTION + "="):
            raw = arg[len(_OPTION) + 1:]
            i += 1
        else:
            rest.append(arg)
            i += 1
            continue
        value = _parse_positive_int(raw)
        if value is None:
            raise _invalid_value_error(raw)
    return value, rest


def apply_thread_limit(n: int, override: bool = True) -> None:
    """Export the thread limit so PhyKIT, numeric libraries, and worker
    processes all see it."""
    os.environ[PHYKIT_THREADS_ENV] = str(n)
    for name in THREAD_ENV_VARS:
        if override:
            os.environ[name] = str(n)
        else:
            os.environ.setdefault(name, str(n))


def thread_limit() -> int | None:
    """Return the active thread limit, or None when unlimited."""
    return _parse_positive_int(os.environ.get(PHYKIT_THREADS_ENV))


def limit_workers(n: int) -> int:
    """Cap a worker count at the active thread limit (never below 1)."""
    limit = thread_limit()
    if limit is not None:
        n = min(n, limit)
    return max(1, n)


def configure_threads(argv: list[str]) -> list[str]:
    """Apply ``--threads`` (or ``PHYKIT_THREADS``) and return argv without it.

    An explicit ``--threads`` overrides any numeric-library variables already
    set; ``PHYKIT_THREADS`` alone only fills in the ones that are unset.
    """
    value, rest = extract_threads_option(argv)
    if value is not None:
        apply_thread_limit(value)
    else:
        env_value = thread_limit()
        if env_value is not None:
            apply_thread_limit(env_value, override=False)
    return rest
