"""Shared PhyKIT exception types."""

from typing import Iterable, List


class PhykitUserError(SystemExit):
    """User-facing error that should be rendered by the CLI with a specific exit code."""

    def __init__(self, messages: Iterable[str], code: int = 2):
        self.messages: List[str] = list(messages)
        super().__init__(code)

