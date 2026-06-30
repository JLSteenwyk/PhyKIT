"""Shared PhyKIT exception types."""


class PhykitUserError(SystemExit):
    """User-facing error that should be rendered by the CLI with a specific exit code."""

    def __init__(self, messages, code: int = 2):
        self.messages: list[str] = list(messages)
        super().__init__(code)
