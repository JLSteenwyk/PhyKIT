from __future__ import annotations

from ._fasta import read_unique_fasta_entries
from .base import Alignment


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class Faidx(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(fasta=parsed["fasta"], entry=parsed["entry"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        entries = self._parse_entries(self.entry)
        records = self._fetch_entries(self.fasta, entries)

        if self.json_output:
            rows = [
                {"entry": entry, "name": entry, "sequence": records[entry]}
                for entry in entries
            ]
            print_json(
                dict(
                    rows=rows,
                    entries=rows,
                )
            )
            return

        blocks = [f">{entry}\n{records[entry]}" for entry in entries]
        if blocks:
            print("\n".join(blocks))

    def process_args(self, args) -> dict[str, str]:
        return dict(
            fasta=args.fasta,
            entry=args.entry,
            json_output=getattr(args, "json", False),
        )

    @staticmethod
    def _fetch_entries(path: str, entries: list[str]) -> dict[str, str]:
        return read_unique_fasta_entries(path, entries)

    @staticmethod
    def _parse_entries(entry_arg: str) -> list[str]:
        if (
            " " not in entry_arg
            and "\t" not in entry_arg
            and "\n" not in entry_arg
            and "\r" not in entry_arg
        ):
            return [entry for entry in entry_arg.split(",") if entry]
        return [entry for entry in map(str.strip, entry_arg.split(",")) if entry]
