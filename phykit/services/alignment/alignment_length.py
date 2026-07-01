from __future__ import annotations

from .base import Alignment


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


_FASTA_HEADER_BYTE = ord(">")
_SPACE_BYTE = ord(" ")
_CARRIAGE_RETURN_BYTE = ord("\r")
_TRAILING_WHITESPACE_BYTES = (9, 11, 12, 13, 32)


class AlignmentLength(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        aln_len = self._get_fasta_alignment_length(self.alignment_file_path)
        if aln_len is None:
            alignment, _, _ = self.get_alignment_and_format()
            aln_len = alignment.get_alignment_length()
        if self.json_output:
            print_json(dict(alignment_length=aln_len))
            return
        print(aln_len)

    def process_args(self, args) -> dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
        )

    @staticmethod
    def _get_fasta_alignment_length(path: str) -> int | None:
        try:
            with open(path, "rb") as handle:
                first_line = handle.readline()
                if not first_line.startswith(b">"):
                    return None

                aln_len = None
                seq_len = 0
                for line in handle:
                    if line[0] == _FASTA_HEADER_BYTE:
                        if aln_len is None:
                            aln_len = seq_len
                        elif seq_len != aln_len:
                            return None
                        seq_len = 0
                        continue

                    if not line.isascii():
                        return None

                    line_len = len(line)
                    if line_len and line[-1] == 10:
                        line_len -= 1

                    if (
                        _SPACE_BYTE in line
                        or _CARRIAGE_RETURN_BYTE in line
                        or (
                            line_len
                            and line[line_len - 1] in _TRAILING_WHITESPACE_BYTES
                        )
                    ):
                        sequence_line = line.rstrip()
                        if _SPACE_BYTE in sequence_line:
                            sequence_line = sequence_line.replace(b" ", b"")
                        if _CARRIAGE_RETURN_BYTE in sequence_line:
                            sequence_line = sequence_line.replace(b"\r", b"")
                        seq_len += len(sequence_line)
                    else:
                        seq_len += line_len

                if aln_len is None:
                    aln_len = seq_len
                elif seq_len != aln_len:
                    return None
        except OSError:
            return None

        return aln_len
