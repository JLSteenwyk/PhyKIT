"""Alignment subsampling: randomly subsample genes, partitions, or sites."""

from __future__ import annotations

import re

from ._fasta import read_fasta_first_token
from .base import Alignment
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


_FASTA_WRITE_CHUNK_ROWS = 8192
_PARTITION_WRITE_CHUNK_ROWS = 8192
_RANGE_SITE_SAMPLE_MIN_LENGTH = 1_000_000
_RANGE_SITE_SAMPLE_MID_MIN_LENGTH = 100_000
_RANGE_SITE_SAMPLE_MAX_COUNT = 50_000
_RANGE_SITE_SAMPLE_MAX_FRACTION = 0.1
_PARTITION_LINE_PATTERN = re.compile(
    r"\s*(\S+)\s*,\s*(\S+)\s*=\s*(\d+)\s*-\s*(\d+)\s*"
)


class AlignmentSubsample(Alignment):
    def __init__(self, args):
        parsed = self.process_args(args)
        super().__init__()
        self.mode = parsed["mode"]
        self.alignment_path = parsed["alignment_path"]
        self.list_path = parsed["list_path"]
        self.partition_path = parsed["partition_path"]
        self.number = parsed["number"]
        self.fraction = parsed["fraction"]
        self.seed = parsed["seed"]
        self.bootstrap = parsed["bootstrap"]
        self.output_prefix = parsed["output_prefix"]
        self.json_output = parsed["json_output"]

    def process_args(self, args):
        number = getattr(args, "number", None)
        fraction = getattr(args, "fraction", None)
        if number is None and fraction is None:
            raise PhykitUserError(
                ["Either --number or --fraction is required."], code=2
            )
        if number is not None and fraction is not None:
            raise PhykitUserError(
                ["--number and --fraction are mutually exclusive."], code=2
            )
        return dict(
            mode=args.mode,
            alignment_path=getattr(args, "alignment", None),
            list_path=getattr(args, "list", None),
            partition_path=getattr(args, "partition", None),
            number=number,
            fraction=fraction,
            seed=getattr(args, "seed", None),
            bootstrap=getattr(args, "bootstrap", False),
            output_prefix=getattr(args, "output", "subsampled"),
            json_output=getattr(args, "json", False),
        )

    def run(self):
        import random

        rng = random.Random(self.seed)

        if self.mode == "genes":
            self._run_genes(rng)
        elif self.mode == "partitions":
            self._run_partitions(rng)
        elif self.mode == "sites":
            self._run_sites(rng)
        else:
            raise PhykitUserError(
                [f"Unknown mode: {self.mode}. Use genes, partitions, or sites."],
                code=2,
            )

    # ------------------------------------------------------------------
    # genes mode
    # ------------------------------------------------------------------

    def _run_genes(self, rng):
        if self.list_path is None:
            raise PhykitUserError(
                ["--list is required for --mode genes."], code=2
            )

        genes = self._read_list_file(self.list_path)
        total = len(genes)
        n = self._compute_n(total, "genes")
        selected = self._sample(rng, genes, n)

        out_file = f"{self.output_prefix}.txt"
        with open(out_file, "w") as fh:
            fh.write("\n".join(selected))
            fh.write("\n")

        self._print_summary("genes", total, n, [out_file])

    # ------------------------------------------------------------------
    # partitions mode
    # ------------------------------------------------------------------

    def _run_partitions(self, rng):
        if self.alignment_path is None or self.partition_path is None:
            raise PhykitUserError(
                ["--alignment and --partition are required for --mode partitions."],
                code=2,
            )

        partitions = self._parse_partition_file(self.partition_path)
        sequences = self._read_alignment(self.alignment_path)
        alignment_length = self._alignment_length(sequences)
        self._validate_partition_bounds(partitions, alignment_length)
        total = len(partitions)
        n = self._compute_n(total, "partitions")
        selected = self._sample(rng, partitions, n)

        out_aln = f"{self.output_prefix}.fa"
        out_part = f"{self.output_prefix}.partition"

        new_sequences, new_partitions = self._assemble_partition_subsample(
            sequences,
            selected,
        )
        self._write_fasta(out_aln, new_sequences)
        self._write_partition_file(out_part, new_partitions)

        self._print_summary("partitions", total, n, [out_aln, out_part])

    # ------------------------------------------------------------------
    # sites mode
    # ------------------------------------------------------------------

    def _run_sites(self, rng):
        if self.alignment_path is None:
            raise PhykitUserError(
                ["--alignment is required for --mode sites."], code=2
            )

        sequences = self._read_alignment(self.alignment_path)
        aln_len = self._alignment_length(sequences)
        n = self._compute_n(aln_len, "sites")

        if not self.bootstrap and n == aln_len:
            new_sequences = sequences
        else:
            selected_indices = self._sample_site_indices(
                rng,
                aln_len,
                n,
                self.bootstrap,
            )
            selected_ranges = self._selected_index_ranges(selected_indices)
            if len(selected_ranges) * 4 < len(selected_indices):
                new_sequences: dict[str, str] = {
                    taxon: self._select_site_ranges(seq, selected_ranges)
                    for taxon, seq in sequences.items()
                }
            else:
                from operator import itemgetter

                site_selector = itemgetter(*selected_indices)
                new_sequences = {
                    taxon: self._select_sites(seq, site_selector)
                    for taxon, seq in sequences.items()
                }

        out_aln = f"{self.output_prefix}.fa"
        self._write_fasta(out_aln, new_sequences)

        self._print_summary("sites", aln_len, n, [out_aln])

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _read_list_file(path: str) -> list[str]:
        entries = []
        with open(path) as fh:
            for line in fh:
                stripped = line.strip()
                if stripped and stripped[0] != "#":
                    entries.append(stripped)
        return entries

    @staticmethod
    def _read_alignment(path: str) -> dict[str, str]:
        return read_fasta_first_token(path)

    @staticmethod
    def _alignment_length(sequences: dict[str, str]) -> int:
        seq_iter = iter(sequences.values())
        try:
            alignment_length = len(next(seq_iter))
        except StopIteration:
            raise PhykitUserError(
                ["All sequences in the alignment must have the same length."],
                code=2,
            )
        for sequence in seq_iter:
            if len(sequence) != alignment_length:
                raise PhykitUserError(
                    ["All sequences in the alignment must have the same length."],
                    code=2,
                )
        return alignment_length

    @staticmethod
    def _parse_partition_file(path: str) -> list[tuple[str, int, int]]:
        """Parse RAxML-style partition file.

        Expected format per line: ``AUTO, name=start-end``
        """
        partitions: list[tuple[str, int, int]] = []
        with open(path) as fh:
            for line_number, line in enumerate(fh, start=1):
                stripped = line.strip()
                if not stripped or stripped[0] == "#":
                    continue
                partition = AlignmentSubsample._parse_partition_line(stripped)
                if partition is None:
                    raise PhykitUserError(
                        [
                            f"Malformed partition row on line {line_number}: "
                            f"'{stripped}'.",
                            "Expected: MODEL, name=start-end",
                        ],
                        code=2,
                    )
                partitions.append(partition)
        return partitions

    @staticmethod
    def _parse_partition_line(line: str) -> tuple[str, int, int] | None:
        match = _PARTITION_LINE_PATTERN.fullmatch(line)
        if match is None:
            return None
        _model, name, start_text, end_text = match.groups()
        return name, int(start_text), int(end_text)

    @staticmethod
    def _validate_partition_bounds(
        partitions: list[tuple[str, int, int]],
        alignment_length: int,
    ) -> None:
        for name, start, end in partitions:
            if start < 1:
                raise PhykitUserError(
                    [
                        f"Partition '{name}' starts at {start}; "
                        "coordinates must be >= 1."
                    ],
                    code=2,
                )
            if end < start:
                raise PhykitUserError(
                    [
                        f"Partition '{name}' ends at {end} before its start "
                        f"at {start}."
                    ],
                    code=2,
                )
            if end > alignment_length:
                raise PhykitUserError(
                    [
                        f"Partition '{name}' ends at {end}, beyond alignment "
                        f"length {alignment_length}."
                    ],
                    code=2,
                )

    def _compute_n(self, total: int, label: str) -> int:
        if self.number is not None:
            n = self.number
        else:
            n = round(self.fraction * total)

        if n <= 0:
            raise PhykitUserError(
                [f"Computed number of {label} to select is 0. Increase --number or --fraction."],
                code=2,
            )
        if not self.bootstrap and n > total:
            raise PhykitUserError(
                [f"Cannot select {n} {label} from {total} without replacement."],
                code=2,
            )
        return n

    def _sample(self, rng, items, n):
        if self.bootstrap:
            return rng.choices(items, k=n)
        else:
            return rng.sample(items, k=n)

    @staticmethod
    def _sample_site_indices(
        rng,
        aln_len: int,
        n: int,
        bootstrap: bool,
    ) -> list[int]:
        if bootstrap:
            return rng.choices(range(aln_len), k=n)

        if (
            aln_len >= _RANGE_SITE_SAMPLE_MIN_LENGTH
            and n <= aln_len * _RANGE_SITE_SAMPLE_MAX_FRACTION
        ) or (
            aln_len >= _RANGE_SITE_SAMPLE_MID_MIN_LENGTH
            and n <= _RANGE_SITE_SAMPLE_MAX_COUNT
        ):
            return sorted(rng.sample(range(aln_len), k=n))

        return sorted(rng.sample(list(range(aln_len)), k=n))

    @staticmethod
    def _select_sites(seq: str, site_selector) -> str:
        selected = site_selector(seq)
        if isinstance(selected, str):
            return selected
        return "".join(selected)

    @staticmethod
    def _selected_index_ranges(indices: list[int]) -> list[tuple[int, int]]:
        iterator = iter(indices)
        try:
            start = previous = next(iterator)
        except StopIteration:
            return []

        ranges: list[tuple[int, int]] = []
        for index in iterator:
            if index == previous + 1:
                previous = index
                continue
            ranges.append((start, previous + 1))
            start = previous = index
        ranges.append((start, previous + 1))
        return ranges

    @staticmethod
    def _select_site_ranges(seq: str, ranges: list[tuple[int, int]]) -> str:
        if len(ranges) == 1:
            start, stop = ranges[0]
            return seq[start:stop]
        return "".join([seq[start:stop] for start, stop in ranges])

    @staticmethod
    def _assemble_partition_subsample(
        sequences: dict[str, str],
        selected: list[tuple[str, int, int]],
    ) -> tuple[dict[str, str], list[tuple[str, int, int]]]:
        selected_slices: list[slice] = []
        new_partitions: list[tuple[str, int, int]] = []
        current_pos = 1
        name_counts: dict[str, int] = {}

        for name, start, end in selected:
            start_idx = start - 1
            selected_slices.append(slice(start_idx, end))
            part_len = end - start_idx

            count = name_counts.get(name, 0)
            if count:
                display_name = f"{name}_dup{count}"
            else:
                display_name = name
            name_counts[name] = count + 1

            new_partitions.append(
                (display_name, current_pos, current_pos + part_len - 1)
            )
            current_pos += part_len

        if not selected_slices:
            new_sequences = {taxon: "" for taxon in sequences}
        elif len(selected_slices) == 1:
            selected_slice = selected_slices[0]
            new_sequences = {
                taxon: seq[selected_slice] for taxon, seq in sequences.items()
            }
        else:
            from operator import itemgetter

            slice_selector = itemgetter(*selected_slices)
            new_sequences = {
                taxon: "".join(slice_selector(seq))
                for taxon, seq in sequences.items()
            }
        return new_sequences, new_partitions

    @staticmethod
    def _write_fasta(path: str, sequences: dict[str, str]) -> None:
        from itertools import islice

        with open(path, "w") as fh:
            write = fh.write
            iterator = iter(sequences.items())
            while True:
                rows = [
                    f">{taxon}\n{seq}\n"
                    for taxon, seq in islice(iterator, _FASTA_WRITE_CHUNK_ROWS)
                ]
                if not rows:
                    break
                write("".join(rows))

    @staticmethod
    def _write_partition_file(
        path: str, partitions: list[tuple[str, int, int]]
    ) -> None:
        from itertools import islice

        with open(path, "w") as fh:
            write = fh.write
            iterator = iter(partitions)
            while True:
                rows = [
                    f"AUTO, {name}={start}-{end}\n"
                    for name, start, end in islice(
                        iterator,
                        _PARTITION_WRITE_CHUNK_ROWS,
                    )
                ]
                if not rows:
                    break
                write("".join(rows))

    def _print_summary(
        self, mode: str, total: int, selected: int, output_files: list[str]
    ):
        seed_display = self.seed if self.seed is not None else "random"
        bootstrap_str = "yes" if self.bootstrap else "no"

        if self.json_output:
            payload = {
                "mode": mode,
                "total": total,
                "selected": selected,
                "bootstrap": self.bootstrap,
                "seed": self.seed,
                "output_files": output_files,
            }
            print_json(payload)
            return

        lines = [
            "Alignment Subsampling",
            f"Mode: {mode}",
            f"Total {mode}: {total}",
            f"Selected: {selected}",
            f"Bootstrap: {bootstrap_str}",
            f"Seed: {seed_display}",
        ]
        lines.extend(f"Output: {f}" for f in output_files)

        try:
            print("\n".join(lines))
        except BrokenPipeError:
            pass
