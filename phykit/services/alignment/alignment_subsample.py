"""Alignment subsampling: randomly subsample genes, partitions, or sites."""

from __future__ import annotations

from ._fasta import read_fasta_first_token
from .base import Alignment
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


_FASTA_WRITE_CHUNK_ROWS = 8192
_PARTITION_WRITE_CHUNK_ROWS = 8192


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
        lengths = set(len(seq) for seq in sequences.values())
        if len(lengths) != 1:
            raise PhykitUserError(
                ["All sequences in the alignment must have the same length."],
                code=2,
            )
        aln_len = lengths.pop()
        n = self._compute_n(aln_len, "sites")

        if not self.bootstrap and n == aln_len:
            new_sequences = dict(sequences)
        else:
            indices = list(range(aln_len))
            if self.bootstrap:
                selected_indices = rng.choices(indices, k=n)
            else:
                selected_indices = sorted(rng.sample(indices, k=n))

            from operator import itemgetter
            site_selector = itemgetter(*selected_indices)
            new_sequences: dict[str, str] = {
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
                if stripped and not stripped.startswith("#"):
                    entries.append(stripped)
        return entries

    @staticmethod
    def _read_alignment(path: str) -> dict[str, str]:
        return read_fasta_first_token(path)

    @staticmethod
    def _parse_partition_file(path: str) -> list[tuple[str, int, int]]:
        """Parse RAxML-style partition file.

        Expected format per line: ``AUTO, name=start-end``
        """
        partitions: list[tuple[str, int, int]] = []
        with open(path) as fh:
            for line in fh:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                partition = AlignmentSubsample._parse_partition_line(stripped)
                if partition is not None:
                    partitions.append(partition)
        return partitions

    @staticmethod
    def _parse_partition_line(line: str) -> tuple[str, int, int] | None:
        try:
            model, rest = line.split(", ", 1)
            name, range_part = rest.split("=", 1)
            start_text, end_text = range_part.split("-", 1)
            if (
                model
                and " " not in model
                and name
                and " " not in name
                and start_text.isdigit()
            ):
                if end_text.isdigit():
                    return name, int(start_text), int(end_text)
                end_digits = AlignmentSubsample._leading_digits(end_text)
                if end_digits.isdigit():
                    return name, int(start_text), int(end_digits)
        except (IndexError, ValueError):
            pass

        try:
            model, rest = line.split(",", 1)
            name_part, range_part = rest.lstrip().split("=", 1)
            start_part, end_part = range_part.lstrip().split("-", 1)
        except ValueError:
            return None

        if len(model.split()) != 1:
            return None

        name = name_part.rstrip()
        if not name or len(name.split()) != 1:
            return None

        start_text = start_part.strip()
        if not start_text.isdigit():
            return None

        end_text = end_part.lstrip()
        end_digits = AlignmentSubsample._leading_digits(end_text)
        if not end_digits.isdigit():
            return None

        return name, int(start_text), int(end_digits)

    @staticmethod
    def _leading_digits(text: str) -> str:
        idx = 0
        text_len = len(text)
        while idx < text_len and text[idx].isdigit():
            idx += 1
        return text[:idx]

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
    def _select_sites(seq: str, site_selector) -> str:
        selected = site_selector(seq)
        if isinstance(selected, str):
            return selected
        return "".join(selected)

    @staticmethod
    def _assemble_partition_subsample(
        sequences: dict[str, str],
        selected: list[tuple[str, int, int]],
    ) -> tuple[dict[str, str], list[tuple[str, int, int]]]:
        selected_ranges: list[tuple[int, int]] = []
        new_partitions: list[tuple[str, int, int]] = []
        current_pos = 1
        name_counts: dict[str, int] = {}

        for name, start, end in selected:
            start_idx = start - 1
            selected_ranges.append((start_idx, end))
            part_len = end - start_idx

            if name in name_counts:
                name_counts[name] += 1
                display_name = f"{name}_dup{name_counts[name]}"
            else:
                name_counts[name] = 0
                display_name = name

            new_partitions.append(
                (display_name, current_pos, current_pos + part_len - 1)
            )
            current_pos += part_len

        new_sequences = {
            taxon: "".join(seq[start:end] for start, end in selected_ranges)
            for taxon, seq in sequences.items()
        }
        return new_sequences, new_partitions

    @staticmethod
    def _write_fasta(path: str, sequences: dict[str, str]) -> None:
        with open(path, "w") as fh:
            rows: list[str] = []
            append = rows.append
            write = fh.write
            count = 0
            for taxon, seq in sequences.items():
                append(f">{taxon}\n{seq}\n")
                count += 1
                if count == _FASTA_WRITE_CHUNK_ROWS:
                    write("".join(rows))
                    rows.clear()
                    count = 0
            if rows:
                write("".join(rows))

    @staticmethod
    def _write_partition_file(
        path: str, partitions: list[tuple[str, int, int]]
    ) -> None:
        with open(path, "w") as fh:
            rows: list[str] = []
            append = rows.append
            write = fh.write
            count = 0
            for name, start, end in partitions:
                append(f"AUTO, {name}={start}-{end}\n")
                count += 1
                if count == _PARTITION_WRITE_CHUNK_ROWS:
                    write("".join(rows))
                    rows.clear()
                    count = 0
            if rows:
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
