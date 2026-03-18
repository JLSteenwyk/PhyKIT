"""Alignment subsampling: randomly subsample genes, partitions, or sites."""

import random
import re
from typing import Dict, List, Tuple

from Bio import SeqIO

from .base import Alignment
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


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
            for path in selected:
                fh.write(f"{path}\n")

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

        # Build new sequences by concatenating selected partition columns
        new_sequences: Dict[str, str] = {taxon: "" for taxon in sequences}
        new_partitions: List[Tuple[str, int, int]] = []
        current_pos = 1

        # Track duplicate names for bootstrap
        name_counts: Dict[str, int] = {}
        for name, start, end in selected:
            part_len = end - start + 1
            for taxon in sequences:
                new_sequences[taxon] += sequences[taxon][start - 1 : end]

            # Handle duplicate partition names (bootstrap)
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

        out_aln = f"{self.output_prefix}.fa"
        out_part = f"{self.output_prefix}.partition"

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

        indices = list(range(aln_len))
        if self.bootstrap:
            selected_indices = rng.choices(indices, k=n)
        else:
            selected_indices = sorted(rng.sample(indices, k=n))

        new_sequences: Dict[str, str] = {}
        for taxon, seq in sequences.items():
            new_sequences[taxon] = "".join(seq[i] for i in selected_indices)

        out_aln = f"{self.output_prefix}.fa"
        self._write_fasta(out_aln, new_sequences)

        self._print_summary("sites", aln_len, n, [out_aln])

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _read_list_file(path: str) -> List[str]:
        entries = []
        with open(path) as fh:
            for line in fh:
                stripped = line.strip()
                if stripped and not stripped.startswith("#"):
                    entries.append(stripped)
        return entries

    @staticmethod
    def _read_alignment(path: str) -> Dict[str, str]:
        sequences: Dict[str, str] = {}
        for record in SeqIO.parse(path, "fasta"):
            sequences[record.id] = str(record.seq)
        return sequences

    @staticmethod
    def _parse_partition_file(path: str) -> List[Tuple[str, int, int]]:
        """Parse RAxML-style partition file.

        Expected format per line: ``AUTO, name=start-end``
        """
        partitions: List[Tuple[str, int, int]] = []
        pattern = re.compile(r"^\s*\S+\s*,\s*(\S+)\s*=\s*(\d+)\s*-\s*(\d+)")
        with open(path) as fh:
            for line in fh:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                m = pattern.match(stripped)
                if m:
                    name = m.group(1)
                    start = int(m.group(2))
                    end = int(m.group(3))
                    partitions.append((name, start, end))
        return partitions

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
    def _write_fasta(path: str, sequences: Dict[str, str]) -> None:
        with open(path, "w") as fh:
            for taxon, seq in sequences.items():
                fh.write(f">{taxon}\n{seq}\n")

    @staticmethod
    def _write_partition_file(
        path: str, partitions: List[Tuple[str, int, int]]
    ) -> None:
        with open(path, "w") as fh:
            for name, start, end in partitions:
                fh.write(f"AUTO, {name}={start}-{end}\n")

    def _print_summary(self, mode: str, total: int, selected: int, output_files: List[str]):
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

        try:
            print("Alignment Subsampling")
            print(f"Mode: {mode}")
            print(f"Total {mode}: {total}")
            print(f"Selected: {selected}")
            print(f"Bootstrap: {bootstrap_str}")
            print(f"Seed: {seed_display}")
            for f in output_files:
                print(f"Output: {f}")
        except BrokenPipeError:
            pass
