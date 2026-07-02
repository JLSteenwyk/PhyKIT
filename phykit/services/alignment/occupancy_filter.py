"""
Filter alignments and/or trees by cross-file taxon occupancy.

Given a list of alignment or tree files, counts how many files each
taxon appears in and retains only taxa meeting a minimum occupancy
threshold. Outputs filtered copies of each input file.
"""
from __future__ import annotations

from collections import Counter
import os
from pathlib import Path

from ._fasta import _clean_sequence, read_fasta_first_tokens
from ...errors import PhykitUserError


_path_exists = os.path.exists


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def _normalize_list_path(path: str) -> str:
    separator = os.sep
    if not (
        separator + separator in path
        or path == "."
        or path.startswith("." + separator)
        or separator + "." + separator in path
        or path.endswith(separator + ".")
    ):
        return path

    is_absolute = os.path.isabs(path)
    if is_absolute:
        if (
            separator == "/"
            and path.startswith("//")
            and not path.startswith("///")
        ):
            prefix = "//"
            rest = path[2:]
        else:
            prefix = separator
            rest = path.lstrip(separator)
    else:
        prefix = ""
        rest = path

    parts = [
        part for part in rest.split(separator)
        if part and part != "."
    ]
    normalized = separator.join(parts)
    if prefix:
        return prefix + normalized if normalized else prefix
    return normalized or "."


class OccupancyFilter:
    """Filter files to retain only taxa meeting an occupancy threshold."""

    def __init__(self, args):
        parsed = self.process_args(args)
        self.list_path = parsed["list_path"]
        self.format = parsed["format"]
        self.threshold = parsed["threshold"]
        self.output_dir = parsed["output_dir"]
        self.suffix = parsed["suffix"]
        self.json_output = parsed["json_output"]

    def process_args(self, args):
        return dict(
            list_path=args.list,
            format=getattr(args, "format", "fasta"),
            threshold=getattr(args, "threshold", 0.5),
            output_dir=getattr(args, "output_dir", None),
            suffix=getattr(args, "suffix", ".filtered"),
            json_output=getattr(args, "json", False),
        )

    def _resolve_threshold(self, total_files: int) -> int:
        """Convert threshold to an integer count.

        If threshold is between 0 and 1 (inclusive), treat as a fraction
        of total files (e.g., 0.5 = 50%, 1.0 = 100%). Values > 1 are
        treated as an absolute count.
        """
        t = self.threshold
        if 0 < t <= 1:
            import math
            return max(1, math.ceil(t * total_files))
        return int(t)

    def run(self):
        file_paths = self._read_file_list(self.list_path)

        if not file_paths:
            raise PhykitUserError(["No files found in list."], code=2)

        # Extract taxa from each file
        file_taxa = {}
        for path in file_paths:
            taxa = self._extract_taxa(path)
            file_taxa[path] = set(taxa)

        # Count occupancy across files
        occupancy = Counter()
        update_occupancy = occupancy.update
        for taxa in file_taxa.values():
            update_occupancy(taxa)

        total_files = len(file_paths)
        min_count = self._resolve_threshold(total_files)

        # Determine which taxa pass the threshold
        kept_taxa = set()
        removed_taxa = set()
        for taxon, count in occupancy.items():
            if count >= min_count:
                kept_taxa.add(taxon)
            else:
                removed_taxa.add(taxon)

        if not kept_taxa:
            raise PhykitUserError(
                [
                    f"No taxa meet the occupancy threshold of {min_count} "
                    f"(from {self.threshold}).",
                    f"Maximum occupancy is {max(occupancy.values())} "
                    f"across {total_files} files.",
                ],
                code=2,
            )

        # Determine output directory
        if self.output_dir:
            out_dir = Path(self.output_dir)
            out_dir.mkdir(parents=True, exist_ok=True)
        else:
            out_dir = None

        # Filter and write each file
        output_paths = []
        files_modified = 0
        taxa_removed_per_file = {}

        for path in file_paths:
            taxa_in_file = file_taxa[path]
            taxa_to_remove = taxa_in_file - kept_taxa

            out_path = self._output_path(path, out_dir)
            output_paths.append(out_path)

            if self.format == "fasta":
                n_kept, n_removed = self._filter_fasta(
                    path, out_path, kept_taxa
                )
            else:
                n_kept, n_removed = self._filter_tree(
                    path, out_path, kept_taxa
                )

            taxa_removed_per_file[path] = sorted(taxa_to_remove)
            if n_removed > 0:
                files_modified += 1

        # Output
        if self.json_output:
            self._print_json(
                file_paths, output_paths, occupancy, kept_taxa,
                removed_taxa, taxa_removed_per_file, total_files,
                min_count,
            )
        else:
            self._print_text(
                file_paths, output_paths, occupancy, kept_taxa,
                removed_taxa, taxa_removed_per_file, total_files,
                files_modified, min_count,
            )

    def _read_file_list(self, path):
        """Read file paths from a list file."""
        source = Path(path)
        if not source.exists():
            raise PhykitUserError(
                [f"{path} corresponds to no such file or directory."],
                code=2,
            )

        paths = []
        append = paths.append
        parent_str = str(source.parent)
        parent_prefix = "" if parent_str == "." else parent_str + os.sep
        separator = os.sep
        double_separator = separator + separator
        dot_prefix = "." + separator
        dot_segment = separator + "." + separator
        dot_suffix = separator + "."
        check_isabs = os.path.isabs if separator == "\\" else None
        with source.open() as handle:
            for line in handle:
                line = line.strip()
                if not line or line[0] == "#":
                    continue
                if line == ".":
                    append(parent_str)
                    continue
                has_separator = separator in line
                if not has_separator:
                    append(parent_prefix + line)
                    continue
                needs_normalization = (
                    line.startswith(dot_prefix)
                    or double_separator in line
                    or dot_segment in line
                    or line.endswith(dot_suffix)
                )
                if line[0] == separator or (
                    check_isabs is not None and check_isabs(line)
                ):
                    append(
                        _normalize_list_path(line)
                        if needs_normalization
                        else line
                    )
                elif needs_normalization:
                    normalized = _normalize_list_path(line)
                    if normalized == ".":
                        append(parent_str)
                    else:
                        append(parent_prefix + normalized)
                else:
                    append(parent_prefix + line)
        return paths

    def _extract_taxa(self, path):
        """Extract taxon names from a tree or FASTA file."""
        if not _path_exists(path):
            raise PhykitUserError(
                [f"{path} corresponds to no such file or directory."],
                code=2,
            )

        if self.format == "fasta":
            try:
                return read_fasta_first_tokens(path)
            except Exception:
                raise PhykitUserError(
                    [f"Could not parse FASTA file: {path}"], code=2
                )
        else:
            from Bio import Phylo
            from ..tree.base import Tree

            try:
                tree = Phylo.read(path, "newick")
                taxa = Tree.calculate_terminal_names_fast(tree)
                if taxa is not None:
                    return taxa
                return [t.name for t in tree.get_terminals()]
            except Exception:
                raise PhykitUserError(
                    [f"Could not parse tree file: {path}"], code=2
                )

    def _output_path(self, input_path, out_dir):
        """Determine output file path."""
        p = Path(input_path)
        stem = p.stem
        ext = p.suffix
        out_name = f"{stem}{self.suffix}{ext}"
        if out_dir:
            return str(out_dir / out_name)
        return str(p.parent / out_name)

    def _filter_fasta(self, input_path, output_path, kept_taxa):
        """Filter a FASTA file to keep only taxa in kept_taxa."""
        total = 0
        kept = 0

        with open(output_path, "w") as f, open(input_path) as handle:
            write = f.write
            title = None
            record_kept = False
            sequence_parts = []

            for line in handle:
                if line[0] == ">":
                    if record_kept:
                        write(f">{title}\n")
                        self._write_wrapped_fasta_sequence(
                            f,
                            _clean_sequence(sequence_parts),
                        )
                    title = line[1:].rstrip()
                    total += 1
                    record_id = title.split(None, 1)[0]
                    record_kept = record_id in kept_taxa
                    sequence_parts = []
                    if record_kept:
                        kept += 1
                elif record_kept:
                    sequence_parts.append(line.rstrip())

            if record_kept:
                write(f">{title}\n")
                self._write_wrapped_fasta_sequence(
                    f,
                    _clean_sequence(sequence_parts),
                )

        return kept, total - kept

    @staticmethod
    def _write_wrapped_fasta_sequence(handle, sequence: str, width: int = 60) -> None:
        if not sequence:
            return
        if len(sequence) <= width:
            handle.write(sequence)
            handle.write("\n")
            return
        handle.write(
            "\n".join(
                [
                    sequence[idx:idx + width]
                    for idx in range(0, len(sequence), width)
                ]
            )
        )
        handle.write("\n")

    def _filter_tree(self, input_path, output_path, kept_taxa):
        """Filter a tree to keep only tips in kept_taxa."""
        from Bio import Phylo
        from ..tree.base import Tree

        tree = Phylo.read(input_path, "newick")
        terminals = []
        tips_to_prune = []
        target_ids = set()
        try:
            root = tree.root
            root.clades
        except AttributeError:
            root = None

        if root is not None:
            stack = [root]
            pop = stack.pop
            append = stack.append
            append_terminal = terminals.append
            append_prune = tips_to_prune.append
            add_target = target_ids.add
            contains = kept_taxa.__contains__
            try:
                while stack:
                    clade = pop()
                    children = clade.clades
                    if children:
                        child_count = len(children)
                        if child_count == 2:
                            append(children[1])
                            append(children[0])
                        else:
                            for index in range(child_count - 1, -1, -1):
                                append(children[index])
                    else:
                        append_terminal(clade)
                        if not contains(clade.name):
                            append_prune(clade)
                            add_target(id(clade))
            except AttributeError:
                terminals = []
                tips_to_prune = []
                target_ids = set()

        if root is None or (not terminals and not tips_to_prune):
            terminals = tree.get_terminals()
            tips_to_prune = [
                tip for tip in terminals if tip.name not in kept_taxa
            ]

        if tips_to_prune and target_ids:
            if not Tree._prune_terminal_objects_batch_standard_tree(
                tree,
                target_ids,
            ):
                for tip in tips_to_prune:
                    try:
                        tree.prune(tip)
                    except Exception:
                        pass
        elif tips_to_prune:
            for tip in tips_to_prune:
                try:
                    tree.prune(tip)
                except Exception:
                    pass
        Phylo.write(tree, output_path, "newick")

        return len(terminals) - len(tips_to_prune), len(tips_to_prune)

    def _print_text(
        self, file_paths, output_paths, occupancy, kept_taxa,
        removed_taxa, taxa_removed_per_file, total_files, files_modified,
        min_count,
    ):
        try:
            threshold_str = (
                f"{self.threshold} = {min_count}/{total_files} files"
                if 0 < self.threshold < 1
                else f"{min_count}/{total_files} files"
            )
            lines = [
                f"Occupancy filter (threshold: {threshold_str})",
                f"Input files: {total_files}",
                f"Total taxa: {len(occupancy)}",
                f"Taxa kept: {len(kept_taxa)}",
                f"Taxa removed: {len(removed_taxa)}",
                f"Files modified: {files_modified}",
                "",
            ]

            if removed_taxa:
                lines.append("Removed taxa (occupancy):")
                for taxon in sorted(removed_taxa):
                    lines.append(f"  {taxon}: {occupancy[taxon]}/{total_files}")
                lines.append("")

            lines.append("Occupancy per taxon:")
            for taxon in sorted(occupancy):
                status = "kept" if taxon in kept_taxa else "REMOVED"
                lines.append(f"  {taxon}: {occupancy[taxon]}/{total_files} ({status})")
            lines.append("")

            lines.append("Output files:")
            for in_path, out_path in zip(file_paths, output_paths):
                n_removed = len(taxa_removed_per_file.get(in_path, []))
                lines.append(f"  {out_path}" + (
                    f" ({n_removed} taxa removed)" if n_removed else ""
                ))

            print("\n".join(lines))

        except BrokenPipeError:
            return

    def _print_json(
        self, file_paths, output_paths, occupancy, kept_taxa,
        removed_taxa, taxa_removed_per_file, total_files, min_count,
    ):
        payload = {
            "threshold": self.threshold,
            "threshold_count": min_count,
            "total_files": total_files,
            "total_taxa": len(occupancy),
            "taxa_kept": len(kept_taxa),
            "taxa_removed": len(removed_taxa),
            "kept_taxa": sorted(kept_taxa),
            "removed_taxa": sorted(removed_taxa),
            "occupancy": {
                taxon: occupancy[taxon]
                for taxon in sorted(occupancy)
            },
            "files": [
                {
                    "input": in_path,
                    "output": out_path,
                    "taxa_removed": taxa_removed_per_file.get(in_path, []),
                }
                for in_path, out_path in zip(file_paths, output_paths)
            ],
        }
        print_json(payload, sort_keys=False)
