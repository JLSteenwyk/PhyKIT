"""
Filter alignments and/or trees by cross-file taxon occupancy.

Given a list of alignment or tree files, counts how many files each
taxon appears in and retains only taxa meeting a minimum occupancy
threshold. Outputs filtered copies of each input file.
"""
from collections import Counter
from pathlib import Path
from typing import Dict, List

from ...errors import PhykitUserError
from ...helpers.json_output import print_json


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
        for taxa in file_taxa.values():
            for taxon in taxa:
                occupancy[taxon] += 1

        total_files = len(file_paths)
        min_count = self._resolve_threshold(total_files)
        all_taxa = sorted(occupancy.keys())

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

        lines = source.read_text().splitlines()
        paths = []
        for line in lines:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = Path(line)
            if not p.is_absolute():
                p = source.parent / p
            paths.append(str(p))
        return paths

    def _extract_taxa(self, path):
        """Extract taxon names from a tree or FASTA file."""
        if not Path(path).exists():
            raise PhykitUserError(
                [f"{path} corresponds to no such file or directory."],
                code=2,
            )

        if self.format == "fasta":
            from Bio import SeqIO
            try:
                records = list(SeqIO.parse(path, "fasta"))
                return [r.id for r in records]
            except Exception:
                raise PhykitUserError(
                    [f"Could not parse FASTA file: {path}"], code=2
                )
        else:
            from Bio import Phylo
            try:
                tree = Phylo.read(path, "newick")
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
        from Bio import SeqIO

        records = list(SeqIO.parse(input_path, "fasta"))
        kept = [r for r in records if r.id in kept_taxa]
        removed = len(records) - len(kept)

        with open(output_path, "w") as f:
            SeqIO.write(kept, f, "fasta")

        return len(kept), removed

    def _filter_tree(self, input_path, output_path, kept_taxa):
        """Filter a tree to keep only tips in kept_taxa."""
        from Bio import Phylo
        import copy

        tree = Phylo.read(input_path, "newick")
        tips = [t.name for t in tree.get_terminals()]
        tips_to_prune = [t for t in tips if t not in kept_taxa]

        if tips_to_prune:
            tree_copy = copy.deepcopy(tree)
            for tip_name in tips_to_prune:
                try:
                    tree_copy.prune(tip_name)
                except Exception:
                    pass
            Phylo.write(tree_copy, output_path, "newick")
        else:
            Phylo.write(tree, output_path, "newick")

        return len(tips) - len(tips_to_prune), len(tips_to_prune)

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
            print(f"Occupancy filter (threshold: {threshold_str})")
            print(f"Input files: {total_files}")
            print(f"Total taxa: {len(occupancy)}")
            print(f"Taxa kept: {len(kept_taxa)}")
            print(f"Taxa removed: {len(removed_taxa)}")
            print(f"Files modified: {files_modified}")
            print()

            if removed_taxa:
                print("Removed taxa (occupancy):")
                for taxon in sorted(removed_taxa):
                    print(f"  {taxon}: {occupancy[taxon]}/{total_files}")
                print()

            print("Occupancy per taxon:")
            for taxon in sorted(occupancy.keys()):
                status = "kept" if taxon in kept_taxa else "REMOVED"
                print(f"  {taxon}: {occupancy[taxon]}/{total_files} ({status})")
            print()

            print("Output files:")
            for in_path, out_path in zip(file_paths, output_paths):
                n_removed = len(taxa_removed_per_file.get(in_path, []))
                print(f"  {out_path}" + (
                    f" ({n_removed} taxa removed)" if n_removed else ""
                ))

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
                for taxon in sorted(occupancy.keys())
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
