from collections import defaultdict
import os

from ._fasta import read_fasta_first_tokens
from ...errors import PhykitUserError
from ._file_list import _normalize_list_path, read_file_list


_path_exists = os.path.exists


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def Path(*args, **kwargs):
    from pathlib import Path as _Path

    return _Path(*args, **kwargs)


class TaxonGroups:
    """Determine which files share the same set of taxa."""

    def __init__(self, args):
        parsed = self.process_args(args)
        self.list_path = parsed["list_path"]
        self.format = parsed["format"]
        self.json_output = parsed["json_output"]

    def process_args(self, args):
        return dict(
            list_path=args.list,
            format=getattr(args, "format", "trees"),
            json_output=getattr(args, "json", False),
        )

    def run(self):
        # 1. Read file list
        file_paths = self._read_file_list(self.list_path)

        if not file_paths:
            raise PhykitUserError(["No files found in list."], code=2)

        # 2. Extract taxa from each file and group by taxon set
        groups = defaultdict(list)
        for path in file_paths:
            taxa = self._extract_taxa(path)
            groups[frozenset(taxa)].append(path)

        # 3. Sort groups by size (largest first)
        sorted_groups = sorted(groups.items(), key=lambda x: len(x[1]), reverse=True)

        # 4. Output
        if self.json_output:
            result = {
                "total_files": len(file_paths),
                "total_groups": len(sorted_groups),
                "groups": [
                    {
                        "group": i,
                        "n_files": len(files),
                        "n_taxa": len(taxa_set),
                        "taxa": sorted(taxa_set),
                        "files": sorted(files),
                    }
                    for i, (taxa_set, files) in enumerate(sorted_groups, 1)
                ],
            }
            print_json(result)
        else:
            lines = [
                f"Total files: {len(file_paths)}",
                f"Total groups: {len(sorted_groups)}",
                "",
            ]
            for i, (taxa_set, files) in enumerate(sorted_groups, 1):
                n_files = len(files)
                n_taxa = len(taxa_set)
                files_str = ", ".join(sorted(files))
                taxa_str = ", ".join(sorted(taxa_set))
                lines.extend(
                    [
                        f"Group {i} ({n_taxa} taxa, {n_files} files):",
                        f"  Files: {files_str}",
                        f"  Taxa: {taxa_str}",
                        "",
                    ]
                )
            print("\n".join(lines))

    def _read_file_list(self, path):
        return read_file_list(
            path, path_factory=Path, normalize_path=_normalize_list_path,
        )

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
            from ..tree.base import Tree

            try:
                taxa = Tree._scan_simple_newick_tip_names(path)
                if taxa is not None:
                    return list(taxa)

                from Bio import Phylo

                tree = Phylo.read(path, "newick")
                taxa = Tree.calculate_terminal_names_fast(tree)
                if taxa is not None:
                    return taxa
                return [t.name for t in tree.get_terminals()]
            except Exception:
                raise PhykitUserError(
                    [f"Could not parse tree file: {path}"], code=2
                )
