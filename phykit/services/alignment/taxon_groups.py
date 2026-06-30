from collections import defaultdict
import os
from pathlib import Path

from ._fasta import read_fasta_first_tokens
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
        sorted_groups = sorted(groups.items(), key=lambda x: -len(x[1]))

        # 4. Output
        if self.json_output:
            result = {
                "total_files": len(file_paths),
                "total_groups": len(sorted_groups),
                "groups": [],
            }
            for i, (taxa_set, files) in enumerate(sorted_groups, 1):
                result["groups"].append({
                    "group": i,
                    "n_files": len(files),
                    "n_taxa": len(taxa_set),
                    "taxa": sorted(taxa_set),
                    "files": sorted(files),
                })
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
