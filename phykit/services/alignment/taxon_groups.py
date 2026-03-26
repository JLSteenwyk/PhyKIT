from collections import defaultdict
from pathlib import Path

from ...errors import PhykitUserError
from ...helpers.json_output import print_json


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

        # 2. Extract taxa from each file
        file_taxa = {}
        for path in file_paths:
            taxa = self._extract_taxa(path)
            file_taxa[path] = frozenset(taxa)

        # 3. Group by taxon set
        groups = defaultdict(list)
        for path, taxa in file_taxa.items():
            groups[taxa].append(path)

        # 4. Sort groups by size (largest first)
        sorted_groups = sorted(groups.items(), key=lambda x: -len(x[1]))

        # 5. Output
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
            print(f"Total files: {len(file_paths)}")
            print(f"Total groups: {len(sorted_groups)}")
            print()
            for i, (taxa_set, files) in enumerate(sorted_groups, 1):
                n_files = len(files)
                n_taxa = len(taxa_set)
                files_str = ", ".join(sorted(files))
                taxa_str = ", ".join(sorted(taxa_set))
                print(f"Group {i} ({n_taxa} taxa, {n_files} files):")
                print(f"  Files: {files_str}")
                print(f"  Taxa: {taxa_str}")
                print()

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
            # Resolve relative paths
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
