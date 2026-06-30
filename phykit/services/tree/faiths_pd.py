"""
Faith's phylogenetic diversity (PD).

Given a tree and a community (list of tip labels), sum the branch
lengths of the minimum subtree connecting the community. When
``include_root`` is True (default), the sum includes the path from
the community's MRCA up to the tree root, matching Faith (1992) and
``picante::pd(..., include.root = TRUE)``.
"""
from __future__ import annotations

from .base import Tree
from ...errors import PhykitUserError


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def read_single_column_file_to_list(*args, **kwargs):
    from ...helpers.files import (
        read_single_column_file_to_list as _read_single_column_file_to_list,
    )

    return _read_single_column_file_to_list(*args, **kwargs)


class FaithsPD(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.taxa_file = parsed["taxa_file"]
        self.include_root = parsed["include_root"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> dict:
        return dict(
            tree_file_path=args.tree,
            taxa_file=args.taxa,
            include_root=not getattr(args, "exclude_root", False),
            json_output=getattr(args, "json", False),
        )

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        taxa = self._load_taxa(self.taxa_file)
        pd_value, n_tips = self.calculate_faiths_pd(
            tree, taxa, include_root=self.include_root
        )
        pd_rounded = round(pd_value, 4)

        if self.json_output:
            print_json(
                dict(
                    faiths_pd=pd_rounded,
                    n_taxa=n_tips,
                    include_root=self.include_root,
                )
            )
            return
        print(pd_rounded)

    @staticmethod
    def _load_taxa(taxa_file: str) -> list[str]:
        raw = read_single_column_file_to_list(taxa_file)
        taxa = [name for name in dict.fromkeys(raw) if name]
        if not taxa:
            raise PhykitUserError(
                [f"No taxa found in {taxa_file}."], code=2,
            )
        return taxa

    def calculate_faiths_pd(
        self, tree, taxa: list[str], include_root: bool = True,
    ) -> tuple[float, int]:
        """Compute Faith's PD for ``taxa`` on ``tree``.

        Returns (pd, n_taxa). n_taxa is the deduplicated community size.
        """
        self.validate_tree(tree, min_tips=2, require_branch_lengths=True,
                           context="Faith's PD")

        taxa = [name for name in dict.fromkeys(taxa) if name]
        if not taxa:
            raise PhykitUserError(
                ["Community must contain at least one taxon."], code=2,
            )

        selected_names = set(taxa)
        clades = []
        parent_by_id = {}
        tip_names = set()
        tip_by_name = {}
        total_branch_length = 0.0
        stack = [(tree.root, None)]
        pop = stack.pop
        append = stack.append
        while stack:
            clade, parent = pop()
            clades.append(clade)
            if parent is not None:
                parent_by_id[id(clade)] = parent
                total_branch_length += clade.branch_length or 0.0
            children = clade.clades
            child_count = len(children)
            if child_count == 2:
                append((children[1], clade))
                append((children[0], clade))
            elif child_count:
                for idx in range(child_count - 1, -1, -1):
                    append((children[idx], clade))
            else:
                tip_names.add(clade.name)
                if clade.name in selected_names:
                    tip_by_name[clade.name] = clade

        missing = [name for name in taxa if name not in tip_names]
        if missing:
            sample = ", ".join(sorted(missing)[:5])
            suffix = f" ... ({len(missing)} total)" if len(missing) > 5 else ""
            raise PhykitUserError(
                [
                    "Taxa not found in tree:",
                    sample + suffix,
                ],
                code=2,
            )
        if len(selected_names) == len(tip_names):
            return total_branch_length, len(taxa)
        if include_root:
            selected_branch_ids = set()
            total = 0.0
            for name in taxa:
                clade = tip_by_name[name]
                while clade is not tree.root:
                    cid = id(clade)
                    if cid in selected_branch_ids:
                        break
                    selected_branch_ids.add(cid)
                    total += clade.branch_length or 0.0
                    clade = parent_by_id[cid]
            return total, len(taxa)

        selected_counts = {}
        for clade in reversed(clades):
            children = clade.clades
            child_count = len(children)
            if child_count == 2:
                selected_counts[id(clade)] = (
                    selected_counts[id(children[0])]
                    + selected_counts[id(children[1])]
                )
            elif child_count:
                count = 0
                for child in children:
                    count += selected_counts[id(child)]
                selected_counts[id(clade)] = count
            else:
                selected_counts[id(clade)] = 1 if clade.name in selected_names else 0

        if len(taxa) == 1:
            # Single-tip community: PD with include_root=True is the path
            # length from the root to that tip; with include_root=False it
            # is 0 (no induced subtree). picante returns NA for the latter;
            # we return 0 for programmatic convenience.
            if not include_root:
                return 0.0, 1
            path = [
                clade for clade in clades
                if clade is not tree.root and selected_counts[id(clade)] > 0
            ]
            return (
                sum((c.branch_length or 0.0) for c in path),
                1,
            )

        excluded_ids = set()
        if include_root:
            pass
        else:
            community_size = len(taxa)
            mrca = next(
                (
                    clade for clade in reversed(clades)
                    if selected_counts[id(clade)] == community_size
                ),
                tree.root,
            )
            # If the MRCA is the root, include.root=False is equivalent
            # to include.root=True because there is no branch leading to
            # the MRCA to subtract (matches picante).
            if mrca is not tree.root:
                clade = mrca
                while True:
                    excluded_ids.add(id(clade))
                    if clade is tree.root:
                        break
                    clade = parent_by_id[id(clade)]

        total = 0.0
        for clade in clades:
            cid = id(clade)
            if clade is tree.root or cid in excluded_ids:
                continue
            if selected_counts[cid] > 0:
                total += clade.branch_length or 0.0

        return total, len(taxa)
