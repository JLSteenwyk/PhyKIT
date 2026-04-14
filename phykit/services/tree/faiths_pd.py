"""
Faith's phylogenetic diversity (PD).

Given a tree and a community (list of tip labels), sum the branch
lengths of the minimum subtree connecting the community. When
``include_root`` is True (default), the sum includes the path from
the community's MRCA up to the tree root, matching Faith (1992) and
``picante::pd(..., include.root = TRUE)``.
"""
from typing import Dict, List, Tuple

from .base import Tree
from ...errors import PhykitUserError
from ...helpers.files import read_single_column_file_to_list
from ...helpers.json_output import print_json


class FaithsPD(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.taxa_file = parsed["taxa_file"]
        self.include_root = parsed["include_root"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            taxa_file=args.taxa,
            include_root=not getattr(args, "exclude_root", False),
            json_output=getattr(args, "json", False),
        )

    def run(self) -> None:
        tree = self.read_tree_file()
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
    def _load_taxa(taxa_file: str) -> List[str]:
        raw = read_single_column_file_to_list(taxa_file)
        seen = set()
        taxa: List[str] = []
        for name in raw:
            if not name or name in seen:
                continue
            seen.add(name)
            taxa.append(name)
        if not taxa:
            raise PhykitUserError(
                [f"No taxa found in {taxa_file}."], code=2,
            )
        return taxa

    def calculate_faiths_pd(
        self, tree, taxa: List[str], include_root: bool = True,
    ) -> Tuple[float, int]:
        """Compute Faith's PD for ``taxa`` on ``tree``.

        Returns (pd, n_taxa). n_taxa is the deduplicated community size.
        """
        self.validate_tree(tree, min_tips=2, require_branch_lengths=True,
                           context="Faith's PD")

        seen: set = set()
        deduped: List[str] = []
        for name in taxa:
            if name and name not in seen:
                seen.add(name)
                deduped.append(name)
        taxa = deduped
        if not taxa:
            raise PhykitUserError(
                ["Community must contain at least one taxon."], code=2,
            )

        tip_map = {t.name: t for t in tree.get_terminals()}
        missing = [name for name in taxa if name not in tip_map]
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

        community_tips = [tip_map[name] for name in taxa]

        if len(community_tips) == 1:
            # Single-tip community: PD with include_root=True is the path
            # length from the root to that tip; with include_root=False it
            # is 0 (no induced subtree). picante returns NA for the latter;
            # we return 0 for programmatic convenience.
            if not include_root:
                return 0.0, 1
            path = tree.get_path(community_tips[0])
            return (
                sum((c.branch_length or 0.0) for c in path),
                1,
            )

        if include_root:
            start_clade = tree.root
        else:
            mrca = tree.common_ancestor(community_tips)
            # If the MRCA is the root, include.root=False is equivalent
            # to include.root=True because there is no branch leading to
            # the MRCA to subtract (matches picante).
            start_clade = mrca

        total = 0.0
        seen_ids = set()
        for tip in community_tips:
            # Skip clades at or above start_clade. For include_root=True,
            # start_clade is the root, which is never in get_path.
            path = tree.get_path(tip)
            if start_clade is not tree.root:
                try:
                    idx = path.index(start_clade)
                    path = path[idx + 1:]
                except ValueError:
                    # start_clade (MRCA) not on path - defensive; shouldn't happen
                    pass
            for clade in path:
                cid = id(clade)
                if cid in seen_ids:
                    continue
                seen_ids.add(cid)
                total += clade.branch_length or 0.0

        return total, len(community_tips)
