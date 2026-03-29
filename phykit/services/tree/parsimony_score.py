"""
Maximum parsimony score of a tree given an alignment.

Computes the Fitch (1971) parsimony score — the minimum number of
character state changes required to explain the alignment on the
given tree topology. Each site is scored independently and the
total is summed.

Cross-validated against R's phangorn::parsimony().
"""
import copy
from typing import Dict, List

from Bio import SeqIO

from .base import Tree
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class ParsimonyScore(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.alignment_path = parsed["alignment_path"]
        self.verbose = parsed["verbose"]
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file()
        tree = copy.deepcopy(tree)
        self.validate_tree(tree, min_tips=3, context="parsimony score")
        self._resolve_polytomies(tree)

        sequences = self._parse_alignment(self.alignment_path)

        # Prune to shared taxa
        tree_tips = set(t.name for t in tree.get_terminals())
        seq_taxa = set(sequences.keys())
        shared = tree_tips & seq_taxa
        if len(shared) < 3:
            raise PhykitUserError(
                [
                    f"Only {len(shared)} shared taxa between tree and alignment.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        tips_to_prune = list(tree_tips - shared)
        if tips_to_prune:
            tree = self.prune_tree_using_taxa_list(tree, tips_to_prune)
            self._resolve_polytomies(tree)

        sequences = {t: sequences[t] for t in shared}
        aln_length = len(next(iter(sequences.values())))

        total_score, per_site = self._fitch_parsimony(tree, sequences, aln_length)

        if self.json_output:
            payload = {
                "parsimony_score": total_score,
                "alignment_length": aln_length,
                "n_taxa": len(shared),
            }
            if self.verbose:
                payload["per_site_scores"] = per_site
            print_json(payload)
        else:
            print(total_score)
            if self.verbose:
                print()
                for i, score in enumerate(per_site, 1):
                    print(f"{i}\t{score}")

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            alignment_path=args.alignment,
            verbose=getattr(args, "verbose", False),
            json_output=getattr(args, "json", False),
        )

    def _resolve_polytomies(self, tree) -> None:
        from Bio.Phylo import Newick
        for clade in tree.find_clades(order="postorder"):
            while len(clade.clades) > 2:
                child1 = clade.clades.pop()
                child2 = clade.clades.pop()
                new_internal = Newick.Clade(branch_length=0.0)
                new_internal.clades = [child1, child2]
                clade.clades.append(new_internal)

    def _parse_alignment(self, path: str) -> Dict[str, str]:
        try:
            records = list(SeqIO.parse(path, "fasta"))
        except Exception:
            raise PhykitUserError(
                [f"Could not parse alignment from {path}."], code=2
            )
        if not records:
            raise PhykitUserError(
                ["Alignment file is empty."], code=2
            )
        sequences = {r.id: str(r.seq).upper() for r in records}

        # Validate same length
        lengths = set(len(s) for s in sequences.values())
        if len(lengths) > 1:
            raise PhykitUserError(
                ["Sequences have different lengths. Provide an aligned FASTA file."],
                code=2,
            )
        return sequences

    def _fitch_parsimony(
        self, tree, sequences: Dict[str, str], aln_length: int
    ) -> tuple:
        """Compute Fitch parsimony score.

        For each site, performs a postorder (downpass) traversal:
        - Terminal nodes: state set = {observed character}
        - Internal nodes: if children share states, intersection;
          otherwise union (and add 1 to score)

        Gap characters (-) and ambiguous characters (N, X, ?) are
        treated as wildcards (matching any state).
        """
        wildcard_chars = {"-", "N", "X", "?", "n", "x"}
        all_states = {"A", "C", "G", "T"}
        per_site = []
        total_score = 0

        for site_idx in range(aln_length):
            site_score = 0
            node_states = {}

            for clade in tree.find_clades(order="postorder"):
                if clade.is_terminal():
                    char = sequences[clade.name][site_idx]
                    if char in wildcard_chars:
                        node_states[id(clade)] = set(all_states)
                    else:
                        node_states[id(clade)] = {char}
                else:
                    if len(clade.clades) != 2:
                        continue
                    left, right = clade.clades
                    left_states = node_states.get(id(left), set(all_states))
                    right_states = node_states.get(id(right), set(all_states))

                    intersection = left_states & right_states
                    if intersection:
                        node_states[id(clade)] = intersection
                    else:
                        node_states[id(clade)] = left_states | right_states
                        site_score += 1

            per_site.append(site_score)
            total_score += site_score

        return total_score, per_site
