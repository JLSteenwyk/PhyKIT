"""
Maximum parsimony score of a tree given an alignment.

Computes the Fitch (1971) parsimony score — the minimum number of
character state changes required to explain the alignment on the
given tree topology. Each site is scored independently and the
total is summed.

Cross-validated against R's phangorn::parsimony().
"""
from __future__ import annotations

from ..alignment._fasta import read_fasta_first_token_upper
from .base import Tree
from ...errors import PhykitUserError


def _json_dumps(*args, **kwargs):
    import json

    return json.dumps(*args, **kwargs)


class _LazyNumpy:
    _module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = _np
            self._module = module
        value = getattr(module, name)
        setattr(self, name, value)
        return value


np = _LazyNumpy()


class ParsimonyScore(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(tree_file_path=parsed["tree_file_path"])
        self.alignment_path = parsed["alignment_path"]
        self.verbose = parsed["verbose"]
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        tree = self.read_tree_file_unmodified()
        self.validate_tree(tree, min_tips=3, context="parsimony score")
        copied_tree = False
        has_polytomies = self._has_polytomies(tree)
        if has_polytomies is not False:
            tree = self._fast_copy(tree)
            copied_tree = True
            self._resolve_polytomies(tree)

        sequences = self._parse_alignment(self.alignment_path)

        shared_count, tips_to_prune, sequences = self._shared_alignment_taxa_setup(
            self.get_tip_names_from_tree(tree),
            sequences,
        )
        if shared_count < 3:
            raise PhykitUserError(
                [
                    f"Only {shared_count} shared taxa between tree and alignment.",
                    "At least 3 shared taxa are required.",
                ],
                code=2,
            )

        if tips_to_prune:
            if not copied_tree:
                tree = self._fast_copy(tree)
                copied_tree = True
            tree = self.prune_tree_using_taxa_list(tree, tips_to_prune)
            self._resolve_polytomies(tree)

        aln_length = len(next(iter(sequences.values())))

        total_score, per_site = self._fitch_parsimony(
            tree,
            sequences,
            aln_length,
            return_per_site=self.verbose,
        )

        if self.json_output:
            payload = {
                "parsimony_score": total_score,
                "alignment_length": aln_length,
                "n_taxa": shared_count,
            }
            if self.verbose:
                payload["per_site_scores"] = per_site
            try:
                print(_json_dumps(payload, sort_keys=True))
            except BrokenPipeError:
                pass
        else:
            if self.verbose:
                lines = [str(total_score), ""]
                lines.extend(
                    f"{i}\t{score}"
                    for i, score in enumerate(per_site, 1)
                )
                print("\n".join(lines))
            else:
                print(total_score)

    @staticmethod
    def _shared_alignment_taxa_setup(tree_tip_names: list[str], sequences: dict):
        seq_count = len(sequences)
        if seq_count >= 3 and len(tree_tip_names) == seq_count:
            if tree_tip_names[0] == next(iter(sequences)):
                if tree_tip_names[-1] == next(reversed(sequences)):
                    seq_names = list(sequences)
                    if seq_names == tree_tip_names:
                        return seq_count, [], sequences

        tree_tips = set(tree_tip_names)
        seq_taxa = set(sequences.keys())
        shared = tree_tips & seq_taxa
        if len(shared) < 3:
            return len(shared), [], sequences
        tips_to_prune = list(tree_tips - shared)
        filtered_sequences = {taxon: sequences[taxon] for taxon in shared}
        return len(shared), tips_to_prune, filtered_sequences

    @staticmethod
    def _has_polytomies(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        stack = [root]
        pop = stack.pop
        extend = stack.extend
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if len(children) > 2:
                    return True
                if children:
                    extend(children)
        except AttributeError:
            return None
        return False

    def process_args(self, args) -> dict:
        return dict(
            tree_file_path=args.tree,
            alignment_path=args.alignment,
            verbose=getattr(args, "verbose", False),
            json_output=getattr(args, "json", False),
        )

    def _resolve_polytomies(self, tree) -> None:
        if self._resolve_standard_tree_polytomies(tree):
            return

        newick_clade = None
        for clade in tree.find_clades(order="postorder"):
            while len(clade.clades) > 2:
                if newick_clade is None:
                    from Bio.Phylo import Newick

                    newick_clade = Newick.Clade
                child1 = clade.clades.pop()
                child2 = clade.clades.pop()
                new_internal = newick_clade(branch_length=0.0)
                new_internal.clades = [child1, child2]
                clade.clades.append(new_internal)

    @staticmethod
    def _postorder_clades_fast(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        clades = []
        stack = [root]
        try:
            while stack:
                clade = stack.pop()
                clades.append(clade)
                stack.extend(clade.clades)
        except AttributeError:
            return None
        clades.reverse()
        return clades

    def _resolve_standard_tree_polytomies(self, tree) -> bool:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return False

        stack = [root]
        newick_clade = None
        try:
            pop = stack.pop
            extend = stack.extend
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    extend(children)
                while len(children) > 2:
                    if newick_clade is None:
                        from Bio.Phylo import Newick

                        newick_clade = Newick.Clade
                    child1 = children.pop()
                    child2 = children.pop()
                    new_internal = newick_clade(branch_length=0.0)
                    new_internal.clades = [child1, child2]
                    children.append(new_internal)
        except AttributeError:
            return False
        return True

    def _parse_alignment(self, path: str) -> dict[str, str]:
        try:
            sequences = read_fasta_first_token_upper(path)
        except Exception:
            raise PhykitUserError(
                [f"Could not parse alignment from {path}."], code=2
            )
        if not sequences:
            raise PhykitUserError(
                ["Alignment file is empty."], code=2
            )

        seq_iter = iter(sequences.values())
        first_length = len(next(seq_iter))
        for seq in seq_iter:
            if len(seq) != first_length:
                raise PhykitUserError(
                    ["Sequences have different lengths. Provide an aligned FASTA file."],
                    code=2,
                )
        return sequences

    def _fitch_parsimony(
        self,
        tree,
        sequences: dict[str, str],
        aln_length: int,
        return_per_site: bool = True,
    ) -> tuple:
        """Compute Fitch parsimony score.

        For each site, performs a postorder (downpass) traversal:
        - Terminal nodes: state set = {observed character}
        - Internal nodes: if children share states, intersection;
          otherwise union (and add 1 to score)

        Gap characters (-) and ambiguous characters (N, X, ?) are
        treated as wildcards (matching any state).
        """
        if self._sequences_are_identical(sequences) and (
            self._tree_terminals_have_sequences(tree, sequences)
        ):
            per_site = [0] * aln_length if return_per_site else []
            return 0, per_site

        if aln_length <= 16:
            return self._fitch_parsimony_small(
                tree,
                sequences,
                aln_length,
                return_per_site=return_per_site,
            )

        wildcard_chars = {"-", "N", "X", "?", "n", "x"}
        all_states = {"A", "C", "G", "T"}

        observed_states = {
            char
            for seq in sequences.values()
            for char in seq
            if char not in wildcard_chars
        }
        state_symbols = sorted(all_states | observed_states)
        if len(state_symbols) > 63:
            return self._fitch_parsimony_sets(
                tree,
                sequences,
                aln_length,
                return_per_site=return_per_site,
            )

        state_masks = {
            char: np.uint64(1 << idx)
            for idx, char in enumerate(state_symbols)
        }
        wildcard_mask = np.uint64(0)
        for char in all_states:
            wildcard_mask |= state_masks[char]

        state_lookup = None
        max_state_ord = max((ord(char) for char in state_symbols), default=0)
        if max_state_ord <= 255:
            state_lookup = np.empty(256, dtype=np.uint64)
            state_lookup.fill(wildcard_mask)
            for char, mask in state_masks.items():
                state_lookup[ord(char)] = mask
            for char in wildcard_chars:
                state_lookup[ord(char)] = wildcard_mask

        node_states = {}
        if return_per_site:
            site_scores = np.zeros(aln_length, dtype=np.int64)
        else:
            total_score = 0
        default_states = np.full(aln_length, wildcard_mask, dtype=np.uint64)
        sequence_state_cache = {}

        clades = self._postorder_clades_fast(tree)
        if clades is None:
            clades = tree.find_clades(order="postorder")

        for clade in clades:
            if not clade.clades:
                seq = sequences[clade.name]
                states = sequence_state_cache.get(seq)
                if states is None:
                    if state_lookup is not None:
                        states = state_lookup[
                            np.frombuffer(seq.encode("latin-1"), dtype=np.uint8)
                        ]
                    else:
                        states = np.fromiter(
                            (
                                wildcard_mask
                                if char in wildcard_chars
                                else state_masks[char]
                                for char in seq
                            ),
                            dtype=np.uint64,
                            count=aln_length,
                        )
                    sequence_state_cache[seq] = states
                node_states[id(clade)] = states
            else:
                if len(clade.clades) != 2:
                    continue
                left, right = clade.clades
                left_states = node_states.get(id(left), default_states)
                right_states = node_states.get(id(right), default_states)

                intersection = left_states & right_states
                changed = intersection == 0
                if return_per_site:
                    site_scores += changed
                else:
                    total_score += int(np.count_nonzero(changed))
                node_states[id(clade)] = np.where(
                    changed,
                    left_states | right_states,
                    intersection,
                )

        if return_per_site:
            total_score = int(site_scores.sum())
            per_site = site_scores.astype(int).tolist()
        else:
            per_site = []
        return total_score, per_site

    @staticmethod
    def _sequences_are_identical(sequences: dict[str, str]) -> bool:
        values = iter(sequences.values())
        try:
            first_sequence = next(values)
        except StopIteration:
            return False
        return all(sequence == first_sequence for sequence in values)

    @staticmethod
    def _tree_terminals_have_sequences(tree, sequences: dict[str, str]) -> bool:
        try:
            root = tree.root
            root.clades
        except AttributeError:
            for clade in tree.find_clades(order="postorder"):
                if not clade.clades and clade.name not in sequences:
                    return False
            return True

        stack = [root]
        pop = stack.pop
        extend = stack.extend
        try:
            while stack:
                clade = pop()
                children = clade.clades
                if children:
                    extend(children)
                elif clade.name not in sequences:
                    return False
        except AttributeError:
            for clade in tree.find_clades(order="postorder"):
                if not clade.clades and clade.name not in sequences:
                    return False
            return True
        return True

    def _fitch_parsimony_small(
        self,
        tree,
        sequences: dict[str, str],
        aln_length: int,
        return_per_site: bool = True,
    ) -> tuple:
        wildcard_chars = {"-", "N", "X", "?", "n", "x"}
        all_states = {"A", "C", "G", "T"}

        observed_states = {
            char
            for seq in sequences.values()
            for char in seq
            if char not in wildcard_chars
        }
        state_symbols = sorted(all_states | observed_states)
        if len(state_symbols) > 63:
            return self._fitch_parsimony_sets(
                tree,
                sequences,
                aln_length,
                return_per_site=return_per_site,
            )

        state_masks = {
            char: 1 << idx
            for idx, char in enumerate(state_symbols)
        }
        wildcard_mask = 0
        for char in all_states:
            wildcard_mask |= state_masks[char]

        state_lookup = {char: mask for char, mask in state_masks.items()}
        for char in wildcard_chars:
            state_lookup[char] = wildcard_mask

        node_states = {}
        if return_per_site:
            site_scores = [0] * aln_length
        else:
            total_score = 0
        default_states = [wildcard_mask] * aln_length
        sequence_states = {}
        sequence_state_cache = {}
        for taxon, seq in sequences.items():
            states = sequence_state_cache.get(seq)
            if states is None:
                states = [state_lookup[char] for char in seq]
                sequence_state_cache[seq] = states
            sequence_states[taxon] = states

        clades = self._postorder_clades_fast(tree)
        if clades is None:
            clades = tree.find_clades(order="postorder")

        for clade in clades:
            children = clade.clades
            if not children:
                node_states[id(clade)] = sequence_states[clade.name]
                continue

            if len(children) != 2:
                continue

            left, right = children
            left_states = node_states.get(id(left), default_states)
            right_states = node_states.get(id(right), default_states)

            states = []
            for idx in range(aln_length):
                left_state = left_states[idx]
                right_state = right_states[idx]
                intersection = left_state & right_state
                if intersection:
                    states.append(intersection)
                else:
                    states.append(left_state | right_state)
                    if return_per_site:
                        site_scores[idx] += 1
                    else:
                        total_score += 1

            node_states[id(clade)] = states

        if return_per_site:
            total_score = sum(site_scores)
            per_site = site_scores
        else:
            per_site = []
        return total_score, per_site

    def _fitch_parsimony_sets(
        self,
        tree,
        sequences: dict[str, str],
        aln_length: int,
        return_per_site: bool = True,
    ) -> tuple:
        wildcard_chars = {"-", "N", "X", "?", "n", "x"}
        all_states = {"A", "C", "G", "T"}
        default_states = set(all_states)
        clades = self._postorder_clades_fast(tree)
        if clades is None:
            clades = list(tree.find_clades(order="postorder"))

        per_site = [] if return_per_site else None
        total_score = 0
        id_ = id

        for site_idx in range(aln_length):
            site_score = 0
            node_states = {}
            node_states_get = node_states.get

            for clade in clades:
                children = clade.clades
                if not children:
                    char = sequences[clade.name][site_idx]
                    if char in wildcard_chars:
                        node_states[id_(clade)] = default_states
                    else:
                        node_states[id_(clade)] = {char}
                else:
                    if len(children) != 2:
                        continue
                    left, right = children
                    left_states = node_states_get(id_(left), default_states)
                    right_states = node_states_get(id_(right), default_states)

                    intersection = left_states & right_states
                    if intersection:
                        node_states[id_(clade)] = intersection
                    else:
                        node_states[id_(clade)] = left_states | right_states
                        site_score += 1

            if return_per_site:
                per_site.append(site_score)
            total_score += site_score

        return total_score, per_site if return_per_site else []
