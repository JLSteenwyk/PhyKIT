"""Patterson's D-statistic (ABBA-BABA test) for detecting introgression.

Supports two modes:
1) Site patterns from a FASTA alignment (-a)
2) Quartet topologies from gene trees (-g)
"""

from __future__ import annotations

import math

from ._fasta import read_fasta_first_token_upper
from .base import Alignment
from ...errors import PhykitUserError


class _LazyNumpy:
    def __init__(self):
        self._module = None

    def __getattr__(self, name):
        module = self._module
        if module is None:
            import numpy as _np

            module = _np
            self._module = module

        attr = getattr(module, name)
        setattr(self, name, attr)
        return attr


np = _LazyNumpy()
_SKIP_CODES = (ord("-"), ord("N"), ord("?"), ord("X"), ord("n"), ord("x"))
_SKIP_BYTES = b"-N?Xnx"
_SCALAR_SKIP_CHARS = "-N?Xnx"
_SKIP_SCAN_BYTES = 4096


def _same_sequence(seq_a: str, seq_b: str) -> bool:
    if seq_a is seq_b:
        return True
    if len(seq_a) != len(seq_b):
        return False
    if not seq_a:
        return True
    if seq_a[0] != seq_b[0] or seq_a[-1] != seq_b[-1]:
        return False
    seq_len = len(seq_a)
    if seq_len > 3 and (
        seq_a[1] != seq_b[1]
        or seq_a[seq_len // 2] != seq_b[seq_len // 2]
        or seq_a[-2] != seq_b[-2]
    ):
        return False
    if len(seq_a) > _SKIP_SCAN_BYTES:
        return (
            seq_a[:_SKIP_SCAN_BYTES] == seq_b[:_SKIP_SCAN_BYTES]
            and seq_a[-_SKIP_SCAN_BYTES:] == seq_b[-_SKIP_SCAN_BYTES:]
            and seq_a == seq_b
        )
    return seq_a == seq_b


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def _sum_squared_deviations(values, center) -> float:
    centered = values - center
    return float(np.dot(centered, centered))


def _normal_two_tailed_p_value(z_score: float) -> float:
    return math.erfc(abs(z_score) / math.sqrt(2.0))


def _chi2_sf_df1(chi2_stat: float) -> float:
    return math.erfc(math.sqrt(chi2_stat / 2.0))


def _discordant_chi2_stat(abba_count: int, baba_count: int) -> float:
    n_informative = abba_count + baba_count
    if n_informative == 0:
        return 0.0
    diff = abba_count - baba_count
    return (diff * diff) / n_informative


class Dstatistic(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_path"])
        self.gene_trees_path = parsed["gene_trees_path"]
        self.p1 = parsed["p1"]
        self.p2 = parsed["p2"]
        self.p3 = parsed["p3"]
        self.outgroup = parsed["outgroup"]
        self.block_size = parsed["block_size"]
        self.support_threshold = parsed["support_threshold"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> dict[str, object]:
        aln = getattr(args, "alignment", None)
        gt = getattr(args, "gene_trees", None)
        if aln is None and gt is None:
            raise PhykitUserError(
                ["Either -a/--alignment or -g/--gene-trees is required."],
                code=2,
            )
        if aln is not None and gt is not None:
            raise PhykitUserError(
                ["-a/--alignment and -g/--gene-trees are mutually exclusive."],
                code=2,
            )
        return dict(
            alignment_path=aln,
            gene_trees_path=gt,
            p1=args.p1,
            p2=args.p2,
            p3=args.p3,
            outgroup=args.outgroup,
            block_size=getattr(args, "block_size", 100),
            support_threshold=getattr(args, "support", None),
            json_output=getattr(args, "json", False),
        )

    def run(self):
        if self.gene_trees_path:
            self._run_gene_tree_mode()
        else:
            self._run_alignment_mode()

    # ------------------------------------------------------------------
    # Gene tree mode
    # ------------------------------------------------------------------

    def _run_gene_tree_mode(self):
        """Count quartet topologies from gene trees."""
        gene_trees = self._parse_gene_trees(self.gene_trees_path)
        quartet = (self.p1, self.p2, self.p3, self.outgroup)

        # For each gene tree, determine the quartet topology
        # Species tree: (((P1,P2),P3),O)
        # Concordant: P1+P2 together → ((P1,P2),(P3,O))
        # ABBA: P2+P3 together → ((P2,P3),(P1,O))
        # BABA: P1+P3 together → ((P1,P3),(P2,O))
        concordant = 0
        abba_count = 0
        baba_count = 0
        unresolved = 0

        for gt in gene_trees:
            topo = self._get_quartet_topology(gt, quartet)
            if topo == "concordant":
                concordant += 1
            elif topo == "abba":
                abba_count += 1
            elif topo == "baba":
                baba_count += 1
            else:
                unresolved += 1

        n_informative = abba_count + baba_count
        n_total = len(gene_trees)

        # D-statistic
        if n_informative == 0:
            d_stat = 0.0
        else:
            d_stat = (abba_count - baba_count) / n_informative

        p_value = None
        chi2_stat = None
        if n_informative > 0:
            chi2_stat = _discordant_chi2_stat(abba_count, baba_count)
            p_value = _chi2_sf_df1(chi2_stat)

        # Output
        if self.json_output:
            payload = {
                "mode": "gene_trees",
                "p1": self.p1,
                "p2": self.p2,
                "p3": self.p3,
                "outgroup": self.outgroup,
                "n_gene_trees": n_total,
                "concordant": concordant,
                "abba_count": abba_count,
                "baba_count": baba_count,
                "unresolved": unresolved,
                "d_statistic": round(d_stat, 4),
                "support_threshold": self.support_threshold,
                "chi2_statistic": round(chi2_stat, 4) if chi2_stat is not None else None,
                "p_value": round(p_value, 6) if p_value is not None else None,
            }
            print_json(payload, sort_keys=False)
            return

        self._print_gene_tree_text_output(
            n_total,
            concordant,
            abba_count,
            baba_count,
            unresolved,
            d_stat,
            chi2_stat,
            p_value,
        )

    def _print_gene_tree_text_output(
        self,
        n_total,
        concordant,
        abba_count,
        baba_count,
        unresolved,
        d_stat,
        chi2_stat,
        p_value,
    ):
        threshold_line = (
            f"Support threshold: {self.support_threshold}\n"
            if self.support_threshold is not None
            else ""
        )
        if chi2_stat is not None:
            tail = (
                f"Chi-squared: {chi2_stat:.4f}\n"
                f"p-value: {p_value:.6f}\n"
                f"\n"
                f"Interpretation: {self._interpret(d_stat, p_value)}"
            )
        else:
            tail = "\nNo informative (discordant) gene trees found."

        try:
            print(
                f"Patterson's D-statistic (Gene Tree Mode)\n"
                f"=========================================\n"
                f"Topology: ((({self.p1}, {self.p2}), "
                f"{self.p3}), {self.outgroup})\n"
                f"P1: {self.p1}\n"
                f"P2: {self.p2}\n"
                f"P3: {self.p3}\n"
                f"Outgroup: {self.outgroup}\n"
                f"\n"
                f"Gene trees: {n_total}\n"
                f"{threshold_line}"
                f"Concordant ((P1,P2),P3): {concordant}\n"
                f"ABBA ((P2,P3),P1): {abba_count}\n"
                f"BABA ((P1,P3),P2): {baba_count}\n"
                f"Unresolved: {unresolved}\n"
                f"D-statistic: {d_stat:.4f}\n"
                f"{tail}"
            )
        except BrokenPipeError:
            pass

    def _parse_gene_trees(self, path: str) -> list:
        """Parse gene trees from a file (one Newick per line)."""
        from Bio import Phylo

        try:
            return list(Phylo.parse(path, "newick"))
        except Exception:
            raise PhykitUserError(
                [f"Could not parse gene trees from {path}."],
                code=2,
            )

    @staticmethod
    def _collect_clade_taxa(tree) -> dict[int, frozenset]:
        clade_taxa, _ = Dstatistic._collect_clade_taxa_and_nonterminals(tree)
        return clade_taxa

    @staticmethod
    def _collect_clade_taxa_and_nonterminals(
        tree,
    ) -> tuple[dict[int, frozenset], list]:
        direct_result = Dstatistic._collect_clade_taxa_and_nonterminals_direct(tree)
        if direct_result is not None:
            return direct_result

        clade_taxa: dict[int, frozenset] = {}
        for clade in tree.find_clades(order="postorder"):
            if clade.is_terminal():
                clade_taxa[id(clade)] = frozenset({clade.name})
            else:
                taxa = frozenset()
                for child in clade.clades:
                    taxa = taxa | clade_taxa.get(id(child), frozenset())
                clade_taxa[id(clade)] = taxa
        return clade_taxa, list(tree.get_nonterminals())

    @staticmethod
    def _collect_clade_taxa_and_nonterminals_direct(tree):
        try:
            root = tree.root
            root.clades
        except AttributeError:
            return None

        preorder = []
        nonterminals = []
        stack = [root]
        append = preorder.append
        append_nonterminal = nonterminals.append
        pop = stack.pop
        push = stack.append
        try:
            while stack:
                clade = pop()
                append(clade)
                children = clade.clades
                if children:
                    append_nonterminal(clade)
                    child_count = len(children)
                    if child_count == 2:
                        push(children[1])
                        push(children[0])
                    elif child_count == 1:
                        push(children[0])
                    else:
                        for idx in range(child_count - 1, -1, -1):
                            push(children[idx])
        except AttributeError:
            return None

        clade_taxa: dict[int, frozenset] = {}
        empty_taxa = frozenset()
        empty_taxa_union = empty_taxa.union
        for clade in reversed(preorder):
            children = clade.clades
            child_count = len(children)
            if child_count == 0:
                clade_taxa[id(clade)] = frozenset({clade.name})
            elif child_count == 2:
                clade_taxa[id(clade)] = (
                    clade_taxa[id(children[0])] | clade_taxa[id(children[1])]
                )
            else:
                clade_taxa[id(clade)] = empty_taxa_union(
                    *(clade_taxa.get(id(child), empty_taxa) for child in children)
                )

        return clade_taxa, nonterminals

    def _get_quartet_topology(self, tree, quartet) -> str:
        """Determine quartet topology from a (possibly multi-taxon) gene tree.

        If support_threshold is set, branches with support below the
        threshold are excluded (treated as collapsed/unresolved).

        Returns 'concordant', 'abba', 'baba', or 'unresolved'.
        """
        p1, p2, p3, outgroup = quartet

        clade_taxa, nonterminals = self._collect_clade_taxa_and_nonterminals(tree)
        tree_taxa = clade_taxa.get(id(tree.root), frozenset())

        # Check all four taxa are present
        if not all(t in tree_taxa for t in quartet):
            return "unresolved"

        # Extract bipartitions from the gene tree
        # Skip branches with support below threshold
        all_taxa = tree_taxa
        for clade in nonterminals:
            # Check support threshold
            if self.support_threshold is not None:
                support = clade.confidence
                if support is not None and support < self.support_threshold:
                    continue  # collapse this branch (skip its bipartition)

            tips = clade_taxa.get(id(clade), frozenset())
            if len(tips) <= 1 or tips == all_taxa:
                continue

            # Check which quartet topology the bipartition supports.
            has_p1 = p1 in tips
            has_p2 = p2 in tips
            has_p3 = p3 in tips
            has_outgroup = outgroup in tips
            if has_p1 + has_p2 + has_p3 + has_outgroup != 2:
                continue

            # Concordant: P1+P2 on one side
            if (has_p1 and has_p2) or (has_p3 and has_outgroup):
                return "concordant"
            # ABBA: P2+P3 on one side
            if (has_p2 and has_p3) or (has_p1 and has_outgroup):
                return "abba"
            # BABA: P1+P3 on one side
            if (has_p1 and has_p3) or (has_p2 and has_outgroup):
                return "baba"

        return "unresolved"

    # ------------------------------------------------------------------
    # Alignment mode
    # ------------------------------------------------------------------

    @staticmethod
    def _read_fasta(path: str) -> dict[str, str]:
        return read_fasta_first_token_upper(path)

    @staticmethod
    def _count_site_patterns(seq_p1, seq_p2, seq_p3, seq_outgroup, block_size):
        """Count ABBA/BABA site patterns and per-block counts."""
        if (
            len(seq_p1)
            == len(seq_p2)
            == len(seq_p3)
            == len(seq_outgroup)
            and (
                _same_sequence(seq_p1, seq_p2)
                or _same_sequence(seq_p3, seq_outgroup)
            )
        ):
            n_blocks = len(seq_p1) // block_size
            if n_blocks == 0:
                return 0, 0, np.empty(0), np.empty(0)
            return 0, 0, np.zeros(n_blocks), np.zeros(n_blocks)

        try:
            p1_bytes = seq_p1.encode("ascii")
            p2_bytes = seq_p2.encode("ascii")
            p3_bytes = seq_p3.encode("ascii")
            out_bytes = seq_outgroup.encode("ascii")
            p1 = np.frombuffer(p1_bytes, dtype=np.uint8)
            p2 = np.frombuffer(p2_bytes, dtype=np.uint8)
            p3 = np.frombuffer(p3_bytes, dtype=np.uint8)
            out = np.frombuffer(out_bytes, dtype=np.uint8)
        except UnicodeEncodeError:
            return Dstatistic._count_site_patterns_scalar(
                seq_p1,
                seq_p2,
                seq_p3,
                seq_outgroup,
                block_size,
            )

        has_skip_code = any(
            code in p1_bytes[:_SKIP_SCAN_BYTES]
            or code in p2_bytes[:_SKIP_SCAN_BYTES]
            or code in p3_bytes[:_SKIP_SCAN_BYTES]
            or code in out_bytes[:_SKIP_SCAN_BYTES]
            or code in p1_bytes[-_SKIP_SCAN_BYTES:]
            or code in p2_bytes[-_SKIP_SCAN_BYTES:]
            or code in p3_bytes[-_SKIP_SCAN_BYTES:]
            or code in out_bytes[-_SKIP_SCAN_BYTES:]
            for code in _SKIP_BYTES
        )
        if not has_skip_code:
            has_skip_code = any(
                code in p1_bytes
                or code in p2_bytes
                or code in p3_bytes
                or code in out_bytes
                for code in _SKIP_BYTES
            )

        if has_skip_code:
            valid = np.ones(len(p1), dtype=bool)
            for code in _SKIP_CODES:
                valid &= (p1 != code) & (p2 != code) & (p3 != code) & (out != code)

            if _same_sequence(seq_p2, seq_p3) or _same_sequence(
                seq_p1,
                seq_outgroup,
            ):
                abba_mask = valid & (p1 == out) & (p2 == p3) & (p2 != out)
                baba_mask = None
            elif _same_sequence(seq_p1, seq_p3) or _same_sequence(
                seq_p2,
                seq_outgroup,
            ):
                abba_mask = None
                baba_mask = valid & (p2 == out) & (p1 == p3) & (p1 != out)
            else:
                abba_mask = valid & (p1 == out) & (p2 == p3) & (p2 != out)
                baba_mask = valid & (p2 == out) & (p1 == p3) & (p1 != out)
        else:
            if _same_sequence(seq_p2, seq_p3) or _same_sequence(
                seq_p1,
                seq_outgroup,
            ):
                abba_mask = (p1 == out) & (p2 == p3) & (p2 != out)
                baba_mask = None
            elif _same_sequence(seq_p1, seq_p3) or _same_sequence(
                seq_p2,
                seq_outgroup,
            ):
                abba_mask = None
                baba_mask = (p2 == out) & (p1 == p3) & (p1 != out)
            else:
                abba_mask = (p1 == out) & (p2 == p3) & (p2 != out)
                baba_mask = (p2 == out) & (p1 == p3) & (p1 != out)

        abba_count = (
            int(np.count_nonzero(abba_mask))
            if abba_mask is not None
            else 0
        )
        baba_count = (
            int(np.count_nonzero(baba_mask))
            if baba_mask is not None
            else 0
        )

        n_blocks = len(p1) // block_size
        if n_blocks == 0:
            return (
                abba_count,
                baba_count,
                np.empty(0),
                np.empty(0),
            )
        block_abba = np.zeros(n_blocks)
        block_baba = np.zeros(n_blocks)
        n_block_sites = n_blocks * block_size
        if abba_mask is not None:
            block_abba = (
                abba_mask[:n_block_sites]
                .reshape(n_blocks, block_size)
                .sum(axis=1)
                .astype(float)
            )
        if baba_mask is not None:
            block_baba = (
                baba_mask[:n_block_sites]
                .reshape(n_blocks, block_size)
                .sum(axis=1)
                .astype(float)
            )

        return abba_count, baba_count, block_abba, block_baba

    @staticmethod
    def _count_site_patterns_scalar(seq_p1, seq_p2, seq_p3, seq_outgroup, block_size):
        skip_chars = _SCALAR_SKIP_CHARS
        aln_length = len(seq_p1)
        n_blocks = aln_length // block_size
        n_block_sites = n_blocks * block_size
        block_abba = [0.0] * n_blocks
        block_baba = [0.0] * n_blocks
        abba_count = 0
        baba_count = 0

        for site, (p1, p2, p3, o) in enumerate(
            zip(seq_p1, seq_p2, seq_p3, seq_outgroup)
        ):
            if (
                p1 in skip_chars
                or p2 in skip_chars
                or p3 in skip_chars
                or o in skip_chars
            ):
                continue

            if p1 == o and p2 != o and p3 != o and p2 == p3:
                abba_count += 1
                if site < n_block_sites:
                    block_abba[site // block_size] += 1.0
            elif p2 == o and p1 != o and p3 != o and p1 == p3:
                baba_count += 1
                if site < n_block_sites:
                    block_baba[site // block_size] += 1.0

        return (
            abba_count,
            baba_count,
            np.asarray(block_abba, dtype=float),
            np.asarray(block_baba, dtype=float),
        )

    @staticmethod
    def _jackknife_d_values(block_abba, block_baba) -> np.ndarray:
        n_blocks = len(block_abba)
        if n_blocks == 0:
            return np.zeros(0, dtype=float)

        total_abba = block_abba.sum()
        total_baba = block_baba.sum()
        loo_abba = total_abba - block_abba
        loo_baba = total_baba - block_baba
        denom = loo_abba + loo_baba

        jackknife_d = np.zeros(n_blocks, dtype=float)
        np.divide(
            loo_abba - loo_baba,
            denom,
            out=jackknife_d,
            where=denom > 0,
        )
        return jackknife_d

    def _run_alignment_mode(self):
        """Count ABBA/BABA site patterns from an alignment."""
        # Read alignment sequences
        sequences = self._read_fasta(self.alignment_file_path)

        # Validate taxa are present
        required = {"p1": self.p1, "p2": self.p2, "p3": self.p3, "outgroup": self.outgroup}
        for label, taxon in required.items():
            if taxon not in sequences:
                raise PhykitUserError(
                    [f"Taxon '{taxon}' ({label}) not found in alignment. "
                     f"Available taxa: {', '.join(sorted(sequences.keys()))}"],
                    code=2,
                )

        seq_p1 = sequences[self.p1]
        seq_p2 = sequences[self.p2]
        seq_p3 = sequences[self.p3]
        seq_outgroup = sequences[self.outgroup]

        # Validate equal lengths
        lengths = {len(seq_p1), len(seq_p2), len(seq_p3), len(seq_outgroup)}
        if len(lengths) != 1:
            raise PhykitUserError(
                ["Sequences have different lengths. All sequences must be aligned."],
                code=2,
            )

        aln_length = len(seq_p1)

        abba_count, baba_count, block_abba, block_baba = self._count_site_patterns(
            seq_p1,
            seq_p2,
            seq_p3,
            seq_outgroup,
            self.block_size,
        )

        informative_sites = abba_count + baba_count

        # Compute D-statistic
        if informative_sites == 0:
            d_stat = 0.0
        else:
            d_stat = (abba_count - baba_count) / informative_sites

        # Block jackknife for significance
        n_blocks = aln_length // self.block_size
        se = None
        z_score = None
        p_value = None

        if n_blocks >= 2:
            jackknife_d = self._jackknife_d_values(block_abba, block_baba)
            mean_d = jackknife_d.mean()
            sum_sq = _sum_squared_deviations(jackknife_d, mean_d)
            se = float(np.sqrt((n_blocks - 1) / n_blocks * sum_sq))

            if se > 0:
                z_score = d_stat / se
                p_value = _normal_two_tailed_p_value(z_score)
            else:
                z_score = float('inf') if d_stat != 0 else 0.0
                p_value = 0.0 if d_stat != 0 else 1.0

        # Output
        if self.json_output:
            payload = {
                "p1": self.p1,
                "p2": self.p2,
                "p3": self.p3,
                "outgroup": self.outgroup,
                "alignment_length": aln_length,
                "informative_sites": informative_sites,
                "abba_count": abba_count,
                "baba_count": baba_count,
                "d_statistic": round(d_stat, 4),
                "block_size": self.block_size,
                "n_blocks": n_blocks if n_blocks >= 2 else n_blocks,
                "standard_error": round(se, 4) if se is not None else None,
                "z_score": round(z_score, 2) if z_score is not None and z_score != float('inf') else z_score,
                "p_value": round(p_value, 6) if p_value is not None else None,
            }
            print_json(payload, sort_keys=False)
            return

        self._print_alignment_text_output(
            aln_length,
            informative_sites,
            abba_count,
            baba_count,
            d_stat,
            se,
            z_score,
            p_value,
        )

    def _print_alignment_text_output(
        self,
        aln_length,
        informative_sites,
        abba_count,
        baba_count,
        d_stat,
        se,
        z_score,
        p_value,
    ):
        if se is not None:
            z_text = "inf" if z_score == float('inf') else f"{z_score:.2f}"
            tail = (
                f"Block jackknife (block size: {self.block_size}):\n"
                f"  Standard error: {se:.4f}\n"
                f"  Z-score: {z_text}\n"
                f"  p-value: {p_value:.6f}\n"
                f"\n"
                f"Interpretation: {self._interpret(d_stat, p_value)}"
            )
        else:
            tail = "\nNot enough blocks for jackknife significance test."

        try:
            print(
                f"Patterson's D-statistic (ABBA-BABA Test)\n"
                f"=========================================\n"
                f"Topology: ((({self.p1}, {self.p2}), "
                f"{self.p3}), {self.outgroup})\n"
                f"P1: {self.p1}\n"
                f"P2: {self.p2}\n"
                f"P3: {self.p3}\n"
                f"Outgroup: {self.outgroup}\n"
                f"\n"
                f"Alignment length: {aln_length}\n"
                f"Informative sites: {informative_sites}\n"
                f"ABBA sites: {abba_count}\n"
                f"BABA sites: {baba_count}\n"
                f"D-statistic: {d_stat:.4f}\n"
                f"{tail}"
            )
        except BrokenPipeError:
            pass

    def _interpret(self, d_stat: float, p_value: float, alpha: float = 0.05) -> str:
        if p_value < alpha:
            if d_stat > 0:
                return (
                    f"Significant excess of ABBA patterns (p < {alpha}) "
                    f"suggests introgression between P2 ({self.p2}) and "
                    f"P3 ({self.p3}). Note: D cannot determine the "
                    f"direction of gene flow."
                )
            else:
                return (
                    f"Significant excess of BABA patterns (p < {alpha}) "
                    f"suggests introgression between P1 ({self.p1}) and "
                    f"P3 ({self.p3}). Note: D cannot determine the "
                    f"direction of gene flow."
                )
        return "No significant evidence of introgression (consistent with ILS)."
