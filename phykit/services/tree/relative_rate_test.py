import itertools
from pathlib import Path
from typing import Dict, List

from .base import Tree
from ...helpers.files import (
    get_alignment_and_format as get_alignment_and_format_helper
)
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


def _chi2_sf(x, df):
    """Chi-squared survival function (lazy scipy import)."""
    from scipy.stats import chi2
    return float(chi2.sf(x, df))


class RelativeRateTest(Tree):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(
            tree_file_path=parsed["tree_file_path"],
            alignment_file_path=parsed.get("alignment_file_path"),
            verbose=parsed["verbose"],
        )
        self.alignment_list = parsed.get("alignment_list")
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> Dict:
        return dict(
            tree_file_path=args.tree,
            alignment_file_path=getattr(args, "alignment", None),
            alignment_list=getattr(args, "alignment_list", None),
            verbose=getattr(args, "verbose", False),
            json_output=getattr(args, "json", False),
        )

    # ------------------------------------------------------------------
    # Core Tajima's test
    # ------------------------------------------------------------------

    @staticmethod
    def _tajima_test(seq_x: str, seq_y: str, seq_out: str) -> Dict:
        """Tajima's relative rate test (Equation 4).

        Counts unique substitutions on each ingroup lineage relative
        to the outgroup.  Returns m1, m2, chi-squared, and p-value.
        """
        m1 = 0  # unique changes on lineage x
        m2 = 0  # unique changes on lineage y
        skip_chars = {"-", "?", "N", "X", "n", "x"}

        for cx, cy, co in zip(seq_x, seq_y, seq_out):
            if cx in skip_chars or cy in skip_chars or co in skip_chars:
                continue
            x_diff = cx.upper() != co.upper()
            y_diff = cy.upper() != co.upper()
            if x_diff and not y_diff:
                m1 += 1
            elif y_diff and not x_diff:
                m2 += 1
            # Both differ or neither differs -> uninformative

        if m1 + m2 == 0:
            return {"m1": m1, "m2": m2, "chi2": 0.0, "p_value": 1.0}

        chi2 = (m1 - m2) ** 2 / (m1 + m2)
        p_value = _chi2_sf(chi2, df=1)

        return {"m1": m1, "m2": m2, "chi2": chi2, "p_value": p_value}

    # ------------------------------------------------------------------
    # Outgroup identification
    # ------------------------------------------------------------------

    @staticmethod
    def _identify_outgroup(tree) -> str:
        """Identify outgroup as the earliest-diverging taxon in a rooted tree.

        The outgroup is the single taxon on the smaller side of the root split.
        """
        root = tree.root
        children = root.clades
        if len(children) < 2:
            raise PhykitUserError(
                ["Tree does not appear to be rooted (root has < 2 children)."]
            )

        # Find the child clade with fewest terminals
        clade_sizes = []
        for child in children:
            terminals = list(child.get_terminals())
            clade_sizes.append((len(terminals), terminals, child))

        clade_sizes.sort(key=lambda x: x[0])
        smallest_clade = clade_sizes[0]

        if smallest_clade[0] != 1:
            raise PhykitUserError(
                [
                    "Cannot determine a single outgroup taxon from the rooted tree. "
                    "The smallest clade at the root contains multiple taxa. "
                    "Please ensure the tree is rooted with a single outgroup taxon."
                ]
            )

        return smallest_clade[1][0].name

    # ------------------------------------------------------------------
    # Multiple testing correction
    # ------------------------------------------------------------------

    @staticmethod
    def _bonferroni(p_values: List[float]) -> List[float]:
        n = len(p_values)
        return [min(p * n, 1.0) for p in p_values]

    @staticmethod
    def _fdr(p_values: List[float]) -> List[float]:
        """Benjamini-Hochberg FDR correction."""
        n = len(p_values)
        if n == 0:
            return []
        indexed = sorted(enumerate(p_values), key=lambda x: x[1])
        corrected = [0.0] * n
        prev = 1.0
        for rank_minus_1 in range(n - 1, -1, -1):
            orig_idx, p = indexed[rank_minus_1]
            rank = rank_minus_1 + 1
            adjusted = min(p * n / rank, prev)
            adjusted = min(adjusted, 1.0)
            corrected[orig_idx] = adjusted
            prev = adjusted
        return corrected

    # ------------------------------------------------------------------
    # Single alignment analysis
    # ------------------------------------------------------------------

    def _run_single(self, alignment_path: str, tree, outgroup: str) -> List[Dict]:
        """Run Tajima's test for all ingroup pairs in one alignment."""
        alignment, _, _ = get_alignment_and_format_helper(alignment_path)

        # Build sequence dict
        seq_dict = {}
        for record in alignment:
            seq_dict[record.id] = str(record.seq)

        if outgroup not in seq_dict:
            raise PhykitUserError(
                [
                    f"Outgroup taxon '{outgroup}' from tree not found in alignment. "
                    f"Alignment taxa: {sorted(seq_dict.keys())}"
                ]
            )

        seq_out = seq_dict[outgroup]
        ingroup = sorted(t for t in seq_dict if t != outgroup)

        if len(ingroup) < 2:
            raise PhykitUserError(
                ["At least 2 ingroup taxa are required (3 total including outgroup)."]
            )

        results = []
        for t1, t2 in itertools.combinations(ingroup, 2):
            test_result = self._tajima_test(seq_dict[t1], seq_dict[t2], seq_out)
            test_result["taxon1"] = t1
            test_result["taxon2"] = t2
            results.append(test_result)

        # Multiple testing correction
        p_values = [r["p_value"] for r in results]
        p_bonf = self._bonferroni(p_values)
        p_fdr = self._fdr(p_values)
        for i, r in enumerate(results):
            r["p_bonf"] = p_bonf[i]
            r["p_fdr"] = p_fdr[i]

        return results

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self):
        tree = self.read_tree_file()

        outgroup = self._identify_outgroup(tree)

        if self.alignment_list:
            # Batch mode
            list_path = Path(self.alignment_list)
            if not list_path.exists():
                raise PhykitUserError([f"{self.alignment_list} not found."])
            aln_paths = [
                line.strip()
                for line in list_path.read_text().splitlines()
                if line.strip() and not line.strip().startswith("#")
            ]
            if not aln_paths:
                raise PhykitUserError(["No alignment paths found in list file."])

            # Run each alignment
            all_results = {}  # (t1, t2) -> list of per-gene results
            n_alignments = 0
            for aln_path in aln_paths:
                resolved = aln_path
                if not Path(aln_path).is_absolute():
                    resolved = str(list_path.parent / aln_path)
                try:
                    results = self._run_single(resolved, tree, outgroup)
                    n_alignments += 1
                    for r in results:
                        key = (r["taxon1"], r["taxon2"])
                        all_results.setdefault(key, []).append(r)
                except (SystemExit, ValueError, FileNotFoundError, KeyError):
                    continue  # skip problematic alignments

            self._output_batch(outgroup, n_alignments, all_results)

        elif self.alignment_file_path:
            # Single mode
            results = self._run_single(self.alignment_file_path, tree, outgroup)
            self._output_single(outgroup, results)
        else:
            raise PhykitUserError(
                ["Provide -a (single alignment) or -l (alignment list)."]
            )

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def _output_single(self, outgroup: str, results: List[Dict]):
        n_ingroup = len(
            set(r["taxon1"] for r in results) | set(r["taxon2"] for r in results)
        )

        if self.json_output:
            print_json(dict(
                outgroup=outgroup,
                n_ingroup_taxa=n_ingroup,
                n_tests=len(results),
                results=[
                    dict(
                        taxon1=r["taxon1"], taxon2=r["taxon2"],
                        m1=r["m1"], m2=r["m2"],
                        chi2=round(r["chi2"], 4),
                        p_value=r["p_value"],
                        p_bonferroni=r["p_bonf"],
                        p_fdr=r["p_fdr"],
                    )
                    for r in results
                ],
            ))
        else:
            try:
                print(f"Outgroup: {outgroup}")
                print(f"Number of ingroup taxa: {n_ingroup}")
                print(f"Number of pairwise tests: {len(results)}")
                print("---")
                print(
                    "taxon1\ttaxon2\tm1\tm2\tchi2\tp_value\tp_bonf\tp_fdr\tsignificant"
                )
                for r in results:
                    sig = "*" if r["p_fdr"] < 0.05 else ""
                    print(
                        f"{r['taxon1']}\t{r['taxon2']}\t{r['m1']}\t{r['m2']}\t"
                        f"{r['chi2']:.4f}\t{r['p_value']:.4e}\t{r['p_bonf']:.4e}\t"
                        f"{r['p_fdr']:.4e}\t{sig}"
                    )
            except BrokenPipeError:
                pass

    def _output_batch(self, outgroup: str, n_alignments: int, all_results: Dict):
        if self.json_output:
            pairs = []
            for (t1, t2), gene_results in sorted(all_results.items()):
                n_reject = sum(1 for r in gene_results if r["p_value"] < 0.05)
                chi2_values = [r["chi2"] for r in gene_results]
                chi2_values.sort()
                median_chi2 = (
                    chi2_values[len(chi2_values) // 2] if chi2_values else 0.0
                )
                pairs.append(dict(
                    taxon1=t1, taxon2=t2,
                    n_reject=n_reject, n_total=len(gene_results),
                    pct_reject=round(100 * n_reject / len(gene_results), 1),
                    median_chi2=round(median_chi2, 4),
                ))
            print_json(dict(
                outgroup=outgroup,
                n_alignments=n_alignments,
                pairs=pairs,
            ))
        else:
            try:
                print(f"Number of alignments: {n_alignments}")
                print(f"Outgroup: {outgroup}")
                print("---")
                print(
                    "taxon1\ttaxon2\tn_reject\tn_total\tpct_reject\tmedian_chi2"
                )
                for (t1, t2), gene_results in sorted(all_results.items()):
                    n_reject = sum(
                        1 for r in gene_results if r["p_value"] < 0.05
                    )
                    chi2_values = sorted(r["chi2"] for r in gene_results)
                    median_chi2 = (
                        chi2_values[len(chi2_values) // 2]
                        if chi2_values
                        else 0.0
                    )
                    pct = 100 * n_reject / len(gene_results)
                    print(
                        f"{t1}\t{t2}\t{n_reject}\t{len(gene_results)}\t"
                        f"{pct:.1f}\t{median_chi2:.4f}"
                    )
            except BrokenPipeError:
                pass
