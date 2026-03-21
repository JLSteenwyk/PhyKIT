"""Phylogenetic GWAS following the Pease et al. (2016) approach.

Performs per-site association tests between alignment columns and a
phenotype, applies Benjamini-Hochberg FDR correction, optionally
classifies significant associations as monophyletic or polyphyletic
using a phylogenetic tree, and produces a Manhattan plot.
"""

import csv
import sys
from collections import Counter
from typing import Dict, List, Optional, Tuple

import numpy as np
from Bio import Phylo, SeqIO

from .base import Alignment
from ...helpers.json_output import print_json
from ...helpers.plot_config import PlotConfig


class PhyloGwas(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__()
        self.alignment_path = parsed["alignment_path"]
        self.phenotype_path = parsed["phenotype_path"]
        self.output_path = parsed["output_path"]
        self.tree_path = parsed["tree_path"]
        self.partition_path = parsed["partition_path"]
        self.alpha = parsed["alpha"]
        self.exclude_monophyletic = parsed["exclude_monophyletic"]
        self.csv_output = parsed["csv_output"]
        self.dot_size = parsed["dot_size"]
        self.json_output = parsed["json_output"]
        self.plot_config = parsed["plot_config"]

    def process_args(self, args) -> Dict[str, object]:
        return dict(
            alignment_path=args.alignment,
            phenotype_path=args.phenotype,
            output_path=args.output,
            tree_path=getattr(args, "tree", None),
            partition_path=getattr(args, "partition", None),
            alpha=getattr(args, "alpha", 0.05),
            exclude_monophyletic=getattr(args, "exclude_monophyletic", False),
            csv_output=getattr(args, "csv", None),
            dot_size=getattr(args, "dot_size", 1.0),
            json_output=getattr(args, "json", False),
            plot_config=PlotConfig.from_args(args),
        )

    # ------------------------------------------------------------------
    # Input parsing
    # ------------------------------------------------------------------

    @staticmethod
    def _read_fasta(path: str) -> Dict[str, str]:
        """Read FASTA alignment into {taxon: sequence} dict."""
        seqs = {}
        for record in SeqIO.parse(path, "fasta"):
            seqs[record.id] = str(record.seq).upper()
        return seqs

    @staticmethod
    def _read_phenotype(path: str) -> Dict[str, str]:
        """Read two-column TSV phenotype file into {taxon: value} dict."""
        pheno = {}
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                pheno[parts[0]] = parts[1]
        return pheno

    @staticmethod
    def _detect_phenotype_type(values: List[str]) -> str:
        """Return 'continuous' if all values castable to float, else 'categorical'."""
        for v in values:
            try:
                float(v)
            except ValueError:
                return "categorical"
        return "continuous"

    @staticmethod
    def _parse_partition_file(path: str) -> List[Tuple[str, int, int]]:
        """Parse RAxML-style partition file.

        Expected format per line (0- or 1-based):
            DNA, gene_name = start-end
        Returns list of (name, start_0based, end_0based_inclusive).
        """
        import re

        partitions = []
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                # Match: TYPE, name = start-end
                m = re.match(r"[^,]*,\s*(\S+)\s*=\s*(\d+)\s*-\s*(\d+)", line)
                if m:
                    name = m.group(1)
                    start = int(m.group(2)) - 1  # convert to 0-based
                    end = int(m.group(3)) - 1
                    partitions.append((name, start, end))
        return partitions

    # ------------------------------------------------------------------
    # Per-site association tests
    # ------------------------------------------------------------------

    @staticmethod
    def _is_ambiguous(char: str) -> bool:
        return char in {"-", "N", "X", "?", "*"}

    def _test_site_categorical(
        self, alleles: List[str], groups: List[str], unique_groups: List[str]
    ) -> Optional[Tuple[float, str, str, Dict[str, float]]]:
        """Fisher's exact or chi2 test for a categorical phenotype at one site.

        Returns (p_value, allele_0, allele_1, group_freqs) or None if site
        should be skipped.
        """
        from scipy.stats import chi2_contingency, fisher_exact

        # Determine unique alleles
        allele_counts = Counter(alleles)
        unique_alleles = sorted(allele_counts.keys())

        if len(unique_alleles) < 2:
            return None  # invariant
        if len(unique_alleles) > 2:
            return None  # skip multiallelic

        allele_0, allele_1 = unique_alleles
        # allele_0 = major, allele_1 = minor by frequency
        if allele_counts[allele_0] < allele_counts[allele_1]:
            allele_0, allele_1 = allele_1, allele_0

        n_groups = len(unique_groups)
        if n_groups == 2:
            # Build 2x2 table
            table = np.zeros((2, 2), dtype=int)
            for a, g in zip(alleles, groups):
                gi = 0 if g == unique_groups[0] else 1
                ai = 0 if a == allele_0 else 1
                table[gi, ai] += 1
            _, p_value = fisher_exact(table)
        else:
            # Build k x 2 table
            table = np.zeros((n_groups, 2), dtype=int)
            group_idx = {g: i for i, g in enumerate(unique_groups)}
            for a, g in zip(alleles, groups):
                gi = group_idx[g]
                ai = 0 if a == allele_0 else 1
                table[gi, ai] += 1
            # Check for zero-sum rows/columns
            if np.any(table.sum(axis=1) == 0) or np.any(table.sum(axis=0) == 0):
                return None
            _, p_value, _, _ = chi2_contingency(table)

        # Compute frequency of allele_1 per group
        group_freqs = {}
        for g in unique_groups:
            group_alleles = [a for a, grp in zip(alleles, groups) if grp == g]
            if group_alleles:
                group_freqs[g] = sum(1 for a in group_alleles if a == allele_1) / len(
                    group_alleles
                )
            else:
                group_freqs[g] = 0.0

        return p_value, allele_0, allele_1, group_freqs

    def _test_site_continuous(
        self, alleles: List[str], phenotype_values: List[float]
    ) -> Optional[Tuple[float, str, str, float]]:
        """Point-biserial correlation for a continuous phenotype at one site.

        Returns (p_value, allele_0, allele_1, correlation_r) or None.
        """
        from scipy.stats import pointbiserialr

        allele_counts = Counter(alleles)
        unique_alleles = sorted(allele_counts.keys())

        if len(unique_alleles) < 2:
            return None  # invariant
        if len(unique_alleles) > 2:
            return None  # skip multiallelic

        allele_0, allele_1 = unique_alleles
        if allele_counts[allele_0] < allele_counts[allele_1]:
            allele_0, allele_1 = allele_1, allele_0

        # Encode: 0 = major allele, 1 = minor allele
        binary = np.array([0 if a == allele_0 else 1 for a in alleles], dtype=int)
        values = np.array(phenotype_values, dtype=float)

        # Need variance in both arrays
        if np.std(binary) == 0 or np.std(values) == 0:
            return None

        r, p_value = pointbiserialr(binary, values)
        return p_value, allele_0, allele_1, float(r)

    # ------------------------------------------------------------------
    # Multiple testing correction
    # ------------------------------------------------------------------

    @staticmethod
    def _benjamini_hochberg(p_values: List[float]) -> np.ndarray:
        """Return BH-adjusted p-values."""
        n = len(p_values)
        if n == 0:
            return np.array([])
        p_arr = np.array(p_values, dtype=float)
        sorted_indices = np.argsort(p_arr)
        sorted_p = p_arr[sorted_indices]

        adjusted = np.zeros(n, dtype=float)
        for i in range(n - 1, -1, -1):
            if i == n - 1:
                adjusted[i] = sorted_p[i]
            else:
                adjusted[i] = min(sorted_p[i] * n / (i + 1), adjusted[i + 1])
        adjusted = np.minimum(adjusted, 1.0)

        result = np.zeros(n, dtype=float)
        result[sorted_indices] = adjusted
        return result

    # ------------------------------------------------------------------
    # Phylogenetic pattern classification
    # ------------------------------------------------------------------

    @staticmethod
    def _classify_phylo_pattern(tree, taxa_with_derived_allele: List[str]) -> str:
        """Check if taxa with the derived allele form a monophyletic group."""
        if len(taxa_with_derived_allele) <= 1:
            return "monophyletic"

        try:
            mrca = tree.common_ancestor(taxa_with_derived_allele)
        except Exception:
            return "polyphyletic"

        mrca_tips = {t.name for t in mrca.get_terminals()}
        derived_set = set(taxa_with_derived_allele)

        if derived_set == mrca_tips:
            return "monophyletic"
        else:
            return "polyphyletic"

    # ------------------------------------------------------------------
    # Partition gene lookup
    # ------------------------------------------------------------------

    @staticmethod
    def _position_to_gene(
        position: int, partitions: List[Tuple[str, int, int]]
    ) -> Optional[str]:
        """Return gene name for a given 0-based position, or None."""
        for name, start, end in partitions:
            if start <= position <= end:
                return name
        return None

    # ------------------------------------------------------------------
    # Main run
    # ------------------------------------------------------------------

    def run(self) -> None:
        # 1. Parse inputs
        seqs = self._read_fasta(self.alignment_path)
        phenotypes = self._read_phenotype(self.phenotype_path)

        shared_taxa = sorted(set(seqs.keys()) & set(phenotypes.keys()))
        if len(shared_taxa) < 4:
            print(
                f"Error: need at least 4 shared taxa between alignment and phenotype "
                f"file, found {len(shared_taxa)}."
            )
            sys.exit(2)

        # Validate alignment lengths
        aln_len = len(next(iter(seqs.values())))
        for taxon, seq in seqs.items():
            if len(seq) != aln_len:
                print(f"Error: sequences have unequal lengths.")
                sys.exit(2)

        # Detect phenotype type
        pheno_values = [phenotypes[t] for t in shared_taxa]
        pheno_type = self._detect_phenotype_type(pheno_values)
        is_continuous = pheno_type == "continuous"

        if is_continuous:
            pheno_float = [float(phenotypes[t]) for t in shared_taxa]
        else:
            unique_groups = sorted(set(pheno_values))
            group_counts = Counter(pheno_values)

        # Parse tree if provided
        tree = None
        if self.tree_path:
            tree = Phylo.read(self.tree_path, "newick")
            # Prune tree to shared taxa
            tree_tips = {t.name for t in tree.get_terminals()}
            tips_to_remove = tree_tips - set(shared_taxa)
            for tip_name in tips_to_remove:
                matching = [t for t in tree.get_terminals() if t.name == tip_name]
                for t in matching:
                    tree.prune(t)

        # Parse partitions if provided
        partitions = []
        if self.partition_path:
            partitions = self._parse_partition_file(self.partition_path)

        # 2. Per-site association tests
        site_results = []
        raw_p_values = []

        for col_idx in range(aln_len):
            # Extract alleles for shared taxa
            alleles = []
            skip = False
            for taxon in shared_taxa:
                char = seqs[taxon][col_idx]
                if self._is_ambiguous(char):
                    skip = True
                    break
                alleles.append(char)

            if skip:
                continue

            if is_continuous:
                result = self._test_site_continuous(alleles, pheno_float)
                if result is None:
                    continue
                p_value, allele_0, allele_1, corr_r = result
                site_results.append(
                    dict(
                        position=col_idx + 1,  # 1-based
                        allele_0=allele_0,
                        allele_1=allele_1,
                        p_value=p_value,
                        correlation_r=corr_r,
                    )
                )
            else:
                groups = [phenotypes[t] for t in shared_taxa]
                result = self._test_site_categorical(alleles, groups, unique_groups)
                if result is None:
                    continue
                p_value, allele_0, allele_1, group_freqs = result
                site_results.append(
                    dict(
                        position=col_idx + 1,  # 1-based
                        allele_0=allele_0,
                        allele_1=allele_1,
                        p_value=p_value,
                        group_freqs=group_freqs,
                    )
                )

            raw_p_values.append(p_value)

        # 3. Multiple testing correction
        if raw_p_values:
            adjusted = self._benjamini_hochberg(raw_p_values)
            for i, r in enumerate(site_results):
                r["fdr_p_value"] = float(adjusted[i])
                r["fdr_significant"] = bool(adjusted[i] < self.alpha)
        else:
            for r in site_results:
                r["fdr_p_value"] = 1.0
                r["fdr_significant"] = False

        # 4. Gene annotation and phylo pattern
        for r in site_results:
            r["gene"] = self._position_to_gene(r["position"] - 1, partitions)

            if tree is not None and r["fdr_significant"]:
                taxa_minor = [
                    t
                    for t, a in zip(shared_taxa, [seqs[t][r["position"] - 1] for t in shared_taxa])
                    if a == r["allele_1"]
                ]
                r["phylo_pattern"] = self._classify_phylo_pattern(tree, taxa_minor)
            elif tree is not None:
                r["phylo_pattern"] = None
            else:
                r["phylo_pattern"] = None

        # 5. Exclude monophyletic if requested
        if self.exclude_monophyletic:
            site_results = [
                r
                for r in site_results
                if not (r["fdr_significant"] and r.get("phylo_pattern") == "monophyletic")
            ]

        # Compute summary stats
        significant = [r for r in site_results if r["fdr_significant"]]
        n_polyphyletic = sum(
            1 for r in significant if r.get("phylo_pattern") == "polyphyletic"
        )
        n_monophyletic = sum(
            1 for r in significant if r.get("phylo_pattern") == "monophyletic"
        )

        # 6. Manhattan plot
        self._create_manhattan_plot(
            site_results, partitions, tree is not None
        )

        # 7. CSV output
        if self.csv_output:
            self._write_csv(site_results, is_continuous)

        # 8. Text / JSON output
        if self.json_output:
            payload = dict(
                n_taxa=len(shared_taxa),
                alignment_length=aln_len,
                phenotype_type=pheno_type,
                sites_tested=len(site_results),
                significant_sites=len(significant),
                alpha=self.alpha,
                results=[
                    {
                        k: v
                        for k, v in r.items()
                    }
                    for r in site_results
                ],
            )
            if not is_continuous:
                payload["groups"] = dict(group_counts)
            if tree is not None:
                payload["polyphyletic"] = n_polyphyletic
                payload["monophyletic"] = n_monophyletic
            print_json(payload)
            return

        # Text output
        lines = []
        lines.append("Phylogenetic GWAS")
        lines.append(f"Alignment: {self.alignment_path}")
        lines.append(f"Taxa: {len(shared_taxa)}")
        lines.append(f"Alignment length: {aln_len}")
        lines.append(f"Phenotype type: {pheno_type}")
        if not is_continuous:
            group_str = ", ".join(
                f"{g} ({group_counts[g]})" for g in unique_groups
            )
            lines.append(f"Groups: {group_str}")
        lines.append(f"Biallelic sites tested: {len(site_results)}")
        lines.append(f"Significant sites (FDR < {self.alpha}): {len(significant)}")
        if tree is not None:
            lines.append(f"  Polyphyletic: {n_polyphyletic}")
            lines.append(f"  Monophyletic: {n_monophyletic}")

        if significant:
            lines.append("")
            # Sort by p-value
            top_sig = sorted(significant, key=lambda r: r["p_value"])[:10]
            if is_continuous:
                lines.append("Top significant sites:")
                lines.append(
                    f"  {'Position':<10}{'Gene':<12}{'Allele':<10}"
                    f"{'r':<10}{'p-value':<12}{'FDR_p':<12}"
                    f"{'Pattern':<14}"
                )
                for r in top_sig:
                    gene = r["gene"] or "."
                    pattern = r["phylo_pattern"] or "."
                    lines.append(
                        f"  {r['position']:<10}{gene:<12}"
                        f"{r['allele_0']}>{r['allele_1']:<8}"
                        f"{r.get('correlation_r', 0):<10.4f}"
                        f"{r['p_value']:<12.4g}"
                        f"{r['fdr_p_value']:<12.4g}"
                        f"{pattern:<14}"
                    )
            else:
                lines.append("Top significant sites:")
                freq_header = "/".join(unique_groups)
                lines.append(
                    f"  {'Position':<10}{'Gene':<12}{'Allele':<10}"
                    f"{freq_header:<20}{'p-value':<12}{'FDR_p':<12}"
                    f"{'Pattern':<14}"
                )
                for r in top_sig:
                    gene = r["gene"] or "."
                    pattern = r["phylo_pattern"] or "."
                    freqs = "/".join(
                        f"{r['group_freqs'].get(g, 0):.2f}" for g in unique_groups
                    )
                    lines.append(
                        f"  {r['position']:<10}{gene:<12}"
                        f"{r['allele_0']}>{r['allele_1']:<8}"
                        f"{freqs:<20}"
                        f"{r['p_value']:<12.4g}"
                        f"{r['fdr_p_value']:<12.4g}"
                        f"{pattern:<14}"
                    )

        for line in lines:
            print(line)

    # ------------------------------------------------------------------
    # Manhattan plot
    # ------------------------------------------------------------------

    def _create_manhattan_plot(
        self,
        results: List[dict],
        partitions: List[Tuple[str, int, int]],
        has_tree: bool,
    ) -> None:
        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.patches import Patch
        except ImportError:
            print(
                "matplotlib is required for phylo_gwas. Install matplotlib and retry."
            )
            sys.exit(2)

        if not results:
            # Create an empty plot
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.set_xlabel("Alignment position")
            ax.set_ylabel("-log10(p-value)")
            ax.text(
                0.5, 0.5, "No biallelic sites tested",
                transform=ax.transAxes, ha="center", va="center",
            )
            fig.savefig(self.output_path, dpi=self.plot_config.dpi, bbox_inches="tight")
            plt.close(fig)
            return

        config = self.plot_config
        config.resolve(n_rows=None, n_cols=None)

        positions = np.array([r["position"] for r in results])
        p_values = np.array([r["p_value"] for r in results])
        # Avoid log10(0)
        p_values = np.maximum(p_values, 1e-300)
        neg_log_p = -np.log10(p_values)

        default_colors = ["#377eb8", "#e41a1c", "#999999"]
        colors = config.merge_colors(default_colors)
        color_nonsig = colors[0]
        color_polyphyletic = colors[1]
        color_monophyletic = colors[2]

        point_colors = []
        for r in results:
            if r["fdr_significant"]:
                if r.get("phylo_pattern") == "monophyletic":
                    point_colors.append(color_monophyletic)
                elif r.get("phylo_pattern") == "polyphyletic":
                    point_colors.append(color_polyphyletic)
                else:
                    point_colors.append(color_polyphyletic)  # no tree: significant = red
            else:
                point_colors.append(color_nonsig)

        fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
        ax.scatter(
            positions, neg_log_p, c=point_colors,
            s=8 * self.dot_size, alpha=0.6, edgecolors="none",
        )

        # Significance threshold line
        sig_pvals = [r["p_value"] for r in results if r["fdr_significant"]]
        if sig_pvals:
            threshold = -np.log10(max(sig_pvals))
            ax.axhline(
                y=threshold,
                color="red",
                linestyle="--",
                lw=0.8,
                label=f"FDR {self.alpha}",
            )

        # Partition annotations
        if partitions:
            for i, (name, start, end) in enumerate(partitions):
                if i % 2 == 0:
                    ax.axvspan(start + 1, end + 1, alpha=0.05, color="gray")
                mid = (start + end) / 2 + 1
                if len(partitions) <= 50:
                    ax.text(
                        mid,
                        ax.get_ylim()[1] * 0.98,
                        name,
                        ha="center",
                        fontsize=4,
                        rotation=90,
                    )

        ax.set_xlabel("Alignment position")
        ax.set_ylabel(r"-log$_{10}$(p-value)")

        # Legend
        legend_elements = [
            Patch(facecolor=color_nonsig, label="Not significant"),
        ]
        if has_tree:
            legend_elements.append(
                Patch(facecolor=color_polyphyletic, label="Significant (polyphyletic)")
            )
            legend_elements.append(
                Patch(facecolor=color_monophyletic, label="Significant (monophyletic)")
            )
        else:
            legend_elements.append(
                Patch(facecolor=color_polyphyletic, label="Significant")
            )
        legend_loc = config.legend_position or "upper right"
        if legend_loc != "none":
            ax.legend(handles=legend_elements, loc=legend_loc, fontsize=7)

        if config.show_title:
            ax.set_title(
                config.title or "Phylogenetic GWAS", fontsize=config.title_fontsize
            )
        if config.axis_fontsize:
            ax.xaxis.label.set_fontsize(config.axis_fontsize)
            ax.yaxis.label.set_fontsize(config.axis_fontsize)

        fig.tight_layout()
        fig.savefig(self.output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    # ------------------------------------------------------------------
    # CSV output
    # ------------------------------------------------------------------

    def _write_csv(self, results: List[dict], is_continuous: bool) -> None:
        fieldnames = [
            "position",
            "gene",
            "allele_0",
            "allele_1",
        ]
        if is_continuous:
            fieldnames.append("correlation_r")
        fieldnames.extend(
            [
                "p_value",
                "fdr_p_value",
                "significant",
                "phylo_pattern",
            ]
        )

        with open(self.csv_output, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            for r in results:
                row = dict(r)
                row["gene"] = row.get("gene") or ""
                row["significant"] = "yes" if row.get("fdr_significant") else "no"
                row["phylo_pattern"] = row.get("phylo_pattern") or ""
                writer.writerow(row)
