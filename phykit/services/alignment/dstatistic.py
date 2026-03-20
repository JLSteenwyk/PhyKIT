"""Patterson's D-statistic (ABBA-BABA test) for detecting introgression."""

from typing import Dict

import numpy as np
from Bio import SeqIO

from .base import Alignment
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


class Dstatistic(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_path"])
        self.p1 = parsed["p1"]
        self.p2 = parsed["p2"]
        self.p3 = parsed["p3"]
        self.outgroup = parsed["outgroup"]
        self.block_size = parsed["block_size"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> Dict[str, object]:
        return dict(
            alignment_path=args.alignment,
            p1=args.p1,
            p2=args.p2,
            p3=args.p3,
            outgroup=args.outgroup,
            block_size=getattr(args, "block_size", 100),
            json_output=getattr(args, "json", False),
        )

    def run(self):
        # Read alignment sequences
        sequences = {}
        for record in SeqIO.parse(self.alignment_file_path, "fasta"):
            sequences[record.id] = str(record.seq).upper()

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
        skip_chars = {"-", "N", "?", "X", "n", "x"}

        # Count site patterns
        abba_count = 0
        baba_count = 0

        for site in range(aln_length):
            p1 = seq_p1[site]
            p2 = seq_p2[site]
            p3 = seq_p3[site]
            o = seq_outgroup[site]

            # Skip sites with gaps or ambiguous characters
            if any(c in skip_chars for c in [p1, p2, p3, o]):
                continue

            # Skip sites that are not biallelic
            alleles = {p1, p2, p3, o}
            if len(alleles) != 2:
                continue

            # ABBA: P1=ancestral, P2=derived, P3=derived, O=ancestral
            if p1 == o and p2 != o and p3 != o and p2 == p3:
                abba_count += 1
            # BABA: P1=derived, P2=ancestral, P3=derived, O=ancestral
            elif p2 == o and p1 != o and p3 != o and p1 == p3:
                baba_count += 1

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
            block_abba = np.zeros(n_blocks)
            block_baba = np.zeros(n_blocks)

            for site in range(aln_length):
                block_idx = site // self.block_size
                if block_idx >= n_blocks:
                    break

                p1 = seq_p1[site]
                p2 = seq_p2[site]
                p3 = seq_p3[site]
                o = seq_outgroup[site]

                if any(c in skip_chars for c in [p1, p2, p3, o]):
                    continue
                alleles = {p1, p2, p3, o}
                if len(alleles) != 2:
                    continue

                if p1 == o and p2 != o and p3 != o and p2 == p3:
                    block_abba[block_idx] += 1
                elif p2 == o and p1 != o and p3 != o and p1 == p3:
                    block_baba[block_idx] += 1

            total_abba = np.sum(block_abba)
            total_baba = np.sum(block_baba)

            jackknife_d = np.zeros(n_blocks)
            for i in range(n_blocks):
                loo_abba = total_abba - block_abba[i]
                loo_baba = total_baba - block_baba[i]
                denom = loo_abba + loo_baba
                if denom > 0:
                    jackknife_d[i] = (loo_abba - loo_baba) / denom
                else:
                    jackknife_d[i] = 0.0

            mean_d = np.mean(jackknife_d)
            se = float(np.sqrt((n_blocks - 1) / n_blocks * np.sum((jackknife_d - mean_d) ** 2)))

            if se > 0:
                z_score = d_stat / se
                from scipy.stats import norm
                p_value = float(2.0 * norm.sf(abs(z_score)))
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

        try:
            print("Patterson's D-statistic (ABBA-BABA Test)")
            print("=========================================")
            print(f"Topology: ((({self.p1}, {self.p2}), {self.p3}), {self.outgroup})")
            print(f"P1: {self.p1}")
            print(f"P2: {self.p2}")
            print(f"P3: {self.p3}")
            print(f"Outgroup: {self.outgroup}")
            print()
            print(f"Alignment length: {aln_length}")
            print(f"Informative sites: {informative_sites}")
            print(f"ABBA sites: {abba_count}")
            print(f"BABA sites: {baba_count}")
            print(f"D-statistic: {d_stat:.4f}")

            if se is not None:
                print(f"Block jackknife (block size: {self.block_size}):")
                print(f"  Standard error: {se:.4f}")
                if z_score == float('inf'):
                    print("  Z-score: inf")
                else:
                    print(f"  Z-score: {z_score:.2f}")
                print(f"  p-value: {p_value:.6f}")
                print()
                print(f"Interpretation: {self._interpret(d_stat, p_value)}")
            else:
                print()
                print("Not enough blocks for jackknife significance test.")
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
