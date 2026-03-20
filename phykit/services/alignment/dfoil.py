"""DFOIL test (Pease & Hahn 2015) for detecting and polarizing introgression
in a 5-taxon symmetric phylogeny.

Topology: ((P1, P2), (P3, P4), Outgroup)
"""

from typing import Dict

from Bio import SeqIO

from .base import Alignment
from ...helpers.json_output import print_json
from ...errors import PhykitUserError


# All 16 binary site patterns for 5 taxa (P1, P2, P3, P4, Outgroup).
# A = matches outgroup (ancestral), B = differs (derived).
PATTERNS = [
    'AAAAA', 'AAABA', 'AABAA', 'AABBA',
    'ABAAA', 'ABABA', 'ABBAA', 'ABBBA',
    'BAAAA', 'BAABA', 'BABAA', 'BABBA',
    'BBAAA', 'BBABA', 'BBBAA', 'BBBBA',
]

# Invariant / uninformative patterns (all ancestral or all derived).
_UNINFORMATIVE = {'AAAAA', 'BBBBA'}

# Sign-pattern interpretation table (DFO, DIL, DFI, DOL).
INTERPRETATIONS = {
    '+++0': 'Introgression: P1 -> P3 (or P3 -> P1)',
    '--0+': 'Introgression: P1 -> P4 (or P4 -> P1)',
    '++-0': 'Introgression: P2 -> P3 (or P3 -> P2)',
    '--0-': 'Introgression: P2 -> P4 (or P4 -> P2)',
    '+0++': 'Introgression: P3 -> P1 (or P1 -> P3)',
    '-0++': 'Introgression: P4 -> P1 (or P1 -> P4)',
    '0+--': 'Introgression: P3 -> P2 (or P2 -> P3)',
    '0---': 'Introgression: P4 -> P2 (or P2 -> P4)',
    '++00': 'Introgression: ancestor of (P1,P2) <-> P3',
    '--00': 'Introgression: ancestor of (P1,P2) <-> P4',
    '0000': 'No significant introgression detected',
}


class Dfoil(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_path"])
        self.p1 = parsed["p1"]
        self.p2 = parsed["p2"]
        self.p3 = parsed["p3"]
        self.p4 = parsed["p4"]
        self.outgroup = parsed["outgroup"]
        self.json_output = parsed["json_output"]

    def process_args(self, args) -> Dict[str, object]:
        return dict(
            alignment_path=args.alignment,
            p1=args.p1,
            p2=args.p2,
            p3=args.p3,
            p4=args.p4,
            outgroup=args.outgroup,
            json_output=getattr(args, "json", False),
        )

    def run(self):
        # Read alignment sequences
        sequences = {}
        for record in SeqIO.parse(self.alignment_file_path, "fasta"):
            sequences[record.id] = str(record.seq).upper()

        # Validate taxa are present
        required = {
            "p1": self.p1,
            "p2": self.p2,
            "p3": self.p3,
            "p4": self.p4,
            "outgroup": self.outgroup,
        }
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
        seq_p4 = sequences[self.p4]
        seq_o = sequences[self.outgroup]

        # Validate equal lengths
        lengths = {len(seq_p1), len(seq_p2), len(seq_p3), len(seq_p4), len(seq_o)}
        if len(lengths) != 1:
            raise PhykitUserError(
                ["Sequences have different lengths. All sequences must be aligned."],
                code=2,
            )

        aln_length = len(seq_p1)
        skip_chars = {"-", "N", "?", "X", "n", "x"}

        # Initialize pattern counts
        counts: Dict[str, int] = {p: 0 for p in PATTERNS}

        for site in range(aln_length):
            p1 = seq_p1[site]
            p2 = seq_p2[site]
            p3 = seq_p3[site]
            p4 = seq_p4[site]
            o = seq_o[site]

            # Skip sites with gaps or ambiguous characters
            if any(c in skip_chars for c in [p1, p2, p3, p4, o]):
                continue

            # Skip sites that are not biallelic
            alleles = {p1, p2, p3, p4, o}
            if len(alleles) != 2:
                continue

            # Encode pattern: A if matches outgroup, B if differs
            pattern = ''.join(
                'A' if c == o else 'B'
                for c in [p1, p2, p3, p4, o]
            )
            counts[pattern] += 1

        # Count informative sites (exclude AAAAA and BBBBA)
        informative_sites = sum(
            v for k, v in counts.items() if k not in _UNINFORMATIVE
        )

        # Compute the four D-statistics
        dfo_left = counts['AAABA'] + counts['ABABA'] + counts['BABAA'] + counts['BBBAA']
        dfo_right = counts['AABAA'] + counts['ABBAA'] + counts['BAABA'] + counts['BBABA']

        dil_left = counts['AAABA'] + counts['ABBAA'] + counts['BAABA'] + counts['BBBAA']
        dil_right = counts['AABAA'] + counts['ABABA'] + counts['BABAA'] + counts['BBABA']

        dfi_left = counts['ABAAA'] + counts['ABABA'] + counts['BABAA'] + counts['BABBA']
        dfi_right = counts['BAAAA'] + counts['ABBAA'] + counts['BAABA'] + counts['ABBBA']

        dol_left = counts['ABAAA'] + counts['ABBAA'] + counts['BAABA'] + counts['BABBA']
        dol_right = counts['BAAAA'] + counts['ABABA'] + counts['BABAA'] + counts['ABBBA']

        DFO = (dfo_left - dfo_right) / (dfo_left + dfo_right) if (dfo_left + dfo_right) > 0 else 0.0
        DIL = (dil_left - dil_right) / (dil_left + dil_right) if (dil_left + dil_right) > 0 else 0.0
        DFI = (dfi_left - dfi_right) / (dfi_left + dfi_right) if (dfi_left + dfi_right) > 0 else 0.0
        DOL = (dol_left - dol_right) / (dol_left + dol_right) if (dol_left + dol_right) > 0 else 0.0

        # Chi-squared significance tests (1 df)
        from scipy.stats import chi2

        def _chi2_test(left, right):
            total = left + right
            if total == 0:
                return 0.0, 1.0
            chi2_stat = (left - right) ** 2 / total
            p_value = float(chi2.sf(chi2_stat, df=1))
            return float(chi2_stat), p_value

        dfo_chi2, dfo_p = _chi2_test(dfo_left, dfo_right)
        dil_chi2, dil_p = _chi2_test(dil_left, dil_right)
        dfi_chi2, dfi_p = _chi2_test(dfi_left, dfi_right)
        dol_chi2, dol_p = _chi2_test(dol_left, dol_right)

        # Sign pattern
        def _get_sign(d_value, p_value, alpha=0.05):
            if p_value >= alpha:
                return '0'
            return '+' if d_value > 0 else '-'

        sign_pattern = (
            _get_sign(DFO, dfo_p)
            + _get_sign(DIL, dil_p)
            + _get_sign(DFI, dfi_p)
            + _get_sign(DOL, dol_p)
        )

        interpretation = INTERPRETATIONS.get(
            sign_pattern, 'Ambiguous or complex introgression pattern'
        )

        # Significance stars helper
        def _stars(p):
            if p < 0.001:
                return ' ***'
            elif p < 0.01:
                return ' **'
            elif p < 0.05:
                return ' *'
            return ''

        # Output
        if self.json_output:
            # Pattern counts excluding AAAAA (always 'AAAAA' key exists)
            pattern_counts = {k: v for k, v in counts.items() if k not in _UNINFORMATIVE}

            payload = {
                "p1": self.p1,
                "p2": self.p2,
                "p3": self.p3,
                "p4": self.p4,
                "outgroup": self.outgroup,
                "alignment_length": aln_length,
                "informative_sites": informative_sites,
                "pattern_counts": pattern_counts,
                "dfo": {
                    "value": round(DFO, 4),
                    "left": dfo_left,
                    "right": dfo_right,
                    "chi2": round(dfo_chi2, 4),
                    "p_value": round(dfo_p, 6),
                },
                "dil": {
                    "value": round(DIL, 4),
                    "left": dil_left,
                    "right": dil_right,
                    "chi2": round(dil_chi2, 4),
                    "p_value": round(dil_p, 6),
                },
                "dfi": {
                    "value": round(DFI, 4),
                    "left": dfi_left,
                    "right": dfi_right,
                    "chi2": round(dfi_chi2, 4),
                    "p_value": round(dfi_p, 6),
                },
                "dol": {
                    "value": round(DOL, 4),
                    "left": dol_left,
                    "right": dol_right,
                    "chi2": round(dol_chi2, 4),
                    "p_value": round(dol_p, 6),
                },
                "sign_pattern": sign_pattern,
                "interpretation": interpretation,
            }
            print_json(payload, sort_keys=False)
            return

        try:
            print("DFOIL Test (Pease & Hahn 2015)")
            print("================================")
            print(f"Topology: (({self.p1}, {self.p2}), ({self.p3}, {self.p4}), {self.outgroup})")
            print(f"P1: {self.p1}, P2: {self.p2}, P3: {self.p3}, P4: {self.p4}, Outgroup: {self.outgroup}")
            print()
            print(f"Alignment length: {aln_length}")
            print(f"Informative sites: {informative_sites}")
            print()
            print("Site pattern counts:")
            # Print informative patterns in a compact layout
            informative_patterns = [p for p in PATTERNS if p not in _UNINFORMATIVE]
            for i in range(0, len(informative_patterns), 4):
                chunk = informative_patterns[i:i + 4]
                parts = [f"{p}: {counts[p]}" for p in chunk]
                print("  " + "  ".join(parts))
            print()
            print("D-statistics:")
            print(f"  DFO:  {DFO:.4f}  (p = {dfo_p:.6f}{_stars(dfo_p)})")
            print(f"  DIL:  {DIL:.4f}  (p = {dil_p:.6f}{_stars(dil_p)})")
            print(f"  DFI:  {DFI:.4f}  (p = {dfi_p:.6f}{_stars(dfi_p)})")
            print(f"  DOL:  {DOL:.4f}  (p = {dol_p:.6f}{_stars(dol_p)})")
            print()
            print(f"Sign pattern: {sign_pattern}")
            print(f"Interpretation: {interpretation}")
        except BrokenPipeError:
            pass
