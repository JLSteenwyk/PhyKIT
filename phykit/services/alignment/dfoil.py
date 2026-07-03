"""DFOIL test (Pease & Hahn 2015) for detecting and polarizing introgression
in a 5-taxon symmetric phylogeny.

Topology: ((P1, P2), (P3, P4), Outgroup)
"""

from __future__ import annotations

import math

from ._fasta import read_fasta_first_token_upper
from .base import Alignment
from ...errors import PhykitUserError


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


# All 16 binary site patterns for 5 taxa (P1, P2, P3, P4, Outgroup).
# A = matches outgroup (ancestral), B = differs (derived).
PATTERNS = [
    'AAAAA', 'AAABA', 'AABAA', 'AABBA',
    'ABAAA', 'ABABA', 'ABBAA', 'ABBBA',
    'BAAAA', 'BAABA', 'BABAA', 'BABBA',
    'BBAAA', 'BBABA', 'BBBAA', 'BBBBA',
]
_ZERO_PATTERN_COUNTS = dict.fromkeys(PATTERNS, 0)

# Invariant / uninformative patterns (all ancestral or all derived).
_UNINFORMATIVE = {'AAAAA', 'BBBBA'}
_INFORMATIVE_PATTERNS = tuple(
    pattern for pattern in PATTERNS if pattern not in _UNINFORMATIVE
)
_INFORMATIVE_PATTERN_GROUPS = tuple(
    _INFORMATIVE_PATTERNS[i:i + 4]
    for i in range(0, len(_INFORMATIVE_PATTERNS), 4)
)
_SKIP_CODES = (ord("-"), ord("N"), ord("?"), ord("X"), ord("n"), ord("x"))
_SKIP_BYTES = b"-N?Xnx"
_SCALAR_SKIP_CHARS = "-N?Xnx"
_SKIP_SCAN_BYTES = 4096
_SKIP_LOOKUP_SMALL_ALIGNMENT_MAX = 8192
_SKIP_LOOKUP = None

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


def _chi2_sf_df1(chi2_stat: float) -> float:
    return math.erfc(math.sqrt(chi2_stat / 2.0))


def _stars(p):
    if p < 0.001:
        return ' ***'
    elif p < 0.01:
        return ' **'
    elif p < 0.05:
        return ' *'
    return ''


def _get_skip_lookup():
    global _SKIP_LOOKUP
    if _SKIP_LOOKUP is None:
        lookup = np.zeros(256, dtype=np.bool_)
        lookup[np.frombuffer(_SKIP_BYTES, dtype=np.uint8)] = True
        _SKIP_LOOKUP = lookup
    return _SKIP_LOOKUP


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

    def process_args(self, args) -> dict[str, object]:
        return dict(
            alignment_path=args.alignment,
            p1=args.p1,
            p2=args.p2,
            p3=args.p3,
            p4=args.p4,
            outgroup=args.outgroup,
            json_output=getattr(args, "json", False),
        )

    @staticmethod
    def _read_fasta(path: str) -> dict[str, str]:
        return read_fasta_first_token_upper(path)

    @staticmethod
    def _count_site_patterns(seq_p1, seq_p2, seq_p3, seq_p4, seq_outgroup):
        """Count DFOIL binary site patterns."""
        if (
            seq_p1 == seq_outgroup
            and seq_p2 == seq_outgroup
            and seq_p3 == seq_outgroup
            and seq_p4 == seq_outgroup
        ):
            return _ZERO_PATTERN_COUNTS.copy()

        try:
            p1_bytes = seq_p1.encode("ascii")
            p2_bytes = seq_p2.encode("ascii")
            p3_bytes = seq_p3.encode("ascii")
            p4_bytes = seq_p4.encode("ascii")
            out_bytes = seq_outgroup.encode("ascii")
            p1 = np.frombuffer(p1_bytes, dtype=np.uint8)
            p2 = np.frombuffer(p2_bytes, dtype=np.uint8)
            p3 = np.frombuffer(p3_bytes, dtype=np.uint8)
            p4 = np.frombuffer(p4_bytes, dtype=np.uint8)
            out = np.frombuffer(out_bytes, dtype=np.uint8)
        except UnicodeEncodeError:
            return Dfoil._count_site_patterns_scalar(
                seq_p1,
                seq_p2,
                seq_p3,
                seq_p4,
                seq_outgroup,
            )

        has_skip_code = any(
            code in p1_bytes[:_SKIP_SCAN_BYTES]
            or code in p2_bytes[:_SKIP_SCAN_BYTES]
            or code in p3_bytes[:_SKIP_SCAN_BYTES]
            or code in p4_bytes[:_SKIP_SCAN_BYTES]
            or code in out_bytes[:_SKIP_SCAN_BYTES]
            or code in p1_bytes[-_SKIP_SCAN_BYTES:]
            or code in p2_bytes[-_SKIP_SCAN_BYTES:]
            or code in p3_bytes[-_SKIP_SCAN_BYTES:]
            or code in p4_bytes[-_SKIP_SCAN_BYTES:]
            or code in out_bytes[-_SKIP_SCAN_BYTES:]
            for code in _SKIP_BYTES
        )
        if not has_skip_code:
            has_skip_code = any(
                code in p1_bytes
                or code in p2_bytes
                or code in p3_bytes
                or code in p4_bytes
                or code in out_bytes
                for code in _SKIP_BYTES
            )

        diff1 = p1 != out
        diff2 = p2 != out
        diff3 = p3 != out
        diff4 = p4 != out
        any_derived = diff1 | diff2 | diff3 | diff4
        same_derived = (
            ((~diff1) | (~diff2) | (p1 == p2))
            & ((~diff1) | (~diff3) | (p1 == p3))
            & ((~diff1) | (~diff4) | (p1 == p4))
            & ((~diff2) | (~diff3) | (p2 == p3))
            & ((~diff2) | (~diff4) | (p2 == p4))
            & ((~diff3) | (~diff4) | (p3 == p4))
        )

        if has_skip_code:
            if len(p1) < _SKIP_LOOKUP_SMALL_ALIGNMENT_MAX:
                skip_lookup = _get_skip_lookup()
                valid = ~(
                    skip_lookup[p1]
                    | skip_lookup[p2]
                    | skip_lookup[p3]
                    | skip_lookup[p4]
                    | skip_lookup[out]
                )
            else:
                valid = np.ones(len(p1), dtype=bool)
                for code in _SKIP_CODES:
                    valid &= (
                        (p1 != code)
                        & (p2 != code)
                        & (p3 != code)
                        & (p4 != code)
                        & (out != code)
                    )
            biallelic = valid & any_derived & same_derived
        else:
            biallelic = any_derived & same_derived

        pattern_codes = (
            (diff1.astype(np.uint8) << 3)
            | (diff2.astype(np.uint8) << 2)
            | (diff3.astype(np.uint8) << 1)
            | diff4.astype(np.uint8)
        )
        bincounts = np.bincount(pattern_codes[biallelic], minlength=len(PATTERNS))
        if not has_skip_code:
            return dict(zip(PATTERNS, map(int, bincounts)))
        return {pattern: int(bincounts[i]) for i, pattern in enumerate(PATTERNS)}

    @staticmethod
    def _count_site_patterns_scalar(seq_p1, seq_p2, seq_p3, seq_p4, seq_outgroup):
        counts: dict[str, int] = {p: 0 for p in PATTERNS}
        patterns = PATTERNS
        skip_chars = _SCALAR_SKIP_CHARS

        for p1, p2, p3, p4, o in zip(
            seq_p1, seq_p2, seq_p3, seq_p4, seq_outgroup
        ):
            if (
                p1 in skip_chars
                or p2 in skip_chars
                or p3 in skip_chars
                or p4 in skip_chars
                or o in skip_chars
            ):
                continue

            diff1 = p1 != o
            diff2 = p2 != o
            diff3 = p3 != o
            diff4 = p4 != o
            if not (diff1 or diff2 or diff3 or diff4):
                continue

            if diff1:
                derived = p1
            elif diff2:
                derived = p2
            elif diff3:
                derived = p3
            else:
                derived = p4

            if (
                (diff1 and p1 != derived)
                or (diff2 and p2 != derived)
                or (diff3 and p3 != derived)
                or (diff4 and p4 != derived)
            ):
                continue

            pattern_code = (
                (8 if diff1 else 0)
                | (4 if diff2 else 0)
                | (2 if diff3 else 0)
                | (1 if diff4 else 0)
            )
            counts[patterns[pattern_code]] += 1

        return counts

    def _print_text_output(
        self,
        aln_length,
        informative_sites,
        counts,
        DFO,
        DIL,
        DFI,
        DOL,
        dfo_p,
        dil_p,
        dfi_p,
        dol_p,
        sign_pattern,
        interpretation,
    ):
        pattern_lines = "\n".join(
            "  " + "  ".join(f"{pattern}: {counts[pattern]}" for pattern in group)
            for group in _INFORMATIVE_PATTERN_GROUPS
        )

        try:
            print(
                f"DFOIL Test (Pease & Hahn 2015)\n"
                f"================================\n"
                f"Topology: (({self.p1}, {self.p2}), "
                f"({self.p3}, {self.p4}), {self.outgroup})\n"
                f"P1: {self.p1}, P2: {self.p2}, P3: {self.p3}, "
                f"P4: {self.p4}, Outgroup: {self.outgroup}\n"
                f"\n"
                f"Alignment length: {aln_length}\n"
                f"Informative sites: {informative_sites}\n"
                f"\n"
                f"Site pattern counts:\n"
                f"{pattern_lines}\n"
                f"\n"
                f"D-statistics:\n"
                f"  DFO:  {DFO:.4f}  (p = {dfo_p:.6f}{_stars(dfo_p)})\n"
                f"  DIL:  {DIL:.4f}  (p = {dil_p:.6f}{_stars(dil_p)})\n"
                f"  DFI:  {DFI:.4f}  (p = {dfi_p:.6f}{_stars(dfi_p)})\n"
                f"  DOL:  {DOL:.4f}  (p = {dol_p:.6f}{_stars(dol_p)})\n"
                f"\n"
                f"Sign pattern: {sign_pattern}\n"
                f"Interpretation: {interpretation}"
            )
        except BrokenPipeError:
            pass

    def run(self):
        # Read alignment sequences
        sequences = self._read_fasta(self.alignment_file_path)

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
        counts = self._count_site_patterns(seq_p1, seq_p2, seq_p3, seq_p4, seq_o)

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

        def _chi2_test(left, right):
            total = left + right
            if total == 0:
                return 0.0, 1.0
            chi2_stat = (left - right) ** 2 / total
            p_value = _chi2_sf_df1(chi2_stat)
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

        self._print_text_output(
            aln_length,
            informative_sites,
            counts,
            DFO,
            DIL,
            DFI,
            DOL,
            dfo_p,
            dil_p,
            dfi_p,
            dol_p,
            sign_pattern,
            interpretation,
        )
