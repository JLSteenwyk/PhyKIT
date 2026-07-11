from __future__ import annotations

import itertools
import math
import statistics

from .base import Alignment
from ...errors import PhykitUserError


_DNA_SYMBOLS = frozenset("ACGTRYSWKMBDHVN-?.U")
_UNAMBIGUOUS_DNA = frozenset("ACGT")
_CALCULATE_DN_DS = None
_BIO_ALIGNMENT = None
_BIO_SEQ = None
_CODON_TABLES = None


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def _load_codon_api():
    global _CALCULATE_DN_DS, _BIO_ALIGNMENT, _BIO_SEQ, _CODON_TABLES

    if _CALCULATE_DN_DS is None:
        from Bio.Align import Alignment as BioAlignment
        from Bio.Align import analysis
        from Bio.Data import CodonTable
        from Bio.Seq import Seq

        _CALCULATE_DN_DS = analysis.calculate_dn_ds
        _BIO_ALIGNMENT = BioAlignment
        _BIO_SEQ = Seq
        _CODON_TABLES = CodonTable.unambiguous_dna_by_id

    return _CALCULATE_DN_DS, _BIO_ALIGNMENT, _BIO_SEQ, _CODON_TABLES


def _mean(values):
    return statistics.fmean(values) if values else None


def _median(values):
    return statistics.median(values) if values else None


def _format_number(value):
    if value is None:
        return "NA"
    return f"{value:.8g}"


class CodonDnDs(Alignment):
    """Estimate pairwise dN, dS, and omega from a codon alignment."""

    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.method = parsed["method"]
        self.genetic_code = parsed["genetic_code"]
        self.kappa = parsed["kappa"]
        self.codon_frequency = parsed["codon_frequency"]
        self.stop_policy = parsed["stop_policy"]
        self.reference = parsed["reference"]
        self.verbose = parsed["verbose"]
        self.json_output = parsed["json_output"]
        self._validate_options()

    @staticmethod
    def process_args(args) -> dict:
        method = getattr(args, "method", "NG86").upper()
        codon_frequency = getattr(args, "codon_frequency", "F3X4").upper()
        return dict(
            alignment_file_path=args.alignment,
            method=method,
            genetic_code=getattr(args, "genetic_code", 1),
            kappa=getattr(args, "kappa", 1.0),
            codon_frequency=codon_frequency.replace("X", "x"),
            stop_policy=getattr(args, "stop_policy", "error"),
            reference=getattr(args, "reference", None),
            verbose=getattr(args, "verbose", False),
            json_output=getattr(args, "json", False),
        )

    def _validate_options(self) -> None:
        if self.method not in {"NG86", "LWL85", "YN00", "ML"}:
            raise PhykitUserError(
                [f"Unknown dN/dS method: {self.method}. Use NG86, LWL85, YN00, or ML."]
            )
        if not math.isfinite(self.kappa) or self.kappa <= 0.0:
            raise PhykitUserError(["--kappa must be a finite value greater than 0."])
        if self.codon_frequency not in {"F1x4", "F3x4", "F61"}:
            raise PhykitUserError(
                ["--codon-frequency must be F1x4, F3x4, or F61."]
            )
        if self.stop_policy not in {"error", "skip"}:
            raise PhykitUserError(["--stop-policy must be error or skip."])

    def run(self) -> None:
        alignment, _, is_protein = self.get_alignment_and_format()
        codon_table = self._get_codon_table()
        records = self._validate_alignment(alignment, is_protein, codon_table)
        pairs = self._select_pairs(records)
        results = [
            self._estimate_pair(left, right, codon_table) for left, right in pairs
        ]
        summary = self._summarize(results)

        if self.json_output:
            payload = {
                "method": self.method,
                "genetic_code": self.genetic_code,
                "kappa": self.kappa,
                "codon_frequency": (
                    self.codon_frequency if self.method == "ML" else None
                ),
                "reference": self.reference,
                "summary": summary,
            }
            if self.verbose:
                payload["pairs"] = results
            print_json(payload)
            return

        if self.verbose:
            self._print_pairs(results)
        else:
            self._print_summary(summary)

    def _get_codon_table(self):
        _, _, _, tables = _load_codon_api()
        try:
            return tables[self.genetic_code]
        except KeyError:
            available = ", ".join(str(code) for code in sorted(tables))
            raise PhykitUserError(
                [
                    f"Unknown NCBI genetic code: {self.genetic_code}.",
                    f"Available genetic-code IDs: {available}",
                ]
            )

    def _validate_alignment(self, alignment, is_protein, codon_table):
        if is_protein:
            raise PhykitUserError(
                ["Codon dN/dS analysis requires a nucleotide alignment."]
            )
        if len(alignment) < 2:
            raise PhykitUserError(
                ["Codon dN/dS analysis requires at least two sequences."]
            )

        alignment_length = alignment.get_alignment_length()
        if alignment_length == 0:
            raise PhykitUserError(["Codon alignment is empty."])
        if alignment_length % 3:
            raise PhykitUserError(
                [
                    f"Codon alignment length ({alignment_length}) is not divisible by 3.",
                    "The alignment must begin at the first codon position and preserve frame.",
                ]
            )

        stop_codons = frozenset(codon_table.stop_codons)
        last_codon_start = alignment_length - 3
        seen_ids = set()
        records = []
        for record in alignment:
            taxon = record.id
            if not taxon:
                raise PhykitUserError(["Codon alignment contains an empty taxon ID."])
            if taxon in seen_ids:
                raise PhykitUserError(
                    [f"Codon alignment contains duplicate taxon ID: {taxon}"]
                )
            seen_ids.add(taxon)

            sequence = str(record.seq).upper().replace("U", "T").replace(".", "-")
            if len(sequence) != alignment_length:
                raise PhykitUserError(
                    [
                        f"Sequence {taxon!r} has length {len(sequence)}, but the alignment length is {alignment_length}.",
                        "All codon-alignment sequences must have equal length.",
                    ]
                )
            invalid = sorted(set(sequence) - _DNA_SYMBOLS)
            if invalid:
                symbols = ", ".join(repr(symbol) for symbol in invalid)
                raise PhykitUserError(
                    [f"Sequence {taxon!r} contains invalid DNA symbol(s): {symbols}"]
                )

            for start in range(0, alignment_length, 3):
                codon = sequence[start : start + 3]
                if "-" in codon and codon != "---":
                    raise PhykitUserError(
                        [
                            f"Sequence {taxon!r} has a frame-disrupting gap in codon {start // 3 + 1}: {codon}",
                            "Codon-alignment gaps must occupy complete codons (---).",
                        ]
                    )
                if (
                    start != last_codon_start
                    and set(codon) <= _UNAMBIGUOUS_DNA
                    and codon in stop_codons
                    and self.stop_policy == "error"
                ):
                    raise PhykitUserError(
                        [
                            f"Sequence {taxon!r} contains an internal stop codon at codon {start // 3 + 1}: {codon}",
                            "Use --stop-policy skip only when internal stops should be excluded pairwise.",
                        ]
                    )

            records.append((taxon, sequence))
        return records

    def _select_pairs(self, records):
        if self.reference is None:
            return itertools.combinations(records, 2)

        reference_record = None
        for record in records:
            if record[0] == self.reference:
                reference_record = record
                break
        if reference_record is None:
            raise PhykitUserError(
                [f"Reference taxon {self.reference!r} is not present in the alignment."]
            )
        return (
            (reference_record, record)
            for record in records
            if record is not reference_record
        )

    @staticmethod
    def _prepare_pair(sequence_a, sequence_b, stop_codons):
        clean_a = []
        clean_b = []
        skipped = 0
        for start in range(0, len(sequence_a), 3):
            codon_a = sequence_a[start : start + 3]
            codon_b = sequence_b[start : start + 3]
            if (
                set(codon_a) <= _UNAMBIGUOUS_DNA
                and set(codon_b) <= _UNAMBIGUOUS_DNA
                and codon_a not in stop_codons
                and codon_b not in stop_codons
            ):
                clean_a.append(codon_a)
                clean_b.append(codon_b)
            else:
                skipped += 1
        return "".join(clean_a), "".join(clean_b), skipped

    def _estimate_pair(self, left, right, codon_table):
        taxon_a, sequence_a = left
        taxon_b, sequence_b = right
        clean_a, clean_b, skipped = self._prepare_pair(
            sequence_a, sequence_b, frozenset(codon_table.stop_codons)
        )
        codons_used = len(clean_a) // 3
        result = {
            "taxon_a": taxon_a,
            "taxon_b": taxon_b,
            "dN": None,
            "dS": None,
            "omega": None,
            "codons_used": codons_used,
            "codons_skipped": skipped,
            "status": "no_valid_codons",
        }
        if not clean_a:
            return result

        calculate_dn_ds, bio_alignment, bio_seq, _ = _load_codon_api()
        pair_alignment = bio_alignment([bio_seq(clean_a), bio_seq(clean_b)])
        kwargs = {
            "method": self.method,
            "codon_table": codon_table,
        }
        if self.method == "NG86":
            kwargs["k"] = self.kappa
        elif self.method == "ML":
            kwargs["cfreq"] = self.codon_frequency

        try:
            d_n, d_s = calculate_dn_ds(pair_alignment, **kwargs)
            d_n = float(d_n)
            d_s = float(d_s)
        except (ArithmeticError, FloatingPointError, RuntimeError, ValueError) as err:
            result["status"] = "estimation_failed"
            result["message"] = str(err) or err.__class__.__name__
            return result

        valid_d_n = math.isfinite(d_n) and d_n >= 0.0
        valid_d_s = math.isfinite(d_s) and d_s >= 0.0
        if valid_d_n:
            result["dN"] = d_n
        if valid_d_s:
            result["dS"] = d_s

        if not valid_d_n or not valid_d_s:
            result["status"] = "saturated_or_undefined"
        elif d_s == 0.0:
            result["status"] = "zero_synonymous_distance"
        else:
            omega = d_n / d_s
            if math.isfinite(omega):
                result["omega"] = omega
                result["status"] = "ok"
            else:
                result["status"] = "saturated_or_undefined"
        return result

    @staticmethod
    def _summarize(results):
        d_n_values = [row["dN"] for row in results if row["dN"] is not None]
        d_s_values = [row["dS"] for row in results if row["dS"] is not None]
        omega_values = [
            row["omega"] for row in results if row["omega"] is not None
        ]
        jointly_estimated = [
            row
            for row in results
            if row["dN"] is not None
            and row["dS"] is not None
            and row["dS"] > 0.0
        ]
        joint_mean_d_n = _mean([row["dN"] for row in jointly_estimated])
        joint_mean_d_s = _mean([row["dS"] for row in jointly_estimated])
        ratio_of_means = None
        if joint_mean_d_s is not None and joint_mean_d_s > 0.0:
            ratio_of_means = joint_mean_d_n / joint_mean_d_s

        status_counts = {}
        for row in results:
            status = row["status"]
            status_counts[status] = status_counts.get(status, 0) + 1

        return {
            "pairs_total": len(results),
            "pairs_with_dN": len(d_n_values),
            "pairs_with_dS": len(d_s_values),
            "pairs_with_omega": len(omega_values),
            "mean_dN": _mean(d_n_values),
            "median_dN": _median(d_n_values),
            "mean_dS": _mean(d_s_values),
            "median_dS": _median(d_s_values),
            "mean_omega": _mean(omega_values),
            "median_omega": _median(omega_values),
            "ratio_of_mean_dN_to_mean_dS": ratio_of_means,
            "status_counts": status_counts,
        }

    @staticmethod
    def _print_pairs(results):
        try:
            print(
                "taxon_a\ttaxon_b\tdN\tdS\tomega\tcodons_used\t"
                "codons_skipped\tstatus"
            )
            for row in results:
                print(
                    f"{row['taxon_a']}\t{row['taxon_b']}\t"
                    f"{_format_number(row['dN'])}\t"
                    f"{_format_number(row['dS'])}\t"
                    f"{_format_number(row['omega'])}\t"
                    f"{row['codons_used']}\t{row['codons_skipped']}\t"
                    f"{row['status']}"
                )
        except BrokenPipeError:
            pass

    def _print_summary(self, summary):
        try:
            lines = [
                f"method: {self.method}",
                f"genetic code: {self.genetic_code}",
                f"sequence pairs: {summary['pairs_total']}",
                f"pairs with omega: {summary['pairs_with_omega']}",
                f"mean dN: {_format_number(summary['mean_dN'])}",
                f"mean dS: {_format_number(summary['mean_dS'])}",
                "dN/dS (ratio of means): "
                f"{_format_number(summary['ratio_of_mean_dN_to_mean_dS'])}",
                f"mean pairwise omega: {_format_number(summary['mean_omega'])}",
                f"median pairwise omega: {_format_number(summary['median_omega'])}",
                "pair statuses: "
                + ", ".join(
                    f"{status}={count}"
                    for status, count in sorted(summary["status_counts"].items())
                ),
            ]
            print("\n".join(lines))
        except BrokenPipeError:
            pass
