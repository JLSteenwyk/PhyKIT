from __future__ import annotations

from collections import Counter

from .base import Alignment


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)

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
_DNA_GAP_CODES = None
_PROTEIN_GAP_CODES = None
_DNA_GAP_BYTES = b"-?*XN"
_PROTEIN_GAP_BYTES = b"-?*X"
_DNA_STANDARD_CODES = (65, 67, 71, 84)  # A, C, G, T


def _get_gap_codes(is_protein: bool):
    global _DNA_GAP_CODES, _PROTEIN_GAP_CODES
    if is_protein:
        if _PROTEIN_GAP_CODES is None:
            _PROTEIN_GAP_CODES = np.frombuffer(b"-?*X", dtype=np.uint8)
        return _PROTEIN_GAP_CODES

    if _DNA_GAP_CODES is None:
        _DNA_GAP_CODES = np.frombuffer(b"-?*XN", dtype=np.uint8)
    return _DNA_GAP_CODES


def _count_ascii_parsimony_informative_sites(
    alignment_array,
    valid_mask=None,
    block_size: int = 8192,
) -> int:
    aln_len = alignment_array.shape[1]
    pi_sites = 0
    for start in range(0, aln_len, block_size):
        stop = min(aln_len, start + block_size)
        block = alignment_array[:, start:stop]
        column_offsets = np.arange(stop - start, dtype=np.uint32) * 256
        encoded = (block + column_offsets).ravel()
        if valid_mask is not None:
            valid_block = valid_mask[:, start:stop]
            encoded = encoded[valid_block.ravel()]
        counts = np.bincount(
            encoded,
            minlength=(stop - start) * 256,
        ).reshape(stop - start, 256)
        recurrent_symbol_counts = np.count_nonzero(counts >= 2, axis=1)
        pi_sites += int(np.count_nonzero(recurrent_symbol_counts >= 2))
    return pi_sites


def _count_clean_dna_parsimony_informative_sites(alignment_array) -> int | None:
    recurrent = np.zeros(alignment_array.shape[1], dtype=np.uint8)
    standard_total = np.zeros(alignment_array.shape[1], dtype=np.intp)
    for code in _DNA_STANDARD_CODES:
        counts = np.count_nonzero(alignment_array == code, axis=0)
        standard_total += counts
        recurrent += counts >= 2

    if int(standard_total.sum()) != alignment_array.size:
        return None
    return int(np.count_nonzero(recurrent >= 2))


class ParsimonyInformative(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.json_output = parsed["json_output"]

    def run(self):
        alignment, _, is_protein = self.get_alignment_and_format()
        pi_sites, aln_len, pi_sites_per = self.calculate_parsimony_informative_sites(
            alignment, is_protein
        )

        if self.json_output:
            print_json(
                dict(
                    parsimony_informative_sites=pi_sites,
                    alignment_length=aln_len,
                    percent_parsimony_informative_sites=round(pi_sites_per, 4),
                )
            )
            return

        print(f"{pi_sites}\t{aln_len}\t{round(pi_sites_per, 4)}")

    def process_args(self, args) -> dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
        )

    def get_number_of_occurrences_per_character(
        self,
        alignment: MultipleSeqAlignment,
        idx: int,
        is_protein: bool = False,
    ) -> Counter:
        gap_chars = {char.upper() for char in self.get_gap_chars(is_protein)}
        counts = {}
        for record in alignment:
            char = record.seq[idx].upper()
            if char not in gap_chars:
                counts[char] = counts.get(char, 0) + 1

        return Counter(counts)

    def is_parsimony_informative(
        self,
        num_occurrences: Counter,
    ) -> bool:
        """
        Check if a site is parsimony informative.
        That is, the site has two characters that appear at least twice.
        """
        informative_char_count = 0
        for count in num_occurrences.values():
            if count >= 2:
                informative_char_count += 1
                if informative_char_count >= 2:
                    return True
        return False

    def calculate_parsimony_informative_sites(
        self,
        alignment: MultipleSeqAlignment,
        is_protein: bool = False,
    ) -> tuple[int, int, float]:
        aln_len = alignment.get_alignment_length()
        num_records = len(alignment)
        if num_records <= 1:
            return 0, aln_len, 0.0

        sequences = []
        first_sequence = None
        all_identical = True
        for record in alignment:
            sequence = str(record.seq).upper()
            if first_sequence is None:
                first_sequence = sequence
            elif all_identical and sequence != first_sequence:
                all_identical = False
            sequences.append(sequence)

        if aln_len == 0:
            return 0, aln_len, 0.0
        if all_identical:
            return 0, aln_len, 0.0

        try:
            alignment_bytes = "".join(sequences).encode("ascii")
            alignment_array = np.frombuffer(
                alignment_bytes,
                dtype=np.uint8,
            ).reshape(len(sequences), aln_len)
            gap_bytes = _PROTEIN_GAP_BYTES if is_protein else _DNA_GAP_BYTES
            if not any(code in alignment_bytes for code in gap_bytes):
                if not is_protein:
                    pi_sites = _count_clean_dna_parsimony_informative_sites(
                        alignment_array,
                    )
                    if pi_sites is not None:
                        pi_sites_per = (pi_sites / aln_len) * 100
                        return pi_sites, aln_len, pi_sites_per
                valid_mask = None
            else:
                invalid_mask = np.zeros(alignment_array.shape, dtype=np.bool_)
                for gap_code in _get_gap_codes(is_protein):
                    invalid_mask |= alignment_array == gap_code
                valid_mask = ~invalid_mask if invalid_mask.any() else None
            pi_sites = _count_ascii_parsimony_informative_sites(
                alignment_array,
                valid_mask,
            )
            pi_sites_per = (pi_sites / aln_len) * 100
            return pi_sites, aln_len, pi_sites_per
        except UnicodeEncodeError:
            gap_chars = {char.upper() for char in self.get_gap_chars(is_protein)}
            alignment_array = np.array([list(seq) for seq in sequences], dtype="U1")
            gap_chars_array = np.array(list(gap_chars), dtype="U1")
            valid_mask = ~np.isin(alignment_array, gap_chars_array)
        if alignment_array.size == 0:
            return 0, aln_len, 0.0

        valid_chars = np.unique(alignment_array[valid_mask])
        chars_appearing_twice = np.zeros(aln_len, dtype=np.uint16)
        for char in valid_chars:
            chars_appearing_twice += np.sum(alignment_array == char, axis=0) >= 2

        pi_sites = int(np.count_nonzero(chars_appearing_twice >= 2))
        pi_sites_per = (pi_sites / aln_len) * 100

        return pi_sites, aln_len, pi_sites_per
