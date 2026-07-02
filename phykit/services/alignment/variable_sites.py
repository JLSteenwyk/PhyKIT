from __future__ import annotations

from .base import Alignment, _all_sequences_identical


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)

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
_DNA_GAP_CODES = None
_PROTEIN_GAP_CODES = None
_DNA_GAP_BYTES = b"-?*XN"
_PROTEIN_GAP_BYTES = b"-?*X"


def _get_gap_codes(is_protein: bool):
    global _DNA_GAP_CODES, _PROTEIN_GAP_CODES
    if is_protein:
        if _PROTEIN_GAP_CODES is None:
            _PROTEIN_GAP_CODES = np.frombuffer(b"-?*X", dtype=np.uint8)
        return _PROTEIN_GAP_CODES

    if _DNA_GAP_CODES is None:
        _DNA_GAP_CODES = np.frombuffer(b"-?*XN", dtype=np.uint8)
    return _DNA_GAP_CODES


class VariableSites(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(alignment_file_path=parsed["alignment_file_path"])
        self.json_output = parsed["json_output"]

    def run(self):
        alignment, _, is_protein = self.get_alignment_and_format()
        var_sites, aln_len, var_sites_per = \
            self.calculate_variable_sites(alignment, is_protein)

        if self.json_output:
            print_json(
                dict(
                    variable_sites=var_sites,
                    alignment_length=aln_len,
                    percent_variable_sites=round(var_sites_per, 4),
                )
            )
            return

        print(f"{var_sites}\t{aln_len}\t{round(var_sites_per, 4)}")

    def process_args(self, args) -> dict[str, str]:
        return dict(
            alignment_file_path=args.alignment,
            json_output=getattr(args, "json", False),
        )

    def calculate_variable_sites(
        self,
        alignment: MultipleSeqAlignment,
        is_protein: bool = False,
    ) -> tuple[int, int, float]:
        aln_len = alignment.get_alignment_length()
        num_records = len(alignment)
        if num_records <= 1:
            return 0, aln_len, 0.0

        if aln_len == 0:
            return 0, aln_len, 0.0

        raw_sequences = [str(record.seq) for record in alignment]
        first_raw_sequence = raw_sequences[0]
        first_sequence = first_raw_sequence.upper()
        all_identical = True
        for idx in range(1, num_records):
            sequence = raw_sequences[idx]
            if sequence != first_raw_sequence and sequence.upper() != first_sequence:
                all_identical = False
                break

        if all_identical:
            return 0, aln_len, 0.0
        sequences = [sequence.upper() for sequence in raw_sequences]

        try:
            alignment_bytes = "".join(sequences).encode("ascii")
            alignment_array = np.frombuffer(
                alignment_bytes,
                dtype=np.uint8,
            ).reshape(len(sequences), aln_len)
            gap_bytes = _PROTEIN_GAP_BYTES if is_protein else _DNA_GAP_BYTES
            if not any(code in alignment_bytes for code in gap_bytes):
                variable_columns = (
                    alignment_array.min(axis=0) != alignment_array.max(axis=0)
                )
            else:
                invalid_mask = np.zeros(alignment_array.shape, dtype=np.bool_)
                for gap_code in _get_gap_codes(is_protein):
                    invalid_mask |= alignment_array == gap_code
                valid_min = np.where(invalid_mask, 255, alignment_array).min(axis=0)
                valid_max = np.where(invalid_mask, 0, alignment_array).max(axis=0)
                variable_columns = (valid_min != 255) & (valid_min != valid_max)
            var_sites = int(np.count_nonzero(variable_columns))
            var_sites_per = (var_sites / aln_len) * 100
            return var_sites, aln_len, var_sites_per
        except UnicodeEncodeError:
            gap_chars = {char.upper() for char in self.get_gap_chars(is_protein)}
            alignment_array = np.array([list(seq) for seq in sequences], dtype="U1")
            gap_chars_array = np.array(list(gap_chars), dtype="U1")
            valid_mask = ~np.isin(alignment_array, gap_chars_array)
        if alignment_array.size == 0:
            return 0, aln_len, 0.0

        valid_chars = np.unique(alignment_array[valid_mask])
        valid_symbol_counts = np.zeros(aln_len, dtype=np.uint16)
        for char in valid_chars:
            valid_symbol_counts += np.any(alignment_array == char, axis=0)

        var_sites = int(np.count_nonzero(valid_symbol_counts > 1))
        var_sites_per = (var_sites / aln_len) * 100

        return var_sites, aln_len, var_sites_per
