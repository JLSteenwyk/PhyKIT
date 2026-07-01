from __future__ import annotations

from .base import Alignment


class _LazyAlignIO:
    def read(self, *args, **kwargs):
        from Bio import AlignIO as _AlignIO

        return _AlignIO.read(*args, **kwargs)


AlignIO = _LazyAlignIO()


class _LazyNumpy:
    def __getattr__(self, name):
        import numpy as _np

        return getattr(_np, name)


np = _LazyNumpy()


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


def _alignment_size(alignment):
    try:
        return len(alignment)
    except TypeError:
        return None


def _alignment_length(alignment):
    get_alignment_length = getattr(alignment, "get_alignment_length", None)
    if get_alignment_length is None:
        return None
    return get_alignment_length()


class ColumnScore(Alignment):
    def __init__(self, args) -> None:
        parsed = self.process_args(args)
        super().__init__(fasta=parsed["fasta"], reference=parsed["reference"])
        self.json_output = parsed["json_output"]

    def run(self) -> None:
        query_records = AlignIO.read(self.fasta, "fasta")
        if self.reference == self.fasta:
            reference_records = query_records
        else:
            reference_records = AlignIO.read(self.reference, "fasta")

        direct_result = self._calculate_matches_between_alignments_direct(
            reference_records, query_records
        )
        if direct_result is None:
            ref_columns, query_columns = self.get_columns_from_alignments(
                reference_records, query_records
            )
            direct_result = self.calculate_matches_between_ref_and_query_columns(
                ref_columns, query_columns
            )
        number_of_matches, number_of_total_columns = direct_result

        score = round(number_of_matches / number_of_total_columns, 4)

        if self.json_output:
            print_json(dict(column_score=score))
            return

        print(score)

    def process_args(self, args) -> dict[str, str]:
        return dict(
            fasta=args.fasta,
            reference=args.reference,
            json_output=getattr(args, "json", False),
        )

    def get_columns_from_alignments(
        self,
        reference_records: MultipleSeqAlignment,
        query_records: MultipleSeqAlignment,
    ) -> tuple[list[str], list[str]]:
        ref_sequences = [str(record.seq).upper() for record in reference_records]
        query_sequences = [str(record.seq).upper() for record in query_records]

        ref_columns = ["".join(column) for column in zip(*ref_sequences)]
        query_columns = ["".join(column) for column in zip(*query_sequences)]

        return ref_columns, query_columns

    def calculate_matches_between_ref_and_query_columns(
        self,
        ref_columns: list[str],
        query_columns: list[str],
    ) -> tuple[int, int]:
        set1 = set(ref_columns)
        set2 = set(query_columns)

        matches = set1.intersection(set2)

        return len(matches), len(query_columns)

    @staticmethod
    def _unique_column_count_ascii(sequences: list[str], seq_len: int) -> int | None:
        try:
            sequence_matrix = np.frombuffer(
                "".join(sequences).encode("ascii"),
                dtype=np.uint8,
            ).reshape(len(sequences), seq_len)
        except UnicodeEncodeError:
            return None

        column_dtype = np.dtype((np.void, len(sequences)))
        unique_columns = np.unique(
            np.ascontiguousarray(sequence_matrix.T).view(column_dtype).ravel()
        )
        return int(unique_columns.size)

    @staticmethod
    def _repeated_sequence_symbols_ascii(sequences: list[str]) -> frozenset[str] | None:
        first_sequence = sequences[0]
        for sequence in sequences:
            if sequence != first_sequence:
                return None
        try:
            first_sequence.encode("ascii")
        except UnicodeEncodeError:
            return None
        return frozenset(first_sequence)

    @staticmethod
    def _calculate_matches_between_alignments_direct(
        reference_records: MultipleSeqAlignment,
        query_records: MultipleSeqAlignment,
    ) -> tuple[int, int] | None:
        if reference_records is query_records:
            query_sequences = [str(record.seq).upper() for record in query_records]
            if not query_sequences:
                return None

            query_len = len(query_sequences[0])
            if any(len(seq) != query_len for seq in query_sequences):
                return None

            query_symbols = ColumnScore._repeated_sequence_symbols_ascii(
                query_sequences
            )
            if query_symbols is not None:
                return len(query_symbols), query_len

            unique_column_count = ColumnScore._unique_column_count_ascii(
                query_sequences,
                query_len,
            )
            if unique_column_count is None:
                return None
            return unique_column_count, query_len

        ref_size = _alignment_size(reference_records)
        query_size = _alignment_size(query_records)
        if (
            ref_size is not None
            and query_size is not None
            and ref_size > 0
            and query_size > 0
            and ref_size != query_size
        ):
            query_len = _alignment_length(query_records)
            if query_len is not None:
                return 0, query_len

        ref_sequences = [str(record.seq).upper() for record in reference_records]
        query_sequences = [str(record.seq).upper() for record in query_records]
        if not ref_sequences or not query_sequences:
            return None

        ref_len = len(ref_sequences[0])
        query_len = len(query_sequences[0])
        if (
            any(len(seq) != ref_len for seq in ref_sequences)
            or any(len(seq) != query_len for seq in query_sequences)
        ):
            return None

        if len(ref_sequences) != len(query_sequences):
            return 0, query_len

        ref_symbols = ColumnScore._repeated_sequence_symbols_ascii(ref_sequences)
        query_symbols = ColumnScore._repeated_sequence_symbols_ascii(query_sequences)
        if ref_symbols is not None and query_symbols is not None:
            return len(ref_symbols.intersection(query_symbols)), query_len

        if ref_sequences == query_sequences:
            unique_column_count = ColumnScore._unique_column_count_ascii(
                query_sequences,
                query_len,
            )
            if unique_column_count is None:
                return None
            return unique_column_count, query_len

        n_taxa = len(ref_sequences)
        try:
            ref_matrix = np.frombuffer(
                "".join(ref_sequences).encode("ascii"),
                dtype=np.uint8,
            ).reshape(n_taxa, ref_len)
            query_matrix = np.frombuffer(
                "".join(query_sequences).encode("ascii"),
                dtype=np.uint8,
            ).reshape(n_taxa, query_len)
        except UnicodeEncodeError:
            return None

        column_dtype = np.dtype((np.void, n_taxa))
        ref_unique = np.unique(
            np.ascontiguousarray(ref_matrix.T).view(column_dtype).ravel()
        )
        query_unique = np.unique(
            np.ascontiguousarray(query_matrix.T).view(column_dtype).ravel()
        )
        matches = np.intersect1d(ref_unique, query_unique, assume_unique=True)
        return int(matches.size), query_len
