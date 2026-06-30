from __future__ import annotations

import sys
from itertools import compress
from operator import itemgetter

from .base import Alignment


def print_json(*args, **kwargs):
    from ...helpers.json_output import print_json as _print_json

    return _print_json(*args, **kwargs)


class _LazySeqIO:
    def parse(self, *args, **kwargs):
        from Bio import SeqIO as _SeqIO

        return _SeqIO.parse(*args, **kwargs)

    def to_dict(self, *args, **kwargs):
        from Bio import SeqIO as _SeqIO

        return _SeqIO.to_dict(*args, **kwargs)


SeqIO = _LazySeqIO()
_FULL_SEGMENT = object()
_KEEP_CODON_MASK = (True, True, True)
_TRIM_CODON_MASK = (False, False, False)


class _GroupEmitters(list):
    __slots__ = (
        "uniform_emitter",
        "uniform_selectors",
        "uniform_prefix_emitter",
        "uniform_prefix_selectors",
        "uniform_prefix_length",
    )

    def __init__(self):
        super().__init__()
        self.uniform_emitter = None
        self.uniform_selectors = None
        self.uniform_prefix_emitter = None
        self.uniform_prefix_selectors = None
        self.uniform_prefix_length = 0


class DNAThreader(Alignment):
    """
    Threads DNA on top of protein alignment
    """
    GAP_CHARS = frozenset({'-', '?', '*', 'X', 'x'})

    def __init__(self, args) -> None:
        self.process_args(args)

    def process_args(self, args):
        self.remove_stop_codon = args.stop
        self.protein_file_path = args.protein
        self.nucleotide_file_path = args.nucleotide
        self.clipkit_log_file = args.clipkit_log_file
        self.json_output = getattr(args, "json", False)

    @property
    def clipkit_log_data(self) -> list[list[str]]:
        if self.clipkit_log_file:
            with open(self.clipkit_log_file) as f:
                return [line.rstrip("\n").split(" ") for line in f]
        return None

    def run(self) -> None:
        prot_records = SeqIO.parse(self.protein_file_path, "fasta")
        pal2nal = self.thread(prot_records)

        if self.json_output:
            rows = [
                dict(taxon=gene_id, sequence=sequence)
                for gene_id, sequence in pal2nal.items()
            ]
            print_json(
                dict(
                    remove_stop_codon=self.remove_stop_codon,
                    clipkit_log_file=self.clipkit_log_file,
                    rows=rows,
                    taxa=rows,
                )
            )
            return

        if pal2nal:
            print("\n".join(
                f">{gene_id}\n{sequence}"
                for gene_id, sequence in pal2nal.items()
            ))

    def create_mask(self, length: int) -> list[bool]:
        return self._create_mask_and_all_true(length)[0]

    def _create_mask_and_all_true(self, length: int) -> tuple[list[bool], bool]:
        if not self.clipkit_log_file:
            return [True] * length, True

        mask = []
        extend = mask.extend
        all_sites_kept = True
        with open(self.clipkit_log_file, "rb") as f:
            for line in f:
                status_start = line.find(b" ") + 1
                status_end = status_start + 4
                keep = line.startswith(b"keep", status_start) and (
                    len(line) == status_end or line[status_end] in (32, 10)
                )
                if keep:
                    extend(_KEEP_CODON_MASK)
                else:
                    all_sites_kept = False
                    extend(_TRIM_CODON_MASK)
        return mask, all_sites_kept

    def normalize_p_seq(self, p_seq: "Seq") -> str:
        # triplicate each amino acid
        return ''.join(c * 3 for c in p_seq)

    def normalize_n_seq(self, n_seq: "Seq", p_seq: "Seq") -> str:
        # Pre-split codons for faster access
        n_seq_str = str(n_seq)
        codons = [n_seq_str[i:i+3] for i in range(0, len(n_seq_str), 3)]
        normalized_n_seq = []
        gap_chars = self.GAP_CHARS

        codon_idx = 0
        for aa in p_seq:
            if aa in gap_chars:
                normalized_n_seq.append("---")
            else:
                if codon_idx < len(codons):
                    normalized_n_seq.append(codons[codon_idx])
                    codon_idx += 1
                else:
                    normalized_n_seq.append("---")  # fallback in case of misalignment

        return ''.join(normalized_n_seq)

    def _thread_sequence(
        self,
        p_seq: "Seq" | str,
        n_seq: "Seq" | str,
        keep_mask: list[bool],
        keep_mask_all_true: bool = False,
        keep_mask_plan=None,
    ) -> str:
        gap_chars = self.GAP_CHARS
        p_seq_str = str(p_seq)
        n_seq_str = str(n_seq)

        if len(n_seq_str) % 3:
            normalized_p_seq = self.normalize_p_seq(p_seq_str)
            normalized_n_seq = self.normalize_n_seq(n_seq_str, normalized_p_seq)

            limit = min(
                len(normalized_p_seq) * 3,
                len(keep_mask),
                len(normalized_n_seq),
            )
            keep_mask_trimmed = keep_mask[:limit]
            result = []

            for idx in range(limit):
                if (
                    keep_mask_trimmed[idx]
                    and normalized_p_seq[idx // 3] not in gap_chars
                ):
                    result.append(normalized_n_seq[idx])
                else:
                    result.append('-')

            if self.remove_stop_codon and p_seq_str[-1] == "*":
                kept_indices = [
                    idx
                    for idx, should_keep in enumerate(keep_mask_trimmed)
                    if should_keep
                ]
                if len(kept_indices) >= 3:
                    for idx in kept_indices[-3:]:
                        if idx < len(result) and idx < len(normalized_n_seq):
                            result[idx] = normalized_n_seq[idx]

            return ''.join(
                result[idx]
                for idx, should_keep in enumerate(keep_mask_trimmed)
                if should_keep
            )

        limit = min(len(p_seq_str) * 9, len(keep_mask))
        normalized_p_chars = (limit + 2) // 3
        full_repeats, extra = divmod(normalized_p_chars, 3)
        if keep_mask_plan is not None and keep_mask_plan[0] == limit:
            all_sites_kept = keep_mask_all_true or keep_mask_plan[1]
            chunk_keep_bits = keep_mask_plan[2]
            grouped_keep_bits = keep_mask_plan[3]
            group_emitters = keep_mask_plan[4] if len(keep_mask_plan) > 4 else None
        else:
            all_sites_kept = keep_mask_all_true or (
                all(keep_mask) if len(keep_mask) == limit else all(keep_mask[:limit])
            )
            chunk_keep_bits = None
            grouped_keep_bits = None
            group_emitters = None
        inspected_aa_count = full_repeats + (1 if extra else 0)
        inspected_p_seq = p_seq_str[:inspected_aa_count]
        has_gap_chars = self._has_gap_char(inspected_p_seq)
        if (
            p_seq_str
            and all_sites_kept
            and len(n_seq_str) >= limit
            and not has_gap_chars
        ):
            return n_seq_str[:limit]

        restore_stop = self.remove_stop_codon and p_seq_str and p_seq_str[-1] == "*"
        if not restore_stop:
            if len(n_seq_str) >= limit and not has_gap_chars:
                return ''.join(compress(n_seq_str[:limit], keep_mask))

            if all_sites_kept:
                return self._thread_all_sites_kept(
                    p_seq_str,
                    n_seq_str,
                    limit,
                )

            threaded = []
            append = threaded.append
            codon_idx = 0
            codon_count = len(n_seq_str) // 3
            reachable_chunks = (limit + 2) // 3
            if group_emitters is not None:
                uniform_selectors = getattr(
                    group_emitters,
                    "uniform_selectors",
                    None,
                )
                if uniform_selectors is not None:
                    return self._thread_uniform_emit_plan(
                        p_seq_str,
                        n_seq_str,
                        group_emitters,
                        uniform_selectors,
                    )
                uniform_prefix_selectors = getattr(
                    group_emitters,
                    "uniform_prefix_selectors",
                    None,
                )
                if uniform_prefix_selectors is not None:
                    return self._thread_uniform_prefix_emit_plan(
                        p_seq_str,
                        n_seq_str,
                        group_emitters,
                        uniform_prefix_selectors,
                    )
                self._append_group_emitters(
                    p_seq_str,
                    n_seq_str,
                    group_emitters,
                    threaded,
                )
            elif grouped_keep_bits is not None:
                for aa, bit_group in zip(p_seq_str, grouped_keep_bits):
                    if aa in gap_chars:
                        for bits in bit_group:
                            if bits & 1:
                                append('-')
                            if bits & 2:
                                append('-')
                            if bits & 4:
                                append('-')
                    else:
                        for bits in bit_group:
                            if codon_idx >= codon_count:
                                c0 = c1 = c2 = '-'
                            else:
                                start = codon_idx * 3
                                c0, c1, c2 = n_seq_str[start:start + 3]
                                codon_idx += 1
                            if bits & 1:
                                append(c0)
                            if bits & 2:
                                append(c1)
                            if bits & 4:
                                append(c2)
            else:
                idx = 0
                for chunk_idx in range(reachable_chunks):
                    aa = p_seq_str[chunk_idx // 3]
                    if aa in gap_chars or codon_idx >= codon_count:
                        c0 = c1 = c2 = '-'
                    else:
                        start = codon_idx * 3
                        c0, c1, c2 = n_seq_str[start:start + 3]
                        codon_idx += 1

                    if keep_mask[idx]:
                        append(c0)
                    idx += 1
                    if idx < limit and keep_mask[idx]:
                        append(c1)
                    idx += 1
                    if idx < limit and keep_mask[idx]:
                        append(c2)
                    idx += 1

            return ''.join(threaded)

        if all_sites_kept:
            return self._thread_all_sites_kept(
                p_seq_str,
                n_seq_str,
                limit,
            )

        codon_idx = 0
        codon_count = len(n_seq_str) // 3
        reachable_chunks = (limit + 2) // 3
        threaded = []
        append = threaded.append
        idx = 0

        if group_emitters is not None:
            uniform_selectors = getattr(
                group_emitters,
                "uniform_selectors",
                None,
            )
            if uniform_selectors is not None:
                return self._thread_uniform_emit_plan(
                    p_seq_str,
                    n_seq_str,
                    group_emitters,
                    uniform_selectors,
                )
            uniform_prefix_selectors = getattr(
                group_emitters,
                "uniform_prefix_selectors",
                None,
            )
            if uniform_prefix_selectors is not None:
                return self._thread_uniform_prefix_emit_plan(
                    p_seq_str,
                    n_seq_str,
                    group_emitters,
                    uniform_prefix_selectors,
                )
            self._append_group_emitters(
                p_seq_str,
                n_seq_str,
                group_emitters,
                threaded,
            )
        elif grouped_keep_bits is not None:
            for aa, bit_group in zip(p_seq_str, grouped_keep_bits):
                if aa in gap_chars:
                    for bits in bit_group:
                        if bits & 1:
                            append('-')
                        if bits & 2:
                            append('-')
                        if bits & 4:
                            append('-')
                else:
                    for bits in bit_group:
                        if codon_idx >= codon_count:
                            c0 = c1 = c2 = '-'
                        else:
                            start = codon_idx * 3
                            c0, c1, c2 = n_seq_str[start:start + 3]
                            codon_idx += 1
                        if bits & 1:
                            append(c0)
                        if bits & 2:
                            append(c1)
                        if bits & 4:
                            append(c2)
        else:
            for chunk_idx in range(reachable_chunks):
                aa = p_seq_str[chunk_idx // 3]
                if aa in gap_chars or codon_idx >= codon_count:
                    c0 = c1 = c2 = '-'
                else:
                    start = codon_idx * 3
                    c0, c1, c2 = n_seq_str[start:start + 3]
                    codon_idx += 1

                if keep_mask[idx]:
                    append(c0)
                idx += 1
                if idx < limit and keep_mask[idx]:
                    append(c1)
                idx += 1
                if idx < limit and keep_mask[idx]:
                    append(c2)
                idx += 1

        return ''.join(threaded)

    def _append_group_emitters(
        self,
        p_seq_str: str,
        n_seq_str: str,
        group_emitters,
        threaded: list[str],
    ) -> None:
        append = threaded.append
        gap_chars = self.GAP_CHARS
        nucleotide_offset = 0
        join = ''.join

        for aa, emitter in zip(p_seq_str, group_emitters):
            keep_indices, getter, gap_fill = emitter
            if aa in gap_chars:
                append(gap_fill)
                continue

            if getter is None:
                nucleotide_offset += 9
                continue
            segment = n_seq_str[nucleotide_offset:nucleotide_offset + 9]
            nucleotide_offset += 9
            if len(segment) == 9:
                if getter is _FULL_SEGMENT:
                    append(segment)
                    continue
                emitted = getter(segment)
                append(join(emitted) if isinstance(emitted, tuple) else emitted)
            else:
                append(join(
                    segment[index] if index < len(segment) else '-'
                    for index in keep_indices
                ))

    def _thread_uniform_emit_plan(
        self,
        p_seq_str: str,
        n_seq_str: str,
        group_emitters,
        selectors,
    ) -> str:
        emitter = group_emitters.uniform_emitter
        keep_indices, getter, gap_fill = emitter
        group_count = min(len(p_seq_str), len(group_emitters))
        threaded = []
        append = threaded.append
        gap_chars = self.GAP_CHARS
        nucleotide_offset = 0
        aa_idx = 0
        join = ''.join

        while aa_idx < group_count:
            is_gap_run = p_seq_str[aa_idx] in gap_chars
            run_start = aa_idx
            aa_idx += 1
            while (
                aa_idx < group_count
                and (p_seq_str[aa_idx] in gap_chars) == is_gap_run
            ):
                aa_idx += 1

            run_groups = aa_idx - run_start
            if is_gap_run:
                if gap_fill:
                    append(gap_fill * run_groups)
                continue

            emit_len = run_groups * 9
            segment = n_seq_str[nucleotide_offset:nucleotide_offset + emit_len]
            nucleotide_offset += emit_len
            if getter is None:
                continue
            if len(segment) == emit_len:
                if getter is _FULL_SEGMENT:
                    append(segment)
                else:
                    append(join(compress(segment, selectors * run_groups)))
                continue

            for group_offset in range(run_groups):
                base_idx = group_offset * 9
                append(join(
                    segment[base_idx + index]
                    if base_idx + index < len(segment)
                    else '-'
                    for index in keep_indices
                ))

        return join(threaded)

    def _thread_uniform_prefix_emit_plan(
        self,
        p_seq_str: str,
        n_seq_str: str,
        group_emitters,
        selectors,
    ) -> str:
        emitter = group_emitters.uniform_prefix_emitter
        keep_indices, getter, gap_fill = emitter
        group_count = min(len(p_seq_str), len(group_emitters))
        prefix_count = min(group_count, group_emitters.uniform_prefix_length)
        threaded = []
        append = threaded.append
        gap_chars = self.GAP_CHARS
        nucleotide_offset = 0
        aa_idx = 0
        join = ''.join

        while aa_idx < prefix_count:
            is_gap_run = p_seq_str[aa_idx] in gap_chars
            run_start = aa_idx
            aa_idx += 1
            while (
                aa_idx < prefix_count
                and (p_seq_str[aa_idx] in gap_chars) == is_gap_run
            ):
                aa_idx += 1

            run_groups = aa_idx - run_start
            if is_gap_run:
                if gap_fill:
                    append(gap_fill * run_groups)
                continue

            emit_len = run_groups * 9
            segment = n_seq_str[nucleotide_offset:nucleotide_offset + emit_len]
            nucleotide_offset += emit_len
            if getter is None:
                continue
            if len(segment) == emit_len:
                if getter is _FULL_SEGMENT:
                    append(segment)
                else:
                    append(join(compress(segment, selectors * run_groups)))
                continue

            for group_offset in range(run_groups):
                base_idx = group_offset * 9
                append(join(
                    segment[base_idx + index]
                    if base_idx + index < len(segment)
                    else '-'
                    for index in keep_indices
                ))

        for group_idx in range(prefix_count, group_count):
            aa = p_seq_str[group_idx]
            keep_indices, getter, gap_fill = group_emitters[group_idx]
            if aa in gap_chars:
                append(gap_fill)
                continue

            if getter is None:
                nucleotide_offset += 9
                continue
            segment = n_seq_str[nucleotide_offset:nucleotide_offset + 9]
            nucleotide_offset += 9
            if len(segment) == 9:
                if getter is _FULL_SEGMENT:
                    append(segment)
                    continue
                emitted = getter(segment)
                append(join(emitted) if isinstance(emitted, tuple) else emitted)
            else:
                append(join(
                    segment[index] if index < len(segment) else '-'
                    for index in keep_indices
                ))

        return join(threaded)

    @staticmethod
    def _has_gap_char(seq: str) -> bool:
        return (
            "-" in seq
            or "?" in seq
            or "*" in seq
            or "X" in seq
            or "x" in seq
        )

    def _thread_all_sites_kept(
        self,
        p_seq_str: str,
        n_seq_str: str,
        limit: int,
    ) -> str:
        threaded = []
        append = threaded.append
        gap_chars = self.GAP_CHARS
        nucleotide_offset = 0
        group_count = (limit + 8) // 9
        aa_idx = 0

        while aa_idx < group_count:
            is_gap_run = p_seq_str[aa_idx] in gap_chars
            run_start = aa_idx
            aa_idx += 1
            while (
                aa_idx < group_count
                and (p_seq_str[aa_idx] in gap_chars) == is_gap_run
            ):
                aa_idx += 1

            run_groups = aa_idx - run_start
            emit_len = min(run_groups * 9, limit - run_start * 9)
            if is_gap_run:
                append('-' * emit_len)
                continue

            segment = n_seq_str[nucleotide_offset:nucleotide_offset + emit_len]
            nucleotide_offset += run_groups * 9
            if len(segment) == emit_len:
                append(segment)
            else:
                append(segment + '-' * (emit_len - len(segment)))

        return ''.join(threaded)

    @staticmethod
    def _create_thread_mask_plan(keep_mask: list[bool], protein_length: int):
        limit = min(protein_length * 9, len(keep_mask))
        all_sites_kept = all(keep_mask) if len(keep_mask) == limit else all(keep_mask[:limit])
        reachable_chunks = (limit + 2) // 3
        chunk_keep_bits = []
        for chunk_idx in range(reachable_chunks):
            idx = chunk_idx * 3
            bits = 0
            if keep_mask[idx]:
                bits |= 1
            if idx + 1 < limit and keep_mask[idx + 1]:
                bits |= 2
            if idx + 2 < limit and keep_mask[idx + 2]:
                bits |= 4
            chunk_keep_bits.append(bits)
        grouped_keep_bits = [
            tuple(chunk_keep_bits[idx:idx + 3])
            for idx in range(0, len(chunk_keep_bits), 3)
        ]
        group_emitters = _GroupEmitters()
        for bit_group in grouped_keep_bits:
            if bit_group == (7, 7, 7):
                keep_indices = (0, 1, 2, 3, 4, 5, 6, 7, 8)
                getter = _FULL_SEGMENT
            elif not any(bit_group):
                keep_indices = ()
                getter = None
            else:
                keep_indices_list = []
                for chunk_offset, bits in enumerate(bit_group):
                    base_idx = chunk_offset * 3
                    if bits & 1:
                        keep_indices_list.append(base_idx)
                    if bits & 2:
                        keep_indices_list.append(base_idx + 1)
                    if bits & 4:
                        keep_indices_list.append(base_idx + 2)
                keep_indices = tuple(keep_indices_list)
                if len(keep_indices) == 1:
                    getter = itemgetter(keep_indices[0])
                else:
                    getter = itemgetter(*keep_indices)
            group_emitters.append((keep_indices, getter, '-' * len(keep_indices)))
        if group_emitters:
            first_emitter = group_emitters[0]
            if all(emitter[0] == first_emitter[0] for emitter in group_emitters):
                keep_indices = first_emitter[0]
                keep_index_set = set(keep_indices)
                group_emitters.uniform_emitter = first_emitter
                group_emitters.uniform_selectors = tuple(
                    index in keep_index_set for index in range(9)
                )
            elif len(group_emitters) > 1:
                prefix_length = len(group_emitters) - 1
                if all(
                    emitter[0] == first_emitter[0]
                    for emitter in group_emitters[:prefix_length]
                ):
                    keep_indices = first_emitter[0]
                    keep_index_set = set(keep_indices)
                    group_emitters.uniform_prefix_emitter = first_emitter
                    group_emitters.uniform_prefix_selectors = tuple(
                        index in keep_index_set for index in range(9)
                    )
                    group_emitters.uniform_prefix_length = prefix_length
        return (
            limit,
            all_sites_kept,
            chunk_keep_bits,
            grouped_keep_bits,
            group_emitters,
        )

    def thread(self, prot_records) -> dict[str, str]:
        pal2nal = dict()
        prot_dict = SeqIO.to_dict(prot_records)

        if not prot_dict:
            print("Protein file is empty or incorrectly formatted.")
            sys.exit(2)

        # Pre-load nucleotide sequences only for proteins we have
        nucl_records = {}
        for record in SeqIO.parse(self.nucleotide_file_path, "fasta"):
            if record.id in prot_dict:
                nucl_records[record.id] = record

        length = len(next(iter(prot_dict.values())).seq)
        mask_length = length * 3
        if (
            self.__class__.create_mask is DNAThreader.create_mask
            and "create_mask" not in self.__dict__
        ):
            keep_mask, keep_mask_all_true = self._create_mask_and_all_true(
                mask_length
            )
        else:
            keep_mask = self.create_mask(mask_length)
            keep_mask_all_true = self.clipkit_log_file is None or all(keep_mask)
        keep_mask_plan = (
            None
            if keep_mask_all_true
            else self._create_thread_mask_plan(keep_mask, length)
        )

        for gene_id, protein_seq_record in prot_dict.items():
            try:
                if gene_id not in nucl_records:
                    print(f"Nucleotide sequence for {gene_id} not found.")
                    sys.exit(2)

                p_seq = protein_seq_record.seq
                n_seq = nucl_records[gene_id].seq

                pal2nal[gene_id] = self._thread_sequence(
                    p_seq,
                    n_seq,
                    keep_mask,
                    keep_mask_all_true=keep_mask_all_true,
                    keep_mask_plan=keep_mask_plan,
                )

            except KeyError:
                print(f"Nucleotide sequence for {gene_id} not found.")
                sys.exit(2)

        return pal2nal
