def _clean_sequence(sequence_parts: list[str]) -> str:
    sequence = "".join(sequence_parts)
    if " " in sequence:
        sequence = sequence.replace(" ", "")
    if "\r" in sequence:
        sequence = sequence.replace("\r", "")
    return sequence


def _clean_upper_sequence(sequence_parts: list[str]) -> str:
    return _clean_sequence(sequence_parts).upper()


def read_fasta_first_token_upper(path: str) -> dict[str, str]:
    sequences = {}
    with open(path) as handle:
        record_id = None
        sequence_parts = []
        for line in handle:
            if line[0] == ">":
                if record_id is not None:
                    sequences[record_id] = _clean_upper_sequence(sequence_parts)
                record_id = line[1:].split(None, 1)[0]
                sequence_parts = []
            elif record_id is not None:
                sequence_parts.append(line.rstrip())

        if record_id is not None:
            sequences[record_id] = _clean_upper_sequence(sequence_parts)
    return sequences


def read_fasta_first_token(path: str) -> dict[str, str]:
    sequences = {}
    with open(path) as handle:
        record_id = None
        sequence_parts = []
        for line in handle:
            if line[0] == ">":
                if record_id is not None:
                    sequences[record_id] = _clean_sequence(sequence_parts)
                record_id = line[1:].split(None, 1)[0]
                sequence_parts = []
            elif record_id is not None:
                sequence_parts.append(line.rstrip())

        if record_id is not None:
            sequences[record_id] = _clean_sequence(sequence_parts)
    return sequences


def read_fasta_first_token_records(path: str, record_factory) -> tuple[set[str], list]:
    records = []
    taxa = set()
    records_append = records.append
    taxa_add = taxa.add
    clean_sequence = _clean_sequence

    with open(path) as handle:
        record_id = None
        sequence_parts = []
        for line in handle:
            if line[0] == ">":
                if record_id is not None:
                    records_append(
                        record_factory(record_id, clean_sequence(sequence_parts))
                    )
                record_id = line[1:].split(None, 1)[0]
                taxa_add(record_id)
                sequence_parts = []
            elif record_id is not None:
                sequence_parts.append(line.rstrip())

        if record_id is not None:
            records_append(
                record_factory(record_id, clean_sequence(sequence_parts))
            )
    return taxa, records


def read_unique_fasta_first_token(path: str) -> dict[str, str]:
    sequences = {}
    with open(path) as handle:
        record_id = None
        sequence_parts = []
        for line in handle:
            if line[0] == ">":
                if record_id is not None:
                    sequences[record_id] = _clean_sequence(sequence_parts)
                record_id = line[1:].split(None, 1)[0]
                if record_id in sequences:
                    raise ValueError(f"Duplicate key {record_id!r}")
                sequence_parts = []
            elif record_id is not None:
                sequence_parts.append(line.rstrip())

        if record_id is not None:
            sequences[record_id] = _clean_sequence(sequence_parts)
    return sequences


def read_unique_fasta_entries(
    path: str,
    entries: list[str],
) -> dict[str, str]:
    wanted = set(entries)
    records = {}
    seen = set()
    with open(path) as handle:
        record_id = None
        wanted_record = False
        sequence_parts = []
        for line in handle:
            if line[0] == ">":
                if record_id is not None and wanted_record:
                    records[record_id] = _clean_sequence(sequence_parts)
                record_id = line[1:].split(None, 1)[0]
                if record_id in seen:
                    raise ValueError(f"Duplicate key '{record_id}'")
                seen.add(record_id)
                wanted_record = record_id in wanted
                sequence_parts = []
            elif record_id is not None and wanted_record:
                sequence_parts.append(line.rstrip())

        if record_id is not None and wanted_record:
            records[record_id] = _clean_sequence(sequence_parts)

    for entry in entries:
        if entry not in records:
            raise KeyError(entry)
    return records


def read_fasta_first_tokens(path: str) -> list[str]:
    taxa = []
    with open(path) as handle:
        for line in handle:
            if line[0] == ">":
                taxa.append(line[1:].split(None, 1)[0])
    return taxa


def read_fasta_first_token_set(path: str) -> set[str]:
    taxa = set()
    with open(path) as handle:
        for line in handle:
            if line[0] == ">":
                taxa.add(line[1:].split(None, 1)[0])
    return taxa


def read_fasta_first_token_upper_with_taxa(
    path: str,
) -> tuple[list[str], dict[str, str]]:
    taxa = []
    sequences = {}
    with open(path) as handle:
        record_id = None
        sequence_parts = []
        for line in handle:
            if line[0] == ">":
                if record_id is not None:
                    sequences[record_id] = _clean_upper_sequence(sequence_parts)
                record_id = line[1:].split(None, 1)[0]
                taxa.append(record_id)
                sequence_parts = []
            elif record_id is not None:
                sequence_parts.append(line.rstrip())

        if record_id is not None:
            sequences[record_id] = _clean_upper_sequence(sequence_parts)
    return taxa, sequences
