def _clean_sequence(sequence_parts: list[str]) -> str:
    if len(sequence_parts) == 1:
        sequence = sequence_parts[0]
    else:
        sequence = "".join(sequence_parts)
    if " " in sequence:
        sequence = sequence.replace(" ", "")
    if "\r" in sequence:
        sequence = sequence.replace("\r", "")
    return sequence


def _clean_sequence_bytes(sequence_parts: list[bytes]) -> str:
    if len(sequence_parts) == 1:
        sequence = sequence_parts[0]
    else:
        sequence = b"".join(sequence_parts)
    if b" " in sequence:
        sequence = sequence.replace(b" ", b"")
    if b"\r" in sequence:
        sequence = sequence.replace(b"\r", b"")
    return sequence.decode()


def _clean_upper_sequence(sequence_parts: list[str]) -> str:
    if len(sequence_parts) == 1:
        sequence = sequence_parts[0]
        if " " in sequence:
            sequence = sequence.replace(" ", "")
        if "\r" in sequence:
            sequence = sequence.replace("\r", "")
        return sequence.upper()
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
    wanted = {entry.encode(): entry for entry in entries}
    records = {}
    seen = set()
    clean_sequence = _clean_sequence_bytes
    with open(path, "rb") as handle:
        record_id = None
        wanted_record = False
        sequence_parts = []
        for line in handle:
            if line and line[0] == 62:
                if record_id is not None and wanted_record:
                    records[record_id] = clean_sequence(sequence_parts)
                token = line[1:].split(None, 1)[0]
                if token in seen:
                    raise ValueError(f"Duplicate key '{token.decode()}'")
                seen.add(token)
                record_id = wanted.get(token)
                wanted_record = record_id is not None
                sequence_parts = []
            elif record_id is not None and wanted_record:
                sequence_parts.append(line.rstrip())

        if record_id is not None and wanted_record:
            records[record_id] = clean_sequence(sequence_parts)

    for entry in entries:
        if entry not in records:
            raise KeyError(entry)
    return records


def read_fasta_first_tokens(path: str) -> list[str]:
    taxa = []
    append = taxa.append
    with open(path, "rb") as handle:
        for line in handle:
            if line[0] == 62:
                append(line[1:].split(None, 1)[0].decode())
    return taxa


def read_fasta_first_token_set(path: str) -> set[str]:
    taxa = set()
    with open(path, "rb") as handle:
        for line in handle:
            if line and line[0] == 62:
                taxa.add(line[1:].split(None, 1)[0].decode())
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
