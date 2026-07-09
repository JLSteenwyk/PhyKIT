from phykit.services.alignment._fasta import (
    _clean_sequence,
    _clean_upper_sequence,
    read_unique_fasta_entries,
    read_fasta_first_tokens,
    read_fasta_first_token_set,
)


def test_clean_sequence_single_line_preserves_clean_sequence():
    sequence = "ACGTACGT"

    assert _clean_sequence([sequence]) == sequence


def test_clean_sequence_removes_spaces_and_carriage_returns():
    assert _clean_sequence(["AC GT\r", "TA"]) == "ACGTTA"


def test_clean_upper_sequence_single_line_cleans_and_uppercases():
    assert _clean_upper_sequence(["ac gt\r"]) == "ACGT"


def test_read_fasta_first_token_set_uses_first_header_tokens(tmp_path):
    fasta = tmp_path / "aln.fa"
    fasta.write_text(
        "\n"
        ">taxon_a description words\n"
        "ACGT\n"
        ">taxon_b\tother words\n"
        "TGCA\n"
    )

    assert read_fasta_first_token_set(str(fasta)) == {"taxon_a", "taxon_b"}


def test_read_fasta_first_tokens_uses_binary_header_scan(tmp_path):
    fasta = tmp_path / "aln.fa"
    fasta.write_bytes(
        b">taxon_a description words\n"
        b"ACGT\xff\n"
        b">taxon_b\tother words\n"
        b"TGCA\n"
    )

    assert read_fasta_first_tokens(str(fasta)) == ["taxon_a", "taxon_b"]


def test_read_unique_fasta_entries_cleans_requested_sequences(tmp_path):
    fasta = tmp_path / "aln.fa"
    fasta.write_bytes(
        b">taxon_a description words\n"
        b"AC GT\r\n"
        b"TA\n"
        b">taxon_b other words\n"
        b"GG\n"
    )

    assert read_unique_fasta_entries(str(fasta), ["taxon_a"]) == {
        "taxon_a": "ACGTTA",
    }


def test_read_unique_fasta_entries_scans_later_duplicates(tmp_path):
    fasta = tmp_path / "aln.fa"
    fasta.write_text(
        ">taxon_a\n"
        "ACGT\n"
        ">taxon_b\n"
        "GGGG\n"
        ">taxon_a duplicate\n"
        "TTTT\n"
    )

    try:
        read_unique_fasta_entries(str(fasta), ["taxon_a"])
    except ValueError as exc:
        assert "Duplicate key 'taxon_a'" in str(exc)
    else:
        raise AssertionError("expected duplicate requested entry to be detected")
