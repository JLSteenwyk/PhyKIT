from phykit.services.alignment._fasta import (
    _clean_sequence,
    read_fasta_first_tokens,
    read_fasta_first_token_set,
)


def test_clean_sequence_single_line_preserves_clean_sequence():
    sequence = "ACGTACGT"

    assert _clean_sequence([sequence]) == sequence


def test_clean_sequence_removes_spaces_and_carriage_returns():
    assert _clean_sequence(["AC GT\r", "TA"]) == "ACGTTA"


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
