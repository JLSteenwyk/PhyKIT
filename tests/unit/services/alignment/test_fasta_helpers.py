from phykit.services.alignment._fasta import _clean_sequence


def test_clean_sequence_single_line_preserves_clean_sequence():
    sequence = "ACGTACGT"

    assert _clean_sequence([sequence]) == sequence


def test_clean_sequence_removes_spaces_and_carriage_returns():
    assert _clean_sequence(["AC GT\r", "TA"]) == "ACGTTA"
