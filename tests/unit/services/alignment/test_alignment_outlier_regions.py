from types import SimpleNamespace

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pytest

from phykit.errors import PhykitUserError
from phykit.services.alignment.alignment_outlier_regions import (
    AlignmentOutlierRegions,
    _column_divergence_scores,
    _jenks_two_class_cutoff,
    _remove_long_regions,
    _segment_outlier_windows,
    _upper_tail_boundary,
)


def _service(**overrides):
    values = dict(
        alignment="alignment.fa",
        cutoff=3.0,
        mask_output=None,
        mask_character="?",
        report=None,
        json=False,
    )
    values.update(overrides)
    return AlignmentOutlierRegions(SimpleNamespace(**values))


def _local_error_alignment():
    records = [
        SeqRecord(Seq("A" * 120), id=f"normal_{index}")
        for index in range(9)
    ]
    records.append(
        SeqRecord(
            Seq(("A" * 50) + ("C" * 10) + ("A" * 60)),
            id="local_error",
        )
    )
    return MultipleSeqAlignment(records)


def test_column_divergence_scores_use_column_frequency():
    matrix = np.frombuffer(b"AAAC", dtype=np.uint8).reshape(4, 1)
    valid = np.ones_like(matrix, dtype=np.bool_)

    observed = _column_divergence_scores(matrix, valid)

    np.testing.assert_allclose(observed[:, 0], [2 / 3, 2 / 3, 2 / 3, 2])


def test_column_divergence_scores_ignore_missing_symbols():
    matrix = np.frombuffer(b"AA-N", dtype=np.uint8).reshape(4, 1)
    valid = np.array([[True], [True], [False], [False]])

    observed = _column_divergence_scores(matrix, valid)

    np.testing.assert_allclose(observed[:, 0], [1, 1, 0, 0])


def test_jenks_cutoff_separates_high_scoring_windows():
    observed = _jenks_two_class_cutoff(np.array([1, 1, 1, 5, 5]))

    assert observed == 1.0


def test_upper_tail_boundary_matches_tail_fraction():
    observed = _upper_tail_boundary(np.array([1, 2, 3, 4]), 0.25)

    assert observed == 3.0


def test_window_segmentation_smooths_a_contiguous_outlier_run():
    windows = np.zeros(60, dtype=np.bool_)
    windows[20:35] = True

    observed = _segment_outlier_windows(windows, window_size=5)

    assert np.flatnonzero(observed).tolist() == list(range(20, 35))


def test_long_regions_are_removed_at_short_window_scales():
    mask = np.zeros(50, dtype=np.bool_)
    mask[5:12] = True
    mask[20:41] = True

    observed = _remove_long_regions(mask, maximum_length=10)

    assert np.flatnonzero(observed).tolist() == list(range(5, 12))


def test_detect_finds_only_the_local_error_sequence():
    result = _service().detect(_local_error_alignment(), is_protein=False)

    assert result["affected_taxa"] == 1
    assert result["masked_residues"] == 10
    assert result["regions"] == [
        {
            "taxon": "local_error",
            "alignment_start": 49,
            "alignment_end": 58,
            "sequence_start": 49,
            "sequence_end": 58,
            "length": 10,
            "mean_divergence_score": 4.2,
            "max_divergence_score": 5.0,
            "window_sizes": [5],
        }
    ]
    assert result["masked_sequences"][0] == "A" * 120
    assert result["masked_sequences"][-1].count("?") == 10


def test_detect_returns_no_regions_for_identical_sequences():
    alignment = MultipleSeqAlignment(
        [SeqRecord(Seq("ACGT" * 30), id=f"taxon_{index}") for index in range(10)]
    )

    result = _service().detect(alignment, is_protein=False)

    assert result["regions"] == []
    assert result["masked_residues"] == 0
    assert result["affected_taxa"] == 0


def test_detect_preserves_gaps_and_ambiguous_symbols_when_masking():
    alignment = _local_error_alignment()
    gapped = list(str(alignment[-1].seq))
    gapped[10] = "-"
    gapped[20] = "N"
    alignment[-1].seq = Seq("".join(gapped))

    result = _service(mask_character="X").detect(alignment, is_protein=False)

    assert result["masked_sequences"][-1][10] == "-"
    assert result["masked_sequences"][-1][20] == "N"
    assert result["masked_sequences"][-1].count("X") > 0


@pytest.mark.parametrize(
    "overrides,message",
    [
        ({"cutoff": 1.0}, "cutoff must be greater than 1.0"),
        ({"mask_character": "XX"}, "exactly one"),
        ({"mask_character": "-"}, "cannot be '-'"),
        ({"mask_output": "alignment.fa"}, "cannot overwrite"),
        (
            {"mask_output": "same.out", "report": "same.out"},
            "different paths",
        ),
    ],
)
def test_validate_options_rejects_invalid_values(overrides, message):
    service = _service(**overrides)

    with pytest.raises(PhykitUserError) as error:
        service._validate_options()

    assert message in " ".join(error.value.messages)


def test_render_tsv_includes_coordinate_contract():
    result = _service().detect(_local_error_alignment(), is_protein=False)

    report = AlignmentOutlierRegions._render_tsv(result["regions"])

    assert report.splitlines()[0].startswith("taxon\talignment_start")
    assert report.splitlines()[1].startswith("local_error\t49\t58\t49\t58")

