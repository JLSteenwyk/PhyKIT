import pytest
from argparse import Namespace

from phykit.services.alignment.mask_alignment import MaskAlignment
import phykit.services.alignment.mask_alignment as mask_alignment_module


@pytest.fixture
def args():
    kwargs = dict(
        alignment="/some/path/to/file.fa",
        max_gap=1.0,
        min_occupancy=0.0,
        max_entropy=None,
    )
    return Namespace(**kwargs)


class TestMaskAlignment(object):
    def test_init_sets_alignment_file_path(self, args):
        masker = MaskAlignment(args)
        assert masker.alignment_file_path == args.alignment

    def test_keep_mask_by_gap_and_occupancy(self, alignment_simple, args):
        args.max_gap = 0.3
        args.min_occupancy = 0.8
        masker = MaskAlignment(args)
        keep_mask = masker.calculate_keep_mask(alignment_simple, is_protein=False)
        assert keep_mask.tolist() == [True, False, True, False, True, True]

    def test_keep_mask_by_entropy(self, alignment_simple, args):
        args.max_entropy = 0.5
        masker = MaskAlignment(args)
        keep_mask = masker.calculate_keep_mask(alignment_simple, is_protein=False)
        assert keep_mask.tolist() == [True, False, False, True, False, False]

    def test_process_args_json_default(self):
        parsed = MaskAlignment(
            Namespace(alignment="x.fa", max_gap=0.8, min_occupancy=0.2, max_entropy=None)
        ).process_args(
            Namespace(alignment="x.fa", max_gap=0.8, min_occupancy=0.2, max_entropy=None)
        )
        assert parsed["json_output"] is False

    def test_validate_thresholds_invalid_values_exit(self, capsys):
        with pytest.raises(SystemExit) as exc:
            MaskAlignment(
                Namespace(alignment="x.fa", max_gap=1.2, min_occupancy=0.0, max_entropy=None)
            )._validate_thresholds()
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "max_gap must be between 0 and 1." in out

        with pytest.raises(SystemExit) as exc:
            MaskAlignment(
                Namespace(alignment="x.fa", max_gap=0.8, min_occupancy=-0.1, max_entropy=None)
            )._validate_thresholds()
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "min_occupancy must be between 0 and 1." in out

        with pytest.raises(SystemExit) as exc:
            MaskAlignment(
                Namespace(alignment="x.fa", max_gap=0.8, min_occupancy=0.1, max_entropy=-0.5)
            )._validate_thresholds()
        assert exc.value.code == 2
        out, _ = capsys.readouterr()
        assert "max_entropy must be >= 0." in out

    def test_apply_mask(self, alignment_simple, args):
        masker = MaskAlignment(args)
        masked = masker.apply_mask(alignment_simple, keep_mask=mask_alignment_module.np.array([True, False, True, False, True, False]))
        assert set(masked.keys()) == {record.id for record in alignment_simple}
        assert all(len(seq) == 3 for seq in masked.values())

    def test_run_json_output(self, alignment_simple, mocker):
        masker = MaskAlignment(
            Namespace(alignment="x.fa", max_gap=1.0, min_occupancy=0.0, max_entropy=None, json=True)
        )
        mocker.patch.object(MaskAlignment, "get_alignment_and_format", return_value=(alignment_simple, "fasta", False))
        mocker.patch.object(MaskAlignment, "calculate_keep_mask", return_value=mask_alignment_module.np.array([True, True, False, False, True, True]))
        mocker.patch.object(MaskAlignment, "apply_mask", return_value={"a": "ACGT", "b": "ACGT"})
        mocked_json = mocker.patch("phykit.services.alignment.mask_alignment.print_json")

        masker.run()

        payload = mocked_json.call_args.args[0]
        assert payload["kept_sites"] == 4
        assert payload["total_sites"] == 6
        assert payload["rows"] == payload["taxa"]

    def test_run_prints_fasta(self, alignment_simple, capsys, mocker):
        masker = MaskAlignment(
            Namespace(alignment="x.fa", max_gap=1.0, min_occupancy=0.0, max_entropy=None, json=False)
        )
        mocker.patch.object(MaskAlignment, "get_alignment_and_format", return_value=(alignment_simple, "fasta", False))
        mocker.patch.object(MaskAlignment, "calculate_keep_mask", return_value=mask_alignment_module.np.array([True, False, True, False, True, False]))
        mocker.patch.object(MaskAlignment, "apply_mask", return_value={"a": "AAA", "b": "CCC"})

        masker.run()

        out, _ = capsys.readouterr()
        assert ">a\nAAA" in out
        assert ">b\nCCC" in out
