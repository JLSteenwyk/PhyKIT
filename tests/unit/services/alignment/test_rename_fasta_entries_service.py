from argparse import Namespace
import subprocess
import sys

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.rename_fasta_entries import RenameFastaEntries


def test_module_import_does_not_import_biopython_seqio():
    code = """
import sys
import phykit.services.alignment.rename_fasta_entries
assert "typing" not in sys.modules
assert "Bio.SeqIO" not in sys.modules
assert "Bio.SeqIO.FastaIO" not in sys.modules
assert "Bio.AlignIO" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], check=True)


@pytest.fixture
def args():
    return Namespace(
        fasta="/some/path/to/file.fa",
        idmap="/some/path/to/idmap.txt",
        output=None,
    )


class TestRenameFastaEntries:
    def test_init_sets_expected_attrs(self, args):
        service = RenameFastaEntries(args)
        assert service.fasta == args.fasta
        assert service.idmap == args.idmap
        assert service.output_file_path == f"{args.fasta}.renamed.fa"
        assert service.json_output is False

    def test_process_args_honors_output_and_json(self):
        args = Namespace(
            fasta="/some/path/to/file.fa",
            idmap="/some/path/to/idmap.txt",
            output="/tmp/out.fa",
            json=True,
        )
        service = RenameFastaEntries(args)
        assert service.output_file_path == "/tmp/out.fa.renamed.fa"
        assert service.json_output is True

    def test_load_idmap_parses_file(self, tmp_path):
        idmap_file = tmp_path / "idmap.txt"
        idmap_file.write_text("a A\nb B\n")
        args = Namespace(fasta="/x.fa", idmap=str(idmap_file), output=None)
        service = RenameFastaEntries(args)
        assert service.load_idmap(str(idmap_file)) == {"a": "A", "b": "B"}

    def test_load_idmap_preserves_split_semantics(self, tmp_path):
        idmap_file = tmp_path / "idmap.txt"
        idmap_file.write_text("a\tA\nb B\na A2\n")
        args = Namespace(fasta="/x.fa", idmap=str(idmap_file), output=None)
        service = RenameFastaEntries(args)

        assert service.load_idmap(str(idmap_file)) == {"a": "A2", "b": "B"}

    def test_load_idmap_missing_exits(self):
        args = Namespace(fasta="/x.fa", idmap="/does/not/exist.txt", output=None)
        service = RenameFastaEntries(args)
        with pytest.raises(SystemExit) as excinfo:
            service.load_idmap(service.idmap)
        assert excinfo.value.code == 2

    def test_replace_ids_and_write(self, tmp_path, args):
        service = RenameFastaEntries(args)
        output_file = tmp_path / "renamed.fa"
        records = iter(
            [
                SeqRecord(Seq("ACGT"), id="a", description="a"),
                SeqRecord(Seq("TTTT"), id="b", description="b"),
            ]
        )
        renamed_count, total_records = service.replace_ids_and_write(
            str(output_file), records, {"a": "A"}
        )
        assert (renamed_count, total_records) == (1, 2)
        content = output_file.read_text()
        assert ">A" in content
        assert ">b" in content

    def test_replace_ids_and_write_preserves_fasta_formatting(self, tmp_path, args):
        service = RenameFastaEntries(args)
        output_file = tmp_path / "renamed.fa"
        records = iter(
            [
                SeqRecord(Seq("A" * 65), id="a", description="a original description"),
                SeqRecord(Seq("C" * 4), id="b", description="b original description"),
            ]
        )

        renamed_count, total_records = service.replace_ids_and_write(
            str(output_file), records, {"a": "A"}
        )

        assert (renamed_count, total_records) == (1, 2)
        assert output_file.read_text() == (
            ">A\n"
            f"{'A' * 60}\n"
            f"{'A' * 5}\n"
            ">b original description\n"
            "CCCC\n"
        )

    def test_replace_ids_in_file_and_write_preserves_fasta_formatting(self, tmp_path, args):
        service = RenameFastaEntries(args)
        input_file = tmp_path / "input.fa"
        output_file = tmp_path / "renamed.fa"
        input_file.write_text(
            ">a original description\n"
            f"{'A' * 40}\n"
            f"{'A' * 25}\n"
            ">b original description\n"
            "CCCC\n"
        )

        renamed_count, total_records = service.replace_ids_in_file_and_write(
            str(output_file),
            str(input_file),
            {"a": "A"},
        )

        assert (renamed_count, total_records) == (1, 2)
        assert output_file.read_text() == (
            ">A\n"
            f"{'A' * 60}\n"
            f"{'A' * 5}\n"
            ">b original description\n"
            "CCCC\n"
        )

    def test_replace_ids_in_file_and_write_preserves_falsey_mapped_ids(
        self, tmp_path, args
    ):
        service = RenameFastaEntries(args)
        input_file = tmp_path / "input.fa"
        output_file = tmp_path / "renamed.fa"
        input_file.write_text(">a original description\nACGT\n")

        renamed_count, total_records = service.replace_ids_in_file_and_write(
            str(output_file),
            str(input_file),
            {"a": ""},
        )

        assert (renamed_count, total_records) == (1, 1)
        assert output_file.read_text() == ">\nACGT\n"

    def test_replace_ids_in_file_and_write_writes_single_line_sequences_directly(
        self, tmp_path, args
    ):
        service = RenameFastaEntries(args)
        input_file = tmp_path / "input.fa"
        output_file = tmp_path / "renamed.fa"
        input_file.write_text(
            ">a original description\n"
            f"{'A' * 60}\n"
            ">b original description\n"
            "CCCC\n"
        )

        renamed_count, total_records = service.replace_ids_in_file_and_write(
            str(output_file),
            str(input_file),
            {"a": ""},
        )

        assert (renamed_count, total_records) == (1, 2)
        assert output_file.read_text() == (
            ">\n"
            f"{'A' * 60}\n"
            ">b original description\n"
            "CCCC\n"
        )

    def test_write_wrapped_fasta_sequence_helper(self, tmp_path):
        output_file = tmp_path / "wrapped.fa"
        with open(output_file, "w") as handle:
            RenameFastaEntries._write_wrapped_fasta_sequence(handle, "A" * 65)
            RenameFastaEntries._write_wrapped_fasta_sequence(handle, "")

        assert output_file.read_text() == f"{'A' * 60}\n{'A' * 5}\n"

    def test_write_wrapped_fasta_sequence_custom_width(self, tmp_path):
        output_file = tmp_path / "wrapped_custom.fa"
        with open(output_file, "w") as handle:
            RenameFastaEntries._write_wrapped_fasta_sequence(
                handle,
                "ABCDEFG",
                width=3,
            )

        assert output_file.read_text() == "ABC\nDEF\nG\n"

    def test_run_json_payload(self, tmp_path, mocker):
        fasta = tmp_path / "input.fa"
        fasta.write_text(">a\nACGT\n")
        args = Namespace(
            fasta=str(fasta),
            idmap="/some/path/to/idmap.txt",
            output="/tmp/out.fa",
            json=True,
        )
        service = RenameFastaEntries(args)
        mocker.patch.object(RenameFastaEntries, "load_idmap", return_value={"a": "A"})
        mocker.patch.object(RenameFastaEntries, "replace_ids_in_file_and_write", return_value=(1, 2))
        mocked_json = mocker.patch("phykit.services.alignment.rename_fasta_entries.print_json")
        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload["input_fasta"] == str(fasta)
        assert payload["idmap"] == "/some/path/to/idmap.txt"
        assert payload["output_file"] == "/tmp/out.fa.renamed.fa"
        assert payload["total_records"] == 2
        assert payload["renamed_records"] == 1

    def test_run_uses_path_exists_guard_before_rename(self, mocker):
        args = Namespace(
            fasta="/some/path/to/file.fa",
            idmap="/some/path/to/idmap.txt",
            output="/tmp/out.fa",
            json=False,
        )
        service = RenameFastaEntries(args)
        exists = mocker.patch(
            "phykit.services.alignment.rename_fasta_entries._path_exists",
            return_value=True,
        )
        load_idmap = mocker.patch.object(service, "load_idmap", return_value={"a": "A"})
        replace = mocker.patch.object(
            service,
            "replace_ids_in_file_and_write",
            return_value=(1, 2),
        )

        service.run()

        exists.assert_called_once_with("/some/path/to/file.fa")
        load_idmap.assert_called_once_with("/some/path/to/idmap.txt")
        replace.assert_called_once_with(
            "/tmp/out.fa.renamed.fa",
            "/some/path/to/file.fa",
            {"a": "A"},
        )
