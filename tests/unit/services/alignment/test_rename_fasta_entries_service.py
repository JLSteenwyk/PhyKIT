from argparse import Namespace

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.rename_fasta_entries import RenameFastaEntries


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

    def test_run_json_payload(self, mocker):
        args = Namespace(
            fasta="/some/path/to/file.fa",
            idmap="/some/path/to/idmap.txt",
            output="/tmp/out.fa",
            json=True,
        )
        service = RenameFastaEntries(args)
        mocker.patch("phykit.services.alignment.rename_fasta_entries.SeqIO.parse", return_value=iter([]))
        mocker.patch.object(RenameFastaEntries, "load_idmap", return_value={"a": "A"})
        mocker.patch.object(RenameFastaEntries, "replace_ids_and_write", return_value=(1, 2))
        mocked_json = mocker.patch("phykit.services.alignment.rename_fasta_entries.print_json")
        service.run()
        payload = mocked_json.call_args.args[0]
        assert payload["input_fasta"] == "/some/path/to/file.fa"
        assert payload["idmap"] == "/some/path/to/idmap.txt"
        assert payload["output_file"] == "/tmp/out.fa.renamed.fa"
        assert payload["total_records"] == 2
        assert payload["renamed_records"] == 1
