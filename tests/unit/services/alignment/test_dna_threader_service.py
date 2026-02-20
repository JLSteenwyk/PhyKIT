from argparse import Namespace

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phykit.services.alignment.dna_threader import DNAThreader


@pytest.fixture
def args():
    return Namespace(
        stop=True,
        protein="/some/path/to/protein.fa",
        nucleotide="/some/path/to/nucleotide.fa",
        clipkit_log_file=None,
    )


class TestDNAThreader:
    def test_init_sets_expected_attrs(self, args):
        service = DNAThreader(args)
        assert service.remove_stop_codon is True
        assert service.protein_file_path == args.protein
        assert service.nucleotide_file_path == args.nucleotide
        assert service.clipkit_log_file is None
        assert service.json_output is False

    def test_create_mask_without_clipkit_log(self, args):
        service = DNAThreader(args)
        assert service.create_mask(6) == [True] * 6

    def test_normalize_p_seq(self, args):
        service = DNAThreader(args)
        assert service.normalize_p_seq(Seq("A-*")) == "AAA---***"

    def test_normalize_n_seq_handles_gaps(self, args):
        service = DNAThreader(args)
        assert service.normalize_n_seq(Seq("AAACCC"), Seq("A-?")) == "AAA------"

    def test_thread_exits_on_empty_proteins(self, mocker, args):
        service = DNAThreader(args)
        mocked_print = mocker.patch("builtins.print")
        with pytest.raises(SystemExit) as excinfo:
            service.thread(iter([]))
        assert excinfo.value.code == 2
        mocked_print.assert_called_once_with("Protein file is empty or incorrectly formatted.")

    def test_thread_basic(self, mocker, args):
        service = DNAThreader(args)
        prot_records = [SeqRecord(Seq("AC"), id="gene1")]
        nucl_records = [SeqRecord(Seq("AAACCC"), id="gene1")]
        mocker.patch("phykit.services.alignment.dna_threader.SeqIO.parse", return_value=iter(nucl_records))

        pal2nal = service.thread(iter(prot_records))
        assert pal2nal["gene1"] == "AAACCC"

    def test_run_json_output(self, mocker):
        args = Namespace(
            stop=True,
            protein="/some/path/to/protein.fa",
            nucleotide="/some/path/to/nucleotide.fa",
            clipkit_log_file=None,
            json=True,
        )
        service = DNAThreader(args)
        mocker.patch("phykit.services.alignment.dna_threader.SeqIO.parse", return_value=iter([]))
        mocker.patch.object(DNAThreader, "thread", return_value={"gene1": "AAACCC"})
        mocked_json = mocker.patch("phykit.services.alignment.dna_threader.print_json")

        service.run()

        payload = mocked_json.call_args.args[0]
        assert payload["remove_stop_codon"] is True
        assert payload["rows"] == payload["taxa"]
        assert payload["rows"][0] == {"taxon": "gene1", "sequence": "AAACCC"}
